#################################
# Custom functions for Beyond the Border Wars
###################################
require(sf)
require(tidyverse)
require(tidycensus)

#Aggregate blocks to block groups
blocks_to_bg <- function(blk){
  blk %>% st_drop_geometry() %>%
    mutate(GEOID = substr(GEOID, 1, 12)) %>%
    group_by(GEOID) %>%
    summarize(across(starts_with("sum_"), \(x) sum(x, na.rm = TRUE)),
              .groups = "drop")
}

#Target count weighting interpolation
tcw_blocks<-function(sf1,sf2,voteCols=c('biden','trump')){
  #sf1 file with units containing value we want to interpolate from (VTDs)
  #sf2 intermediate unit AND target output unit (blocks)
  sf1_sf2<-st_intersection(sf1,sf2) #intersection of sf1 and sf2 two geographies
  a<-st_area(sf1_sf2) #calculate the area of each polygon
  sf1_sf2<-sf1_sf2[as.numeric(a) > 20,] #only look at units with an area greater than 20 square feet

  #For each observation in sf1 (voting districts) we need the votes by candidate
  sf1_sf2 <- sf1_sf2 %>%
    group_by(GEOID20) %>%
    mutate(vtd_vap = sum(VAP)) %>% #vtd_vap = VAP (from blocks) aggregated to vtd
    ungroup() %>%
    mutate(weight = ifelse(vtd_vap == 0, 0, VAP/vtd_vap)) %>% #VAP weight; 0 where no voting population
    mutate(across(all_of(voteCols), function(x) x*weight)) #replace voteCols with the weighted values
  #Now summarize vote columns back to blocks (sf2). 
  tmp_sf2<-sf1_sf2 %>% st_drop_geometry() %>% group_by(GEOID) %>% summarize(across(all_of(c(voteCols,"TotalPop","VAP","BVAP")), sum, .names = "sum_{.col}"))
  test<-left_join(sf2,tmp_sf2,by="GEOID")
  return(test)
}


# BG-level SMC + shortburst pipeline helpers
list_smc_batches <- function(ndists){
  d_dir <- file.path(bg_smc_dir, sprintf("d%02d", ndists))
  if(!dir.exists(d_dir)) return(character(0))
  fs <- sort(list.files(d_dir, pattern = "^batch_\\d+\\.Rdata$", full.names = TRUE))
  fs[!grepl("\\.tmp$", fs)]
}

#Pull the integer plans matrix out of a redist_plans 
#conflict with built in function threw errors, this is a kludgy workaround
.plans_matrix <- function(plans){
  m <- attr(plans, "plans_matrix")
  if(is.null(m))
    m <- tryCatch(redist::get_plans_matrix(plans),
                  error = function(e) NULL)
  if(is.null(m))
    stop("Cannot extract plans matrix from object of class ",
         paste(class(plans), collapse = "/"))
  m
}

#One contiguous, pop-balanced plan from any available batch — used to
#seed redist_shortburst()
get_one_smc_init_plan <- function(ndists){
  bs <- list_smc_batches(ndists)
  if(!length(bs)) stop("no SMC batches found for d=", ndists)
  e <- new.env(); load(bs[1], envir = e)
  .plans_matrix(e$plans_obj)[, 2]
}

# Adjacency / redist_map helpers
#NOTE on indexing: redist::redist.adjacency() returns 0-indexed neighbor
#IDs but igraph::graph_from_adj_list() expects 1-indexed.
build_pa_map_bg <- function(bgs, ndists, seed_col, pop_tol = 0.01){
  shp_local <- bgs[, c("GEOID","sum_TotalPop","sum_VAP","sum_BVAP",
                       "sum_biden","sum_trump",
                       "county_fips","cousub_id", seed_col), drop = FALSE]
  shp_local$.seed <- shp_local[[seed_col]]
  redist_map(shp_local, total_pop = sum_TotalPop, pop_tol = pop_tol,
             existing_plan = .seed)
}

#SMC ensemble generation

#Counts the number of sims currently on disk.
smc_nsims_accumulated <- function(ndists, batch_size = 10000L){
  length(list_smc_batches(ndists)) * batch_size
}

#Run through a single config (e.g. 6 districts)
#set up to be callable in parallel and used by multiple machines at once
#complicated code to make seeds replicable but not in conflict
ensure_smc_one <- function(ndists, target_nsims = 250000L, batch_size = 10000L,
                            ncores = max(1L, parallel::detectCores() - 1L)){
  current <- smc_nsims_accumulated(ndists, batch_size)
  if(current >= target_nsims) return(invisible())
  if(!file.exists(dataset_cache)) stop("dataset cache missing — re-run stage 1.")
  load(dataset_cache)
  
  #Seed plan via appropriate district base (17,18, or 50)
  bgs$seed <- 1L
  if(ndists == 3){
    bgs$seed[bgs$RoughCD116 > 6]  <- 2L
    bgs$seed[bgs$RoughCD116 > 12] <- 3L
  }else if(ndists == 6){
    bgs$seed[bgs$RoughCD116 > 3]  <- 2L
    bgs$seed[bgs$RoughCD116 > 6]  <- 3L
    bgs$seed[bgs$RoughCD116 > 9]  <- 4L
    bgs$seed[bgs$RoughCD116 > 12] <- 5L
    bgs$seed[bgs$RoughCD116 > 15] <- 6L
  }else if(ndists == 9){
    bgs$seed[bgs$RoughCD116 > 2]  <- 2L
    bgs$seed[bgs$RoughCD116 > 4]  <- 3L
    bgs$seed[bgs$RoughCD116 > 6]  <- 4L
    bgs$seed[bgs$RoughCD116 > 8]  <- 5L
    bgs$seed[bgs$RoughCD116 > 10] <- 6L
    bgs$seed[bgs$RoughCD116 > 12] <- 7L
    bgs$seed[bgs$RoughCD116 > 14] <- 8L
    bgs$seed[bgs$RoughCD116 > 16] <- 9L
  }else if(ndists == 17){
    bgs$seed <- as.integer(bgs$RoughCD118)
  }else if(ndists == 18){
    bgs$seed <- as.integer(bgs$RoughCD116)
  }else if(ndists == 50){
    bgs$seed <- as.integer(bgs$RoughSLDU)
  }else{
    stop("Unsupported ndists: ", ndists)
  }
  bgs$seed[is.na(bgs$seed)] <- 1L
  bgs$seed[bgs$seed < 1 | bgs$seed > ndists] <- 1L

  pa_map <- build_pa_map_bg(bgs, ndists, seed_col = "seed")

  d_dir <- file.path(bg_smc_dir, sprintf("d%02d", ndists))
  dir.create(d_dir, showWarnings = FALSE, recursive = TRUE)

  while(current < target_nsims){
    this_nsims <- as.integer(min(batch_size, target_nsims - current))
    next_idx <- length(list_smc_batches(ndists)) + 1L
    out_f    <- file.path(d_dir, sprintf("batch_%03d.Rdata", next_idx))
    #Per-batch seed: wall time + index so two batches launched in the
    #same R session (or same second) cannot share a seed.
    batch_seed <- (as.integer(Sys.time()) + next_idx * 7919L) %% .Machine$integer.max

    message(sprintf("d=%d batch %03d: redist_smc nsims=%d ncores=%d  (have %d/%d, seed=%d) ...",
                    ndists, next_idx, this_nsims, ncores, current, target_nsims, batch_seed))
    t0 <- Sys.time()
    set.seed(batch_seed)
    plans_obj <- redist_smc(pa_map, nsims = this_nsims, compactness = 1,
                            ncores = ncores)
    wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  batch %03d done in %.1f min", next_idx, wall_min))

    tmp_f <- paste0(out_f, ".tmp")
    save(plans_obj, batch_seed, this_nsims, wall_min,
         file = tmp_f, compress = "xz")
    file.rename(tmp_f, out_f)
    message("  wrote ", out_f)
    rm(plans_obj); gc(verbose = FALSE)
    current <- current + this_nsims
  }
}

#Loops ensure_smc_one across every config in district_specs.
ensure_smc_all <- function(target_nsims = 250000L, batch_size = 10000L){
  for(spec in district_specs)
    ensure_smc_one(spec$ndists, target_nsims, batch_size)
  invisible()
}

#Admin-preserving SMC variant. Same SMC sampler as ensure_smc_one, but adds
#a redist::redist_constr() that penalizes plans for splitting administrative
#units (counties or municipalities).
ensure_smc_one_admin <- function(ndists, admin_col = "county_fips",
                                  target_nsims = 10000L, batch_size = 10000L,
                                  strength = 5, pop_temper = 0,
                                  ncores = max(1L, parallel::detectCores() - 1L)){
  d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, admin_col))
  dir.create(d_dir, showWarnings = FALSE, recursive = TRUE)

  list_batches <- function(){
    fs <- sort(list.files(d_dir, pattern = "^batch_\\d+\\.Rdata$", full.names = TRUE))
    fs[!grepl("\\.tmp$", fs)]
  }
  current <- length(list_batches()) * batch_size
  if(current >= target_nsims) return(invisible())
  if(!file.exists(dataset_cache)) stop("dataset cache missing — re-run stage 1.")
  load(dataset_cache)

  #Seed plan via the same RoughCD116/RoughCD118/RoughSLDU logic as ensure_smc_one.
  bgs$seed <- 1L
  if(ndists == 3){
    bgs$seed[bgs$RoughCD116 > 6]  <- 2L
    bgs$seed[bgs$RoughCD116 > 12] <- 3L
  }else if(ndists == 6){
    bgs$seed[bgs$RoughCD116 > 3]  <- 2L
    bgs$seed[bgs$RoughCD116 > 6]  <- 3L
    bgs$seed[bgs$RoughCD116 > 9]  <- 4L
    bgs$seed[bgs$RoughCD116 > 12] <- 5L
    bgs$seed[bgs$RoughCD116 > 15] <- 6L
  }else if(ndists == 9){
    bgs$seed[bgs$RoughCD116 > 2]  <- 2L
    bgs$seed[bgs$RoughCD116 > 4]  <- 3L
    bgs$seed[bgs$RoughCD116 > 6]  <- 4L
    bgs$seed[bgs$RoughCD116 > 8]  <- 5L
    bgs$seed[bgs$RoughCD116 > 10] <- 6L
    bgs$seed[bgs$RoughCD116 > 12] <- 7L
    bgs$seed[bgs$RoughCD116 > 14] <- 8L
    bgs$seed[bgs$RoughCD116 > 16] <- 9L
  }else if(ndists == 17){
    bgs$seed <- as.integer(bgs$RoughCD118)
  }else if(ndists == 18){
    bgs$seed <- as.integer(bgs$RoughCD116)
  }else if(ndists == 50){
    bgs$seed <- as.integer(bgs$RoughSLDU)
  }else stop("Unsupported ndists: ", ndists)
  bgs$seed[is.na(bgs$seed)] <- 1L
  bgs$seed[bgs$seed < 1 | bgs$seed > ndists] <- 1L

  pa_map <- build_pa_map_bg(bgs, ndists, seed_col = "seed")

  #Build the constraint. 
  constr <- redist::redist_constr(pa_map)
  if(admin_col == "county_fips"){
    constr <- redist::add_constr_splits(constr, strength = strength, admin = county_fips)
  }else if(admin_col == "cousub_id"){
    constr <- redist::add_constr_splits(constr, strength = strength, admin = cousub_id)
  }else{
    stop("admin_col must be 'county_fips' or 'cousub_id'")
  }

  while(current < target_nsims){
    this_nsims <- as.integer(min(batch_size, target_nsims - current))
    next_idx <- length(list_batches()) + 1L
    out_f    <- file.path(d_dir, sprintf("batch_%03d.Rdata", next_idx))
    batch_seed <- (as.integer(Sys.time()) + next_idx * 7919L) %% .Machine$integer.max

    message(sprintf("d=%d (%s, str=%.1f) batch %03d: redist_smc nsims=%d ncores=%d  (have %d/%d, seed=%d) ...",
                    ndists, admin_col, strength, next_idx, this_nsims, ncores,
                    current, target_nsims, batch_seed))
    t0 <- Sys.time()
    set.seed(batch_seed)
    plans_obj <- redist::redist_smc(pa_map, nsims = this_nsims, compactness = 1,
                                     constraints = constr, ncores = ncores,
                                     pop_temper = pop_temper)
    wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  batch %03d done in %.1f min", next_idx, wall_min))

    tmp_f <- paste0(out_f, ".tmp")
    save(plans_obj, batch_seed, this_nsims, wall_min, strength,
         file = tmp_f, compress = "xz")
    file.rename(tmp_f, out_f)
    rm(plans_obj); gc(verbose = FALSE)
    current <- current + this_nsims
  }
}

#Loops ensure_smc_one_admin across every config in district_specs.
ensure_smc_all_admin <- function(admin_col = "county_fips",
                                  target_nsims = 10000L, batch_size = 10000L,
                                  strength = 5, pop_temper = 0){
  resolve <- function(param, nd){
    if(is.null(names(param))) return(param)
    key <- as.character(nd)
    if(key %in% names(param))       return(unname(param[key]))
    if("default" %in% names(param)) return(unname(param["default"]))
    stop("per-config parameter has no 'default' entry and no override for ndists=", nd)
  }
  for(spec in district_specs)
    ensure_smc_one_admin(spec$ndists, admin_col, target_nsims, batch_size,
                          resolve(strength,   spec$ndists),
                          resolve(pop_temper, spec$ndists))
  invisible()
}

#Combined county+cousub-preserving SMC variant. built in relaxation of pop temper to account for failures in plans with more districts (esp. 50 district) 
ensure_smc_one_combined <- function(ndists,
                                     target_nsims    = 250000L,
                                     batch_size      = 10000L,
                                     county_strength = 3,
                                     cousub_strength = 5,
                                     pop_tol         = 0.01,
                                     pop_temper      = 0,
                                     fallback_ladder = list(
                                       list(pop_temper = 0.01),
                                       list(pop_temper = 0.02),
                                       list(pop_temper = 0.03),
                                       list(pop_temper = 0.03, county_strength = 3, cousub_strength = 3)
                                     ),
                                     dir_suffix = "county_cousub",
                                     ncores = max(1L, parallel::detectCores() - 1L)){
  #dir_suffix isolates variant ensembles (e.g. a higher-penalty sweep) in their
  #own d{NN}_<suffix>/ directory without touching the main county_cousub run.
  #The constraints below are always county_fips + cousub_id; only the strengths
  #and the output directory change.
  d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, dir_suffix))
  dir.create(d_dir, showWarnings = FALSE, recursive = TRUE)

  worker_pid <- Sys.getpid()

  list_batches <- function(){
    fs <- sort(list.files(d_dir, pattern = "^batch_.*\\.Rdata$", full.names = TRUE))
    fs[!grepl("\\.tmp$", fs)]
  }
  current <- length(list_batches()) * batch_size
  if(current >= target_nsims) return(invisible())
  if(!file.exists(dataset_cache)) stop("dataset cache missing — re-run stage 1.")
  load(dataset_cache)

  #Seed plan via the same RoughCD116/RoughCD118/RoughSLDU logic as ensure_smc_one.
  bgs$seed <- 1L
  if(ndists == 3){
    bgs$seed[bgs$RoughCD116 > 6]  <- 2L
    bgs$seed[bgs$RoughCD116 > 12] <- 3L
  }else if(ndists == 6){
    bgs$seed[bgs$RoughCD116 > 3]  <- 2L
    bgs$seed[bgs$RoughCD116 > 6]  <- 3L
    bgs$seed[bgs$RoughCD116 > 9]  <- 4L
    bgs$seed[bgs$RoughCD116 > 12] <- 5L
    bgs$seed[bgs$RoughCD116 > 15] <- 6L
  }else if(ndists == 9){
    bgs$seed[bgs$RoughCD116 > 2]  <- 2L
    bgs$seed[bgs$RoughCD116 > 4]  <- 3L
    bgs$seed[bgs$RoughCD116 > 6]  <- 4L
    bgs$seed[bgs$RoughCD116 > 8]  <- 5L
    bgs$seed[bgs$RoughCD116 > 10] <- 6L
    bgs$seed[bgs$RoughCD116 > 12] <- 7L
    bgs$seed[bgs$RoughCD116 > 14] <- 8L
    bgs$seed[bgs$RoughCD116 > 16] <- 9L
  }else if(ndists == 17){
    bgs$seed <- as.integer(bgs$RoughCD118)
  }else if(ndists == 18){
    bgs$seed <- as.integer(bgs$RoughCD116)
  }else if(ndists == 50){
    bgs$seed <- as.integer(bgs$RoughSLDU)
  }else stop("Unsupported ndists: ", ndists)
  bgs$seed[is.na(bgs$seed)] <- 1L
  bgs$seed[bgs$seed < 1 | bgs$seed > ndists] <- 1L

  #pa_map only depends on pop_tol, which is fixed here, so we build it once.
  pa_map <- build_pa_map_bg(bgs, ndists, seed_col = "seed", pop_tol = pop_tol)

  #Full attempt ladder
  initial_attempt <- list(pop_temper      = pop_temper,
                          county_strength = county_strength,
                          cousub_strength = cousub_strength)
  fill_defaults <- function(step){
    for(nm in c("pop_temper","county_strength","cousub_strength")){
      if(is.null(step[[nm]])) step[[nm]] <- initial_attempt[[nm]]
    }
    step
  }
  attempts <- c(list(initial_attempt), lapply(fallback_ladder, fill_defaults))

  #Infer the floor from existing batches so we don't redo work we know fails
  match_step <- function(used){
    for(k in seq_along(attempts)){
      a <- attempts[[k]]
      if(isTRUE(all.equal(a$pop_temper,      used$pop_temper))      &&
         isTRUE(all.equal(a$county_strength, used$county_strength)) &&
         isTRUE(all.equal(a$cousub_strength, used$cousub_strength))) return(k)
    }
    return(1L)
  }
  attempt_floor <- 1L
  for(bf in list_batches()){
    e <- new.env(); tryCatch(load(bf, envir = e), error = function(err) NULL)
    used <- list(pop_temper      = if(exists("used_pop_temper",      envir = e)) e$used_pop_temper      else pop_temper,
                 county_strength = if(exists("used_county_strength", envir = e)) e$used_county_strength else county_strength,
                 cousub_strength = if(exists("used_cousub_strength", envir = e)) e$used_cousub_strength else cousub_strength)
    k <- match_step(used)
    if(k > attempt_floor) attempt_floor <- k
  }

  worker_idx <- 0L   # worker-local counter; increments per successful batch
  while(current < target_nsims){
    worker_idx <- worker_idx + 1L
    this_nsims <- as.integer(min(batch_size, target_nsims - current))
    out_f      <- file.path(d_dir, sprintf("batch_pid%d_%03d.Rdata",
                                            worker_pid, worker_idx))
    batch_seed <- as.integer(
      (as.numeric(Sys.time()) +
         worker_idx * 7919 +
         worker_pid * 104729) %% .Machine$integer.max
    )

    plans_obj <- NULL
    used_pop_temper      <- NA_real_
    used_county_strength <- NA_real_
    used_cousub_strength <- NA_real_
    wall_min             <- NA_real_

    for(k in seq(attempt_floor, length(attempts))){
      a <- attempts[[k]]

      constr <- redist::redist_constr(pa_map)
      constr <- redist::add_constr_splits(constr, strength = a$county_strength,
                                          admin = county_fips)
      constr <- redist::add_constr_splits(constr, strength = a$cousub_strength,
                                          admin = cousub_id)

      message(sprintf("d=%d combined pid=%d batch %03d step %d/%d (pop_temper=%g cnty_str=%g cousub_str=%g pop_tol=%g): redist_smc nsims=%d ncores=%d  (have %d/%d, seed=%d) ...",
                      ndists, worker_pid, worker_idx, k, length(attempts),
                      a$pop_temper, a$county_strength, a$cousub_strength, pop_tol,
                      this_nsims, ncores, current, target_nsims, batch_seed))
      t0 <- Sys.time()
      set.seed(batch_seed)
      result <- tryCatch(
        redist::redist_smc(pa_map, nsims = this_nsims, compactness = 1,
                           constraints = constr, ncores = ncores,
                           pop_temper = a$pop_temper),
        error = function(e){
          message("    SMC failed at step ", k, ": ", conditionMessage(e))
          NULL
        })

      if(!is.null(result)){
        plans_obj            <- result
        used_pop_temper      <- a$pop_temper
        used_county_strength <- a$county_strength
        used_cousub_strength <- a$cousub_strength
        wall_min             <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        message(sprintf("  pid=%d batch %03d done in %.1f min (step %d/%d)",
                        worker_pid, worker_idx, wall_min, k, length(attempts)))
        #Future batches in this run start here.
        attempt_floor <- k
        break
      }
    }

    if(is.null(plans_obj)){
      stop(sprintf("ensure_smc_one_combined: exhausted fallback ladder (%d steps) for d=%d pid=%d batch %03d",
                   length(attempts), ndists, worker_pid, worker_idx))
    }

    tmp_f <- paste0(out_f, ".tmp")
    save(plans_obj, batch_seed, this_nsims, wall_min,
         used_pop_temper, used_county_strength, used_cousub_strength,
         pop_tol, worker_pid, worker_idx,
         file = tmp_f, compress = "xz")
    file.rename(tmp_f, out_f)
    message("  wrote ", out_f,
            " (pop_temper=", used_pop_temper,
            " cnty_str=",    used_county_strength,
            " cousub_str=",  used_cousub_strength, ")")
    rm(plans_obj); gc(verbose = FALSE)
    #Re-read directory in case other workers contributed since last check.
    current <- length(list_batches()) * batch_size
  }
  invisible()
}

#Heler for function below...
resolve <- function(param, nd){
  if(is.list(param) && !is.null(names(param)) && all(sapply(param, is.list))){
    key <- as.character(nd)
    if(key %in% names(param))       return(param[[key]])
    if("default" %in% names(param)) return(param[["default"]])
    stop("per-config list parameter has no 'default' entry and no override for ndists=", nd)
  }
  if(is.atomic(param)){
    if(is.null(names(param))) return(param)
    key <- as.character(nd)
    if(key %in% names(param))       return(unname(param[key]))
    if("default" %in% names(param)) return(unname(param["default"]))
    stop("per-config parameter has no 'default' entry and no override for ndists=", nd)
  }
  param   # raw list (e.g. the default fallback_ladder)
}
#Loops ensure_smc_one_combined across every config in district_specs.
ensure_smc_all_combined <- function(target_nsims    = 250000L,
                                    batch_size      = 10000L,
                                    county_strength = 3,
                                    cousub_strength = 5,
                                    pop_tol         = 0.01,
                                    pop_temper      = 0,
                                    fallback_ladder = list(
                                      list(pop_temper = 0.01),
                                      list(pop_temper = 0.02),
                                      list(pop_temper = 0.03),
                                      list(pop_temper = 0.03, county_strength = 3, cousub_strength = 3)
                                    )){

  for(spec in district_specs)
    ensure_smc_one_combined(spec$ndists,
                            target_nsims    = target_nsims,
                            batch_size      = batch_size,
                            county_strength = resolve(county_strength, spec$ndists),
                            cousub_strength = resolve(cousub_strength, spec$ndists),
                            pop_tol         = resolve(pop_tol,         spec$ndists),
                            pop_temper      = resolve(pop_temper,      spec$ndists),
                            fallback_ladder = resolve(fallback_ladder, spec$ndists))
  invisible()
}

#---- Unconstrained tract-resolution ensemble (Supplementary Materials Sec. 2) ----
#The tract-vs-block-group resolution contrast needs an SMC ensemble built from
#census tracts at pop_tol = 0.01 with NO administrative-preservation constraints.
#These helpers (lifted from the former standalone Run_Tract1Pct.R) let
#BBW - Data Prep BG.R regenerate any missing tract plans the same cache-gated way
#it builds the block-group ensembles. Plans land in
#./Intermediate Data/tract_1pct/smc_plans/pa_plans_tract1pct_d{NN}.Rdata, which
#the tract-metrics block in Data Prep then rescores.

#Build the tract shapefile with the multi-member seed columns (CD3/CD6/CD9 carve
#the 18 RoughCD116 districts into 3/6/9-district nesting; CD/RoughCD116/SLDU are
#the single-member seeds).
build_tract_shp <- function(){
  load("./Intermediate Data/tracts.interpolated.Rdata")  # tracts.interpolated
  shp <- tracts.interpolated[, c("GEOID","sum_TotalPop","sum_biden","sum_trump",
                                 "POP18Black","POP18Plus_PL",
                                 "RoughCD118","RoughSLDU","RoughCD116")]
  shp$FIPS       <- as.numeric(shp$GEOID)
  shp$CD         <- as.numeric(shp$RoughCD118)
  shp$SLDU       <- as.numeric(shp$RoughSLDU)
  shp$RoughCD116 <- as.numeric(shp$RoughCD116)
  shp$two_party  <- shp$sum_biden + shp$sum_trump
  shp$CD3 <- 1L
  shp[shp$RoughCD116 > 6,  "CD3"] <- 2L
  shp[shp$RoughCD116 > 12, "CD3"] <- 3L
  shp$CD6 <- 1L
  shp[shp$RoughCD116 > 3,  "CD6"] <- 2L
  shp[shp$RoughCD116 > 6,  "CD6"] <- 3L
  shp[shp$RoughCD116 > 9,  "CD6"] <- 4L
  shp[shp$RoughCD116 > 12, "CD6"] <- 5L
  shp[shp$RoughCD116 > 15, "CD6"] <- 6L
  shp$CD9 <- 1L
  shp[shp$RoughCD116 > 2,  "CD9"] <- 2L
  shp[shp$RoughCD116 > 4,  "CD9"] <- 3L
  shp[shp$RoughCD116 > 6,  "CD9"] <- 4L
  shp[shp$RoughCD116 > 8,  "CD9"] <- 5L
  shp[shp$RoughCD116 > 10, "CD9"] <- 6L
  shp[shp$RoughCD116 > 12, "CD9"] <- 7L
  shp[shp$RoughCD116 > 14, "CD9"] <- 8L
  shp[shp$RoughCD116 > 16, "CD9"] <- 9L
  shp
}

#Per-config seed columns for the tract ensemble. Kept as its own list (not the
#block-group `district_specs`) so the two never collide in a shared session.
tract_district_specs <- function(){
  specs <- list(
    list(ndists =  3, seed_col = "CD3"),
    list(ndists =  6, seed_col = "CD6"),
    list(ndists =  9, seed_col = "CD9"),
    list(ndists = 17, seed_col = "CD"),
    list(ndists = 18, seed_col = "RoughCD116"),
    list(ndists = 50, seed_col = "SLDU")
  )
  names(specs) <- vapply(specs, function(s) sprintf("d%02d", s$ndists), character(1))
  specs
}

#Generate (cache-gated) the tract SMC plans for one configuration. No-ops if the
#plan cache already exists; otherwise builds the map at pop_tol = 0.01 with no
#admin constraints and runs redist_smc, writing atomically (.tmp then rename).
ensure_tract_smc_one <- function(ndists, nsims = 50000L,
                                 ncores  = max(1L, parallel::detectCores() - 1L),
                                 smc_dir = "./Intermediate Data/tract_1pct/smc_plans"){
  dir.create(smc_dir, showWarnings = FALSE, recursive = TRUE)
  cache_f <- file.path(smc_dir, sprintf("pa_plans_tract1pct_d%02d.Rdata", ndists))
  if(file.exists(cache_f)){ message("tract smc cache exists: ", cache_f); return(invisible()) }
  spec <- tract_district_specs()[[sprintf("d%02d", ndists)]]
  if(is.null(spec)) stop("unknown tract ndists: ", ndists)

  shp <- build_tract_shp()
  shp$.seed <- shp[[spec$seed_col]]
  pa_map <- redist::redist_map(shp, total_pop = sum_TotalPop, pop_tol = 0.01,
                               existing_plan = .seed)
  pa_map <- dplyr::filter(pa_map, sum_TotalPop > 0)

  message(sprintf("running tract redist_smc d=%d pop_tol=0.01 nsims=%d ncores=%d ...",
                  ndists, nsims, ncores))
  t0 <- Sys.time()
  pa_plans <- redist::redist_smc(pa_map, nsims = nsims, compactness = 1, ncores = ncores)
  message(sprintf("tract smc done in %.1f min",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  tmp_f <- paste0(cache_f, ".tmp")
  save(pa_plans, file = tmp_f, compress = "xz")
  file.rename(tmp_f, cache_f)
  message("wrote ", cache_f)
  invisible()
}

#Loops ensure_tract_smc_one across every tract configuration.
ensure_tract_smc_all <- function(nsims = 50000L, ...){
  for(spec in tract_district_specs()) ensure_tract_smc_one(spec$ndists, nsims = nsims, ...)
  invisible()
}

#Shortburst metric generation
flip_prob_q <- function(delta, q, sigma){
  pnorm(-delta / sigma) + pnorm(-(q - delta) / sigma)
}

#Helper: pull the integer plans matrix out of a redist_shortburst result.
.shortburst_plan_matrix <- function(sb){
  m <- attr(sb, "plans_matrix")
  if(is.null(m)) m <- tryCatch(redist::get_plans_matrix(sb),
                                error = function(e) NULL)
  if(is.null(m)) stop("could not extract plans matrix from shortburst result")
  m
}

#Helper: for each column of `pm`, compute the k-th highest district-level
#share of group_col / total_col
.score_plans_topk <- function(pm, bgs_dt, group_col, total_col, k){
  g <- bgs_dt[[group_col]]; t <- bgs_dt[[total_col]]
  vapply(seq_len(ncol(pm)), function(j){
    plan <- pm[, j]
    gd <- as.numeric(rowsum(g, plan))
    td <- as.numeric(rowsum(t, plan))
    pcts <- gd / td
    sort(pcts, decreasing = TRUE)[k]
  }, numeric(1))
}

#Seats-counting scorer. Counts total seats won by `group_col` across the
#plan
make_seats_scorer <- function(map, group_col, total_col, num_members, m_offset){
  group_vec <- as.numeric(map[[group_col]])
  total_vec <- as.numeric(map[[total_col]])
  ndists    <- as.integer(attr(map, "ndists"))
  factor    <- num_members + m_offset
  group_pct_c <- get("group_pct", envir = asNamespace("redist"))
  function(plans){
    pm <- if(inherits(plans, "redist_plans")) redist::get_plans_matrix(plans)
          else if(is.matrix(plans))           plans
          else                                as.matrix(plans)
    pct_mat <- group_pct_c(pm, group_vec, total_vec, ndists)
    colSums(floor(pct_mat * factor))
  }
}

#Backward-compatible wrapper for the original Black-seats-only entry point.
make_black_seats_scorer <- function(map, num_members, m_offset){
  make_seats_scorer(map, "sum_BVAP", "sum_VAP", num_members, m_offset)
}

#Competitiveness scorer
make_competitiveness_scorer <- function(map, num_members, metric){
  biden_vec <- as.numeric(map$sum_biden)
  twop_vec  <- as.numeric(map$two_party)
  q         <- 1 / (num_members + 1)
  function(plans){
    pm <- if(inherits(plans, "redist_plans")) redist::get_plans_matrix(plans)
          else if(is.matrix(plans))           plans
          else                                as.matrix(plans)
    nplans <- ncol(pm)
    out <- numeric(nplans)
    for(j in seq_len(nplans)){
      plan <- pm[, j]
      bd <- as.numeric(rowsum(biden_vec, plan))
      td <- as.numeric(rowsum(twop_vec,  plan))
      dp <- bd / td
      v_modq <- dp %% q
      delta  <- pmin(v_modq, q - v_modq)
      out[j] <- if(identical(metric, "share_marginal")) mean(delta < 0.05)
               else mean(flip_prob_q(delta, q, 0.03))
    }
    out
  }
}

#Re-scorer matching make_competitiveness_scorer so the shortburst visit
#log can be re-ranked independently of redist's internal ordering.
.score_plans_competitiveness <- function(pm, bgs_dt, num_members, metric){
  biden <- as.numeric(bgs_dt$sum_biden); biden[is.na(biden)] <- 0
  twop  <- as.numeric(bgs_dt$two_party); twop[is.na(twop)]   <- 0
  q     <- 1 / (num_members + 1)
  vapply(seq_len(ncol(pm)), function(j){
    plan <- pm[, j]
    bd <- as.numeric(rowsum(biden, plan))
    td <- as.numeric(rowsum(twop,  plan))
    dp <- bd / td
    v_modq <- dp %% q
    delta  <- pmin(v_modq, q - v_modq)
    if(identical(metric, "share_marginal")) mean(delta < 0.05)
    else mean(flip_prob_q(delta, q, 0.03))
  }, numeric(1))
}

#Re-scorer matching make_seats_scorer's objective. Recomputes total
#group_col seats for every plan in the shortburst visit log so we can
#pick the actual best independent of redist's internal ordering.
.score_plans_seats <- function(pm, bgs_dt, group_col, total_col, num_members, m_offset, ndists){
  group_vec <- as.numeric(bgs_dt[[group_col]])
  total_vec <- as.numeric(bgs_dt[[total_col]])
  group_pct_c <- get("group_pct", envir = asNamespace("redist"))
  pct_mat <- group_pct_c(pm, group_vec, total_vec, as.integer(ndists))
  colSums(floor(pct_mat * (num_members + m_offset)))
}

#Backward-compatible wrapper for the original Black-seats-only entry point.
.score_plans_black_seats <- function(pm, bgs_dt, num_members, m_offset, ndists){
  .score_plans_seats(pm, bgs_dt, "sum_BVAP", "sum_VAP", num_members, m_offset, ndists)
}

#Seed-picker: scan every plan in the SMC ensemble for this (ndists,
#admin_col) and return, per shortburst target, the single plan from the
#ensemble that maximizes (or minimizes) the target's score. Used as
#`init_plan` for the shortburst search so it starts at a known good
#local region rather than a random sampled plan. 
find_best_seeds <- function(smc_files, bgs_dt, ndists, num_members){
  bvap  <- as.numeric(bgs_dt$sum_BVAP);  bvap[is.na(bvap)]   <- 0
  vap   <- as.numeric(bgs_dt$sum_VAP);   vap[is.na(vap)]     <- 0
  biden <- as.numeric(bgs_dt$sum_biden); biden[is.na(biden)] <- 0
  twop  <- as.numeric(bgs_dt$two_party); twop[is.na(twop)]   <- 0
  group_pct_c <- get("group_pct", envir = asNamespace("redist"))

  empty <- function() list(score = NA_real_, plan = NULL)
  best <- list(
    black_seats_max     = empty(),
    black_seats37_max   = empty(),
    dem_seats_max       = empty(),
    rep_seats_max       = empty(),
    marg5pt_max         = empty(),
    flipprob_max        = empty(),
    dem_max = replicate(ndists, empty(), simplify = FALSE),
    dem_min = replicate(ndists, empty(), simplify = FALSE)
  )
  q_drp <- 1 / (num_members + 1)

  upd_scalar <- function(slot, scores, pm){
    bi   <- which.max(scores)
    cand <- unname(scores[bi])
    if(is.na(best[[slot]]$score) || cand > best[[slot]]$score){
      best[[slot]] <<- list(score = cand, plan = as.integer(pm[, bi]))
    }
  }

  for(f in smc_files){
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
    if(!ok || !exists("plans_obj", envir = e)) next
    pm <- redist::get_plans_matrix(e$plans_obj)
    pm <- pm[, -1, drop = FALSE]   # drop SMC seed plan column
    if(ncol(pm) == 0L){ rm(e); next }

    pct_bvap <- group_pct_c(pm, bvap,  vap,  as.integer(ndists))
    pct_dem  <- group_pct_c(pm, biden, twop, as.integer(ndists))

    upd_scalar("black_seats_max",   colSums(floor(pct_bvap * (num_members + 1L))),    pm)
    upd_scalar("black_seats37_max", colSums(floor(pct_bvap * (num_members + 2L))),    pm)
    upd_scalar("dem_seats_max",     colSums(floor(pct_dem  * (num_members + 1L))),    pm)
    upd_scalar("rep_seats_max",     colSums(floor((1 - pct_dem) * (num_members + 1L))), pm)

    #Competitiveness targets: distance from each district's Dem share to
    #the nearest Droop boundary, summarized per plan as share-of-districts
    #within 5pp (marg5pt_max) and mean per-district flip probability under
    #a normal swing of sigma=0.03 (flipprob_max). Matches the algebra in
    #generate_metrics_bg lines 1408-1413.
    v_modq <- pct_dem %% q_drp
    delta  <- pmin(v_modq, q_drp - v_modq)
    upd_scalar("marg5pt_max",  colMeans(delta < 0.05),                  pm)
    upd_scalar("flipprob_max", colMeans(flip_prob_q(delta, q_drp, 0.03)), pm)

    #k-sweep targets: sort each plan's district shares descending, then
    #for each k pick best (max for dem_max, min for dem_min).
    sorted_dem <- apply(pct_dem, 2, sort, decreasing = TRUE)
    if(!is.matrix(sorted_dem)) sorted_dem <- matrix(sorted_dem, nrow = ndists)
    for(k in seq_len(ndists)){
      kth <- sorted_dem[k, ]
      #dem_max k: maximize the k-th-largest share
      bmax <- which.max(kth); cmax <- unname(kth[bmax])
      if(is.na(best$dem_max[[k]]$score) || cmax > best$dem_max[[k]]$score){
        best$dem_max[[k]] <- list(score = cmax, plan = as.integer(pm[, bmax]))
      }
      #dem_min k: minimize the k-th-largest share
      bmin <- which.min(kth); cmin <- unname(kth[bmin])
      if(is.na(best$dem_min[[k]]$score) || cmin < best$dem_min[[k]]$score){
        best$dem_min[[k]] <- list(score = cmin, plan = as.integer(pm[, bmin]))
      }
    }
    rm(e); gc(verbose = FALSE)
  }
  best
}

ensure_shortburst_one <- function(ndists,
                                   admin_col  = NULL,
                                   max_bursts = 2000L,
                                   burst_size = 10L){
  suffix  <- if(is.null(admin_col)) "" else paste0("_", admin_col)
  cache_f <- file.path(bg_sb_dir, sprintf("shortburst_bg_d%02d%s.Rdata",
                                          ndists, suffix))
  out <- list()
  if(file.exists(cache_f)){
    e <- new.env()
    load_ok <- tryCatch({ load(cache_f, envir = e); TRUE }, error = function(err) FALSE)
    if(load_ok && exists("out", envir = e) && is.list(e$out)){
      out <- e$out
      message(sprintf("d=%02d%s loaded existing cache (%d entries)",
                      ndists, suffix, length(out)))
    }
  }

  #All expected target keys. If everything is already present we exit
  expected_scalar <- c("black_seats_max", "black_seats37_max",
                       "dem_seats_max",   "rep_seats_max",
                       "marg5pt_max",     "flipprob_max")
  expected_kkeys  <- as.vector(outer(c("dem_max", "dem_min"),
                                      seq_len(ndists),
                                      function(t, k) paste0(t, "_k", k)))
  expected_all    <- c(expected_scalar, expected_kkeys)
  if(all(expected_all %in% names(out))){
    message(sprintf("d=%02d%s already complete (%d entries); skipping",
                    ndists, suffix, length(out)))
    return(invisible())
  }

  smc_dir <- if(is.null(admin_col))
    file.path(bg_smc_dir, sprintf("d%02d", ndists))
  else
    file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, admin_col))
  smc_files <- sort(list.files(smc_dir, pattern = "^batch_.*\\.Rdata$",
                                full.names = TRUE))
  smc_files <- smc_files[!grepl("\\.tmp$", smc_files)]
  if(length(smc_files) == 0L)
    stop("no SMC plans in ", smc_dir, "; SMC must run first")

  #Seed with column 2 (= first sampled plan; column 1 is the SMC seed plan).
  e <- new.env(); load(smc_files[1], envir = e)
  init_plan <- as.integer(.plans_matrix(e$plans_obj)[, 2])
  rm(e); gc(verbose = FALSE)

  load(dataset_cache)  # bgs
  stopifnot(nrow(bgs) == length(init_plan))
  bgs$init_plan <- init_plan
  pa_map <- redist_map(
    bgs[, c("GEOID","sum_TotalPop","sum_VAP","sum_BVAP",
             "sum_biden","sum_trump",
             "county_fips","cousub_id","init_plan")],
    total_pop     = sum_TotalPop,
    pop_tol       = 0.01,
    existing_plan = init_plan)
  pa_map$two_party <- pa_map$sum_biden + pa_map$sum_trump

  #Admin constraints — match the ensemble being optimized against so the
  #shortburst envelope is drawn from the SAME constrained map space as the SMC
  #ensemble whose boxplots it overlays.
  #  county_cousub (OLD main): SOFT penalty on BOTH levels (county 1, cousub 3).
  #  ccstruct_m5  (NEW main):  HARD county constraint via redist_shortburst's
  #                            counties= arg (caps county splits at ndists-1,
  #                            the same lever as redist_smc(counties=)) + a SOFT
  #                            cousub penalty at strength 5.
  constraints <- redist::redist_constr(pa_map)
  sb_counties <- NULL   # vector for redist_shortburst(counties=); NULL = no hard county constraint
  if(identical(admin_col, "county_cousub")){
    constraints <- redist::add_constr_splits(constraints, strength = 1,
                                              admin = county_fips)
    constraints <- redist::add_constr_splits(constraints, strength = 3,
                                              admin = cousub_id)
  }else if(identical(admin_col, "ccstruct_m5")){
    sb_counties <- pa_map$county_fips
    constraints <- redist::add_constr_splits(constraints, strength = 5,
                                              admin = cousub_id)
  }

  bgs_dt <- sf::st_drop_geometry(bgs)
  bgs_dt$two_party <- bgs_dt$sum_biden + bgs_dt$sum_trump

  num_members <- 18 / ndists; if(num_members < 2) num_members <- 1

  #Per-target seed selection: scan the SMC ensemble once and identify
  #the best ensemble plan for every shortburst target.
  message(sprintf("d=%02d%s scanning %d batches for per-target seeds ...",
                  ndists, suffix, length(smc_files)))
  seeds <- find_best_seeds(smc_files, bgs_dt, ndists, num_members)
  fallback_seed <- init_plan
  pick_seed <- function(slot, k = NULL){
    s <- if(is.null(k)) seeds[[slot]] else seeds[[slot]][[k]]
    if(is.null(s$plan)) fallback_seed else s$plan
  }

  save_now <- function(){
    tmp_f <- paste0(cache_f, ".tmp")
    save(out, file = tmp_f, compress = "xz")
    file.rename(tmp_f, cache_f)
  }

  one_sweep_topk <- function(sf_score, group_col, total_col, k, maximize, target, init_plan_use){
    sb <- redist_shortburst(pa_map,
                             counties    = sb_counties,
                             score_fn    = sf_score,
                             maximize    = maximize,
                             max_bursts  = max_bursts,
                             burst_size  = burst_size,
                             return_all  = TRUE,
                             init_plan   = init_plan_use,
                             constraints = constraints,
                             verbose     = FALSE)
    pm     <- .shortburst_plan_matrix(sb)
    scores <- .score_plans_topk(pm, bgs_dt, group_col, total_col, k)
    best_j <- if(maximize) which.max(scores) else which.min(scores)
    list(plan       = as.integer(pm[, best_j]),
         score      = unname(scores[best_j]),
         seed_score = unname(scores[1]),
         n_plans    = ncol(pm),
         maximize   = maximize,
         target     = target,
         k          = k)
  }

  one_sweep_seats <- function(sf_score, group_col, total_col, m_offset, target, init_plan_use){
    sb <- redist_shortburst(pa_map,
                             counties    = sb_counties,
                             score_fn    = sf_score,
                             maximize    = TRUE,
                             max_bursts  = max_bursts,
                             burst_size  = burst_size,
                             return_all  = TRUE,
                             init_plan   = init_plan_use,
                             constraints = constraints,
                             verbose     = FALSE)
    pm     <- .shortburst_plan_matrix(sb)
    scores <- .score_plans_seats(pm, bgs_dt, group_col, total_col,
                                  num_members, m_offset, ndists)
    best_j <- which.max(scores)
    list(plan       = as.integer(pm[, best_j]),
         score      = unname(scores[best_j]),
         seed_score = unname(scores[1]),
         n_plans    = ncol(pm),
         maximize   = TRUE,
         target     = target,
         k          = NA_integer_)
  }

  #Scalar (total-seats) objectives. 
  scalar_specs <- list(
    list(target = "black_seats_max",   group_col = "sum_BVAP",  total_col = "sum_VAP",   m_offset = 1L),
    list(target = "black_seats37_max", group_col = "sum_BVAP",  total_col = "sum_VAP",   m_offset = 2L),
    list(target = "dem_seats_max",     group_col = "sum_biden", total_col = "two_party", m_offset = 1L),
    list(target = "rep_seats_max",     group_col = "sum_trump", total_col = "two_party", m_offset = 1L)
  )
  for(spec in scalar_specs){
    if(spec$target %in% names(out)) next
    sf_score <- make_seats_scorer(pa_map, spec$group_col, spec$total_col,
                                   num_members, spec$m_offset)
    sd_seed  <- pick_seed(spec$target)
    message(sprintf("d=%02d%s %s (m_offset=%d, group=%s/%s, seed_score=%s) ...",
                    ndists, suffix, spec$target, spec$m_offset,
                    spec$group_col, spec$total_col,
                    format(seeds[[spec$target]]$score, digits = 4)))
    out[[spec$target]] <- one_sweep_seats(sf_score, spec$group_col,
                                           spec$total_col, spec$m_offset,
                                           spec$target, sd_seed)
    save_now()
  }

  #Competitiveness objectives: per-district distance to the nearest Droop
  #boundary, summarized as share-of-districts within 5pp (marg5pt_max) and
  #mean per-district flip probability under a normal swing of sigma=0.03
  #(flipprob_max). Same shortburst budget / constraints as the seat
  #scalars above; results feed sb_overlay and surface as the gold-marker
  #shortburst overlay in fig-cd-box.
  one_sweep_competitiveness <- function(sf_score, metric, target, init_plan_use){
    sb <- redist_shortburst(pa_map,
                             counties    = sb_counties,
                             score_fn    = sf_score,
                             maximize    = TRUE,
                             max_bursts  = max_bursts,
                             burst_size  = burst_size,
                             return_all  = TRUE,
                             init_plan   = init_plan_use,
                             constraints = constraints,
                             verbose     = FALSE)
    pm     <- .shortburst_plan_matrix(sb)
    scores <- .score_plans_competitiveness(pm, bgs_dt, num_members, metric)
    best_j <- which.max(scores)
    list(plan       = as.integer(pm[, best_j]),
         score      = unname(scores[best_j]),
         seed_score = unname(scores[1]),
         n_plans    = ncol(pm),
         maximize   = TRUE,
         target     = target,
         k          = NA_integer_)
  }
  comp_specs <- list(
    list(target = "marg5pt_max",  metric = "share_marginal"),
    list(target = "flipprob_max", metric = "flip_prob")
  )
  for(spec in comp_specs){
    if(spec$target %in% names(out)) next
    sf_score <- make_competitiveness_scorer(pa_map, num_members, spec$metric)
    sd_seed  <- pick_seed(spec$target)
    message(sprintf("d=%02d%s %s (metric=%s, seed_score=%s) ...",
                    ndists, suffix, spec$target, spec$metric,
                    format(seeds[[spec$target]]$score, digits = 4)))
    out[[spec$target]] <- one_sweep_competitiveness(sf_score, spec$metric,
                                                    spec$target, sd_seed)
    save_now()
  }

  #Vote-share frontier objectives
  share_specs <- list(
    list(target = "dem_max", group_col = "sum_biden", total_col = "two_party", maximize = TRUE),
    list(target = "dem_min", group_col = "sum_biden", total_col = "two_party", maximize = FALSE)
  )
  for(k in seq_len(ndists)){
    saved_this_k <- FALSE
    for(spec in share_specs){
      key <- paste0(spec$target, "_k", k)
      if(key %in% names(out)) next
      sf_score <- scorer_group_pct(pa_map, group_pop = sum_biden,
                                    total_pop = two_party, k = k)
      sd_seed  <- pick_seed(spec$target, k = k)
      message(sprintf("d=%02d%s k=%02d/%02d %s (seed_score=%s) ...",
                      ndists, suffix, k, ndists, spec$target,
                      format(seeds[[spec$target]][[k]]$score, digits = 4)))
      out[[key]] <- one_sweep_topk(sf_score, spec$group_col, spec$total_col,
                                    k, spec$maximize, spec$target, sd_seed)
      saved_this_k <- TRUE
    }
    #Incremental atomic save after each k 
    if(saved_this_k) save_now()
  }
  message("wrote ", cache_f)
}

#Loops ensure_shortburst_one across every config in district_specs.
ensure_shortburst_all <- function(admin_col = NULL, ...){
  for(spec in district_specs) ensure_shortburst_one(spec$ndists, admin_col, ...)
  invisible()
}

# Supplementary shortburst sweeps
# A motivated BlackSeats / BlackSeats37 sweep applied uniformly to every
# configuration, beyond what the regular ensure_shortburst_all covers:
ensure_extra_shortburst_black <- function(ndists,
                                          out_dir    = "./Output Data",
                                          max_bursts = 5000L,
                                          burst_size = 20L,
                                          admin_col  = "county_cousub",
                                          overwrite  = FALSE,
                                          log_file   = NULL){
  #Marker is per-(config, ensemble). Legacy county_cousub markers were written
  #WITHOUT an admin tag; keep that name for back-compat so existing caches still
  #gate (and a fresh county_cousub render stays fast). Other variants
  #(e.g. ccstruct_m5) get a tagged name so they don't collide.
  marker_name <- if(is.null(admin_col) || identical(admin_col, "county_cousub"))
                   sprintf("extra_shortburst_d%02d_black.Rdata", ndists)
                 else
                   sprintf("extra_shortburst_d%02d_%s_black.Rdata", ndists, admin_col)
  out_f <- file.path(out_dir, marker_name)
  if(file.exists(out_f) && !overwrite){
    message(sprintf("ensure_extra_shortburst_black: cache present at %s -- skipping.",
                    out_f))
    return(invisible(NULL))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  con <- NULL
  if(!is.null(log_file)){
    con <- file(log_file, open = "a")
    sink(con, type = "output",  split = TRUE)
    sink(con, type = "message", append = TRUE)
    on.exit({ sink(type = "message"); sink(type = "output"); close(con) }, add = TRUE)
  }

  num_members <- 18 / ndists
  if(num_members < 2) num_members <- 1
  num_members <- as.integer(num_members)
  message(sprintf("ensure_extra_shortburst_black: start d=%02d (M=%d) %s pid=%d",
                  ndists, num_members, format(Sys.time()), Sys.getpid()))

  smc_dir   <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, admin_col))
  smc_files <- sort(list.files(smc_dir, pattern = "^batch_.*\\.Rdata$",
                                full.names = TRUE))
  smc_files <- smc_files[!grepl("\\.tmp$", smc_files)]
  if(length(smc_files) == 0L)
    stop("ensure_extra_shortburst_black: no SMC batches under ", smc_dir)
  message(sprintf("  found %d SMC batch files", length(smc_files)))

  e <- new.env(); load(smc_files[1], envir = e)
  init_plan <- as.integer(.plans_matrix(e$plans_obj)[, 2])
  rm(e); gc(verbose = FALSE)

  load(dataset_cache)  # bgs
  bgs$init_plan <- init_plan
  pa_map <- redist::redist_map(
    bgs[, c("GEOID","sum_TotalPop","sum_VAP","sum_BVAP",
             "sum_biden","sum_trump",
             "county_fips","cousub_id","init_plan")],
    total_pop     = sum_TotalPop,
    pop_tol       = 0.01,
    existing_plan = init_plan)
  pa_map$two_party <- pa_map$sum_biden + pa_map$sum_trump

  #Match the ensemble's admin constraints (see ensure_shortburst_one).
  constraints <- redist::redist_constr(pa_map)
  sb_counties <- NULL
  if(identical(admin_col, "ccstruct_m5")){
    #NEW main: HARD county (counties=) + SOFT cousub strength 5.
    sb_counties <- pa_map$county_fips
    constraints <- redist::add_constr_splits(constraints, strength = 5,
                                              admin = cousub_id)
  }else{
    #county_cousub (OLD main) and any soft-both-levels variant: county 1 / cousub 3.
    constraints <- redist::add_constr_splits(constraints, strength = 1,
                                              admin = county_fips)
    constraints <- redist::add_constr_splits(constraints, strength = 3,
                                              admin = cousub_id)
  }

  bgs_dt <- sf::st_drop_geometry(bgs)
  bgs_dt$two_party <- bgs_dt$sum_biden + bgs_dt$sum_trump

  bvap <- as.numeric(bgs$sum_BVAP); bvap[is.na(bvap)] <- 0
  vap  <- as.numeric(bgs$sum_VAP);  vap[is.na(vap)]  <- 0
  group_pct_c <- get("group_pct", envir = asNamespace("redist"))

  message(sprintf("  scoring %d batches to pick 5 diverse seeds ...",
                  length(smc_files)))
  all_pm     <- vector("list", length(smc_files))
  all_scores <- list()
  for(i in seq_along(smc_files)){
    f <- smc_files[i]
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
    if(!ok || !exists("plans_obj", envir = e)){ rm(e); next }
    pm <- redist::get_plans_matrix(e$plans_obj)
    pm <- pm[, -1, drop = FALSE]
    all_pm[[i]] <- pm
    pct_bvap   <- group_pct_c(pm, bvap, vap, as.integer(ndists))
    sorted_bvap<- apply(pct_bvap, 2, sort, decreasing = TRUE)
    if(!is.matrix(sorted_bvap)) sorted_bvap <- matrix(sorted_bvap, nrow = ndists)
    all_scores[[i]] <- data.frame(
      batch_idx = i,
      col       = seq_len(ncol(pm)),
      bs1       = colSums(floor(pct_bvap * (num_members + 1L))),
      bs2       = colSums(floor(pct_bvap * (num_members + 2L))),
      top1      = sorted_bvap[1, ],
      top2      = sorted_bvap[2, ],
      sum_top2  = sorted_bvap[1, ] + sorted_bvap[2, ]
    )
    rm(e, pct_bvap, sorted_bvap); gc(verbose = FALSE)
  }
  S <- do.call(rbind, all_scores)
  message(sprintf("  scored %d plans total", nrow(S)))

  pick_seed <- function(criterion) S[which.max(S[[criterion]]), , drop = FALSE]
  seed_rows <- list(
    best_bs   = pick_seed("bs1"),
    best_bs2  = pick_seed("bs2"),
    best_top1 = pick_seed("top1"),
    best_top2 = pick_seed("top2"),
    best_sum  = pick_seed("sum_top2")
  )
  for(nm in names(seed_rows)){
    s <- seed_rows[[nm]]
    message(sprintf("    %-10s  batch=%d col=%d  bs1=%d  bs2=%d  top1=%.3f  top2=%.3f",
                    nm, s$batch_idx, s$col, s$bs1, s$bs2, s$top1, s$top2))
  }

  run_one <- function(init_plan_use, m_offset, label){
    sf_score <- make_seats_scorer(pa_map, "sum_BVAP", "sum_VAP",
                                   num_members, m_offset)
    t0 <- Sys.time()
    sb <- redist::redist_shortburst(pa_map,
                                     counties    = sb_counties,
                                     score_fn    = sf_score,
                                     maximize    = TRUE,
                                     max_bursts  = max_bursts,
                                     burst_size  = burst_size,
                                     return_all  = TRUE,
                                     init_plan   = init_plan_use,
                                     constraints = constraints,
                                     verbose     = FALSE)
    pm_sb   <- .shortburst_plan_matrix(sb)
    scores  <- .score_plans_seats(pm_sb, bgs_dt, "sum_BVAP", "sum_VAP",
                                   num_members, m_offset, ndists)
    best_j  <- which.max(scores)
    wall    <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("    %s: seed=%d  best=%d  n_plans=%d  wall=%.1f min",
                    label, scores[1], scores[best_j], ncol(pm_sb), wall))
    list(plan       = as.integer(pm_sb[, best_j]),
         score      = unname(scores[best_j]),
         seed_score = unname(scores[1]),
         n_plans    = ncol(pm_sb),
         maximize   = TRUE,
         target     = sub("_m1$|_m2$", "", label),
         k          = NA_integer_)
  }

  results <- list()
  for(seed_name in names(seed_rows)){
    s  <- seed_rows[[seed_name]]
    ip <- as.integer(all_pm[[s$batch_idx]][, s$col])
    message(sprintf("  === seed=%s (bs1=%d, bs2=%d) ===",
                    seed_name, s$bs1, s$bs2))
    results[[paste0(seed_name, "_bs")]]  <- run_one(ip, 1L, paste0(seed_name, "_bs"))
    results[[paste0(seed_name, "_bs2")]] <- run_one(ip, 2L, paste0(seed_name, "_bs2"))
  }

  pick_best <- function(suffix){
    hits <- results[grepl(paste0("_", suffix, "$"), names(results))]
    scores <- vapply(hits, function(x) x$score, numeric(1))
    best_nm <- names(hits)[which.max(scores)]
    list(name = best_nm, entry = hits[[best_nm]])
  }
  best_bs  <- pick_best("bs")
  best_bs2 <- pick_best("bs2")
  message(sprintf("  overall best BlackSeats:    %s -> score=%d",
                  best_bs$name,  best_bs$entry$score))
  message(sprintf("  overall best BlackSeats37:  %s -> score=%d",
                  best_bs2$name, best_bs2$entry$score))

  cache_f <- file.path(bg_sb_dir, sprintf("shortburst_bg_d%02d_%s.Rdata",
                                           ndists, admin_col))
  e_cache <- new.env(); load(cache_f, envir = e_cache)
  out <- e_cache$out
  prev_bs  <- out$black_seats_max$score
  prev_bs2 <- out$black_seats37_max$score
  message(sprintf("  previous cache: black_seats_max=%s, black_seats37_max=%s",
                  prev_bs %||% NA, prev_bs2 %||% NA))
  if(is.null(prev_bs)  || best_bs$entry$score  > prev_bs){
    out$black_seats_max   <- best_bs$entry
    out$black_seats_max$target <- "black_seats_max"
    message("  UPDATED black_seats_max in main cache.")
  }
  if(is.null(prev_bs2) || best_bs2$entry$score > prev_bs2){
    out$black_seats37_max <- best_bs2$entry
    out$black_seats37_max$target <- "black_seats37_max"
    message("  UPDATED black_seats37_max in main cache.")
  }
  tmp_f <- paste0(cache_f, ".tmp")
  save(out, file = tmp_f, compress = "xz")
  file.rename(tmp_f, cache_f)

  marker <- list(
    ndists      = ndists,
    num_members = num_members,
    admin_col   = admin_col,
    max_bursts  = max_bursts,
    burst_size  = burst_size,
    n_seeds     = length(seed_rows),
    seed_rows   = seed_rows,
    results     = results,
    best_bs     = best_bs,
    best_bs2    = best_bs2,
    built_at    = format(Sys.time()))
  tmp_f <- paste0(out_f, ".tmp")
  save(marker, file = tmp_f, compress = "xz")
  file.rename(tmp_f, out_f)
  message(sprintf("ensure_extra_shortburst_black: wrote %s", out_f))
  invisible(marker)
}

#Apply the motivated BlackSeats / BlackSeats37 sweep uniformly to every
#configuration in district_specs
ensure_extra_shortburst_all_black <- function(...){
  for(spec in district_specs)
    ensure_extra_shortburst_black(spec$ndists, ...)
  invisible()
}

# Metrics functions
sd_metric <- function(num, den, target){
  r <- num / den
  sqrt(colMeans((r - target)^2))
}

#Compute plan-level + per-unit metrics across an SMC ensemble.
generate_metrics_bg <- function(metrics, metric_vals, ndists, bgs.full,
                                admin_col = NULL,
                                cache = TRUE, overwrite = FALSE){
  num_districts <- ndists
  num_members   <- 18/num_districts
  if(num_members < 2) num_members <- 1
  
  bgs <- st_drop_geometry(bgs.full)
  bgs$Voters <- bgs$sum_biden + bgs$sum_trump
  
  #Stack every numerator/denominator into one nbgs × K matrix so each
  #rowsum() call produces all per-district sums
  col_names <- c("sum_biden","Voters","sum_BVAP","sum_VAP",
                 "NHBlack","TotalPop_ACS","InPoverty","PovDeterm",
                 "POP65Plus","sum_TotalPop","HHwithChld","Households",
                 "EmpManufac","EmpTotal")
  X <- as.matrix(bgs[, col_names])
  X[is.na(X)] <- 0   # mirrors na.rm = TRUE in the per-district sums
  ci   <- setNames(seq_along(col_names), col_names)
  nbgs <- nrow(X)
  
  #Target index lookups
  t_PctBlack         <- metric_vals[match("PctBlack",         metrics)]
  t_PctPoverty       <- metric_vals[match("PctPoverty",       metrics)]
  t_Pct65Plus        <- metric_vals[match("Pct65Plus",        metrics)]
  t_PctHHwChild      <- metric_vals[match("PctHHwChild",      metrics)]
  t_PctManufacturing <- metric_vals[match("PctManufacturing", metrics)]
  t_ShareDem         <- metric_vals[match("ShareDem",         metrics)]
  
  bg_measures <- data.frame(GEOID                = bgs$GEOID,
                            CompetitiveDistricts = 0,
                            CompetitiveSeats     = 0,
                            SeniorsSD            = 0,
                            PovertySD            = 0,
                            DemShare             = 0,
                            BlackShare           = 0,
                            Wasted               = 0,
                            BlackWasted          = 0,
                            BlackMember          = 0,
                            BlackMember37        = 0)
  
  
  
  #Admin-split vectors for the per-plan admin_splits computation
  county_id <- bgs$county_fips
  muni_id   <- bgs$cousub_id
  ok_county <- !is.na(county_id)
  ok_muni   <- !is.na(muni_id)
  
  #Droop quota for this configuration
  drp_q <- 1 / (num_members + 1)
  
  plan_level  <- list()
  total_plans <- 0L
  
  #SMC batches. Dispatch by admin_col so we can score either the canonical
  #unconstrained ensemble or one of the admin-preserving comparison runs.
  if(is.null(admin_col)){
    smc_files <- list_smc_batches(ndists)
  }else{
    d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, admin_col))
    #Permissive pattern: matches both legacy batch_NNN.Rdata (single-writer
    #ensembles like the county_fips and cousub_id 10k runs) and PID-tagged
    #batch_pidPPP_NNN.Rdata produced by parallel workers on the combined
    #county_cousub ensemble.
    smc_files <- sort(list.files(d_dir, pattern = "^batch_.*\\.Rdata$", full.names = TRUE))
    smc_files <- smc_files[!grepl("\\.tmp$", smc_files)]
  }

  #Per-batch metrics cache. Bumped when a stored column's *definition*
  #changes so old caches automatically invalidate on next source.
  #  v2: prob_seat_change_sigma3 switched from 2 * Phi(-delta/sigma) to
  #      the exact two-boundary form via flip_prob_q().
  CACHE_VERSION <- 2L
  cache_dir <- NULL
  if(isTRUE(cache)){
    suffix <- if(is.null(admin_col)) "" else paste0("_", admin_col)
    cache_dir <- file.path(bg_root, "metrics", "om",
                           sprintf("d%02d%s", ndists, suffix))
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }

  #use new environment for easy cleanup afterwards and to keep
  #memory use to a minimum. Had some bugs introduced while running
  #on multiple machines connected by OneDrive. Probably extraneous
  #but it worked...
  for(bf in smc_files){
    cache_f <- if(!is.null(cache_dir)) file.path(cache_dir, basename(bf)) else NULL

    #Cache hit path: load (pl_row, bg_delta, nplans) and accumulate
    #without re-scoring the batch.
    if(!is.null(cache_f) && file.exists(cache_f) && !isTRUE(overwrite)){
      ce <- new.env()
      load_ok <- tryCatch({ load(cache_f, envir = ce); TRUE },
                          error = function(err) FALSE)
      if(load_ok &&
         all(c("pl_row","bg_delta","nplans","cache_version") %in% ls(ce)) &&
         identical(ce$cache_version, CACHE_VERSION)){
        plan_level[[length(plan_level)+1L]] <- ce$pl_row
        for(col in setdiff(colnames(ce$bg_delta), "GEOID"))
          bg_measures[[col]] <- bg_measures[[col]] + ce$bg_delta[[col]]
        total_plans <- total_plans + ce$nplans
        message(sprintf("  d=%02d batch %s: cached (%d plans)",
                        ndists, tools::file_path_sans_ext(basename(bf)), ce$nplans))
        next
      }
      #fall through to recompute on load failure or schema mismatch
    }

    e <- new.env()
    tryCatch(load(bf, envir = e),
             error = function(err){ message("    skipping (load error): ", basename(bf)); NULL })
    plans_obj <- if(exists("plans_obj", envir = e)) e$plans_obj
    else if(exists("pa_plans_bg", envir = e)) e$pa_plans_bg
    else NULL
    if(is.null(plans_obj)) next
    M <- .plans_matrix(plans_obj)
    M <- M[, -1, drop = FALSE]
    if(nrow(M) != nbgs)
      stop(sprintf("plan matrix nrow (%d) != BG rows (%d) for %s",
                   nrow(M), nbgs, basename(bf)))
    nplans <- ncol(M)
    t0 <- Sys.time()
    message(sprintf("  d=%02d batch %s: %d plans ...",
                    ndists, tools::file_path_sans_ext(basename(bf)), nplans))
    
    #Per-district sums for all K cols, packed into ndists × K × nplans
    PD <- array(0, dim = c(num_districts, length(col_names), nplans))
    for(p in seq_len(nplans)){
      rs <- rowsum(X, M[, p], reorder = TRUE)
      d_idx <- as.integer(rownames(rs))
      PD[d_idx, , p] <- rs
    }
    
    #ndists × nplans slices keyed by column name
    sb    <- PD[, ci["sum_biden"],    ]
    vt    <- PD[, ci["Voters"],       ]
    pb    <- PD[, ci["sum_BVAP"],     ]
    pvap  <- PD[, ci["sum_VAP"],      ]
    nhb   <- PD[, ci["NHBlack"],      ]
    tpacs <- PD[, ci["TotalPop_ACS"], ]
    pov   <- PD[, ci["InPoverty"],    ]
    povd  <- PD[, ci["PovDeterm"],    ]
    p65   <- PD[, ci["POP65Plus"],    ]
    tp    <- PD[, ci["sum_TotalPop"], ]
    hhc   <- PD[, ci["HHwithChld"],   ]
    hh    <- PD[, ci["Households"],   ]
    empm  <- PD[, ci["EmpManufac"],   ]
    empt  <- PD[, ci["EmpTotal"],     ]
    
    #Plan-level scalars, vectorized over plans
    #Two-party Dem share per district
    dp <- sb / vt
    CompetitiveDistricts <- colSums((dp > 0.45) & (dp < 0.55))
    
    #Wasted votes / Dem seats — single-expression form of getWasted.
    dem_votes  <- sb
    rep_votes  <- vt - sb
    toWin      <- vt / (num_members + 1)
    dem_wins   <- floor(dem_votes / toWin)
    rep_wins   <- floor(rep_votes / toWin)
    dem_wasted <- dem_votes - dem_wins * toWin
    rep_wasted <- rep_votes - rep_wins * toWin
    
    DemSeats         <- colSums(dem_wins)
    RepSeats         <- colSums(rep_wins)
    ShareDemSeats    <- DemSeats / pmax(DemSeats + RepSeats, 1)
    WastedVotes      <- (colSums(dem_wasted) - colSums(rep_wasted)) / colSums(vt)
    CompetitiveSeats <- if(num_members == 1) CompetitiveDistricts else
      colSums((dem_wins != 0) & (rep_wins != 0))
    
    #Black vote metrics — same algebra, two thresholds (toWin, toWin2).
    black_votes  <- pb
    other_votes  <- pvap - pb
    bk_toWin     <- pvap / (num_members + 1)
    bk_toWin2    <- pvap / (num_members + 2)
    black_wins   <- floor(black_votes / bk_toWin)
    black_wins37 <- floor(black_votes / bk_toWin2)
    other_wins   <- floor(other_votes / bk_toWin)
    black_wasted <- black_votes - black_wins * bk_toWin
    other_wasted <- other_votes - other_wins * bk_toWin
    
    BlackSeats     <- colSums(black_wins)
    BlackSeats37   <- colSums(black_wins37)
    BlackWasted    <- (colSums(black_wasted) - colSums(other_wasted)) / colSums(pvap)
    BlackPctWasted <- (colSums(black_wasted) / colSums(black_votes)) -
      (colSums(other_wasted) / colSums(other_votes))
    
    #Heterogeneity SDs (one number per plan per metric).
    het_PctBlack         <- sd_metric(nhb,  tpacs, t_PctBlack)
    het_PctPoverty       <- sd_metric(pov,  povd,  t_PctPoverty)
    het_Pct65Plus        <- sd_metric(p65,  tp,    t_Pct65Plus)
    het_PctHHwChild      <- sd_metric(hhc,  hh,    t_PctHHwChild)
    het_PctManufacturing <- sd_metric(empm, empt,  t_PctManufacturing)
    het_ShareDem         <- sd_metric(sb,   vt,    t_ShareDem)
    
    #Competitiveness
    #q = 1/(M+1); delta = distance from each district's Dem share to
    #the nearest Droop boundary. Plan-level summaries:
    #  mean_competitiveness_q  in [0,1], 1 = at threshold
    #  prob_seat_change_sigma3 = expected number of seats in play
    #                           under a normal swing of sigma=0.03
    v_modq <- dp %% drp_q
    delta  <- pmin(v_modq, drp_q - v_modq)
    n_marginal_2pt          <- colSums(delta < 0.02) #'really close' to flipping
    n_marginal_5pt          <- colSums(delta < 0.05) #conventional 'competitive'
    mean_competitiveness_q  <- colMeans(1 - 2*(num_members+1)*delta)
    #Exact two-boundary flip probability; see flip_prob_q() definition.
    prob_seat_change_sigma3 <- colMeans(flip_prob_q(delta, drp_q, 0.03))
    
    #Admin-unit splits
    #Count of counties / municipalities that span more than one district per plan.
    n_county_splits <- integer(nplans)
    n_muni_splits   <- integer(nplans)
    for(p in seq_len(nplans)){
      if(any(ok_county)){
        tab_c <- table(county_id[ok_county], M[ok_county, p])
        n_county_splits[p] <- sum(rowSums(tab_c > 0) > 1)
      }else n_county_splits[p] <- NA_integer_
      if(any(ok_muni)){
        tab_m <- table(muni_id[ok_muni], M[ok_muni, p])
        n_muni_splits[p] <- sum(rowSums(tab_m > 0) > 1)
      }else n_muni_splits[p] <- NA_integer_
    }
    
    pl_row <- data.frame(
      Plan                 = total_plans + seq_len(nplans),
      PctBlack             = het_PctBlack,
      PctPoverty           = het_PctPoverty,
      Pct65Plus            = het_Pct65Plus,
      PctHHwChild          = het_PctHHwChild,
      PctManufacturing     = het_PctManufacturing,
      ShareDem             = het_ShareDem,
      DemSeats             = DemSeats,
      ShareDemSeats        = ShareDemSeats,
      CompetitiveSeats     = CompetitiveSeats,
      CompetitiveDistricts = CompetitiveDistricts,
      WastedVotes          = WastedVotes,
      BlackSeats           = BlackSeats,
      BlackSeats37         = BlackSeats37,
      BlackWasted          = BlackWasted,
      BlackPctWasted       = BlackPctWasted,
      #NEW: Droop competitiveness + admin splits
      n_marginal_2pt          = n_marginal_2pt,
      n_marginal_5pt          = n_marginal_5pt,
      mean_competitiveness_q  = mean_competitiveness_q,
      prob_seat_change_sigma3 = prob_seat_change_sigma3,
      n_county_splits         = n_county_splits,
      n_muni_splits           = n_muni_splits
    )
    plan_level[[length(plan_level)+1L]] <- pl_row

    #Per-unit accumulators
    seniors_abs <- abs(p65 / tp     - t_Pct65Plus)
    poverty_abs <- abs(pov / povd   - t_PctPoverty)
    wasted_per  <- (dem_wasted   - rep_wasted)   / vt
    blkwst_per  <- (black_wasted - other_wasted) / pvap
    bshare_per  <- pb / pvap
    comp_per    <- (dp > 0.45) & (dp < 0.55)
    multi_per   <- (dem_wins != 0) & (rep_wins != 0)

    #Build this batch's bg_delta separately from the running bg_measures
    #total so it can be cached and reloaded standalone on later sources.
    bg_delta <- data.frame(
      GEOID                = bgs$GEOID,
      CompetitiveDistricts = 0,
      CompetitiveSeats     = 0,
      SeniorsSD            = 0,
      PovertySD            = 0,
      DemShare             = 0,
      BlackShare           = 0,
      Wasted               = 0,
      BlackWasted          = 0,
      BlackMember          = 0,
      BlackMember37        = 0
    )
    for(p in seq_len(nplans)){
      assign <- M[, p]
      bg_delta$DemShare             <- bg_delta$DemShare             + dp[assign, p]
      bg_delta$CompetitiveDistricts <- bg_delta$CompetitiveDistricts + comp_per[assign, p]
      bg_delta$CompetitiveSeats     <- bg_delta$CompetitiveSeats +
        if(num_members == 1) comp_per[assign, p] else multi_per[assign, p]
      bg_delta$Wasted        <- bg_delta$Wasted        + wasted_per[assign, p]
      bg_delta$BlackWasted   <- bg_delta$BlackWasted   + blkwst_per[assign, p]
      bg_delta$BlackShare    <- bg_delta$BlackShare    + bshare_per[assign, p]
      bg_delta$BlackMember   <- bg_delta$BlackMember   + black_wins[assign, p]
      bg_delta$BlackMember37 <- bg_delta$BlackMember37 + black_wins37[assign, p]
      bg_delta$SeniorsSD     <- bg_delta$SeniorsSD     + seniors_abs[assign, p]
      bg_delta$PovertySD     <- bg_delta$PovertySD     + poverty_abs[assign, p]
    }
    for(col in setdiff(colnames(bg_delta), "GEOID"))
      bg_measures[[col]] <- bg_measures[[col]] + bg_delta[[col]]

    total_plans <- total_plans + nplans

    #Cache this batch's contribution 
    if(!is.null(cache_f)){
      cache_version <- CACHE_VERSION
      tmp_f <- paste0(cache_f, ".tmp")
      save(pl_row, bg_delta, nplans, cache_version,
           file = tmp_f, compress = "xz")
      file.rename(tmp_f, cache_f)
    }

    message(sprintf("    done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    rm(plans_obj, M, PD, sb, vt, pb, pvap, nhb, tpacs, pov, povd, p65, tp,
       hhc, hh, empm, empt, dp, dem_wins, rep_wins, dem_wasted, rep_wasted,
       black_wins, black_wins37, other_wins, black_wasted, other_wasted,
       seniors_abs, poverty_abs, wasted_per, blkwst_per, bshare_per,
       comp_per, multi_per)
    gc(verbose = FALSE)
  }
  
  output_metrics <- bind_rows(plan_level)
  output_metrics$Total <- rowSums(output_metrics[, c("PctPoverty","Pct65Plus","PctHHwChild","PctManufacturing")])

  bg_measures[, -1] <- bg_measures[, -1] / total_plans
  list(output_metrics, bg_measures)
}

#Polsby-Popper compactness scoring across an SMC ensemble.
generate_polsby_bg <- function(ndists, bgs.full, admin_col = NULL,
                                perim_df = NULL,
                                cache = TRUE, overwrite = FALSE){
  if(is.null(perim_df))
    perim_df <- redistmetrics::prep_perims(shp = bgs.full)
  if(is.null(admin_col)){
    smc_files <- list_smc_batches(ndists)
  }else{
    d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, admin_col))
    smc_files <- sort(list.files(d_dir, pattern = "^batch_.*\\.Rdata$",
                                  full.names = TRUE))
    smc_files <- smc_files[!grepl("\\.tmp$", smc_files)]
  }
  if(length(smc_files) == 0L) return(data.frame())

  #Per-batch cache. Necessary to deal with clumsy parallelization across multiple computers.
  CACHE_VERSION <- 1L
  cache_dir <- NULL
  if(isTRUE(cache)){
    suffix <- if(is.null(admin_col)) "" else paste0("_", admin_col)
    cache_dir <- file.path(bg_root, "metrics", "polsby",
                            sprintf("d%02d%s", ndists, suffix))
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }

  pp_rows <- list()
  for(bf in smc_files){
    cache_f <- if(!is.null(cache_dir)) file.path(cache_dir, basename(bf)) else NULL

    #Cache hit path
    if(!is.null(cache_f) && file.exists(cache_f) && !isTRUE(overwrite)){
      ce <- new.env()
      load_ok <- tryCatch({ load(cache_f, envir = ce); TRUE },
                          error = function(err) FALSE)
      if(load_ok &&
         all(c("pp_batch","nplans","cache_version") %in% ls(ce)) &&
         identical(ce$cache_version, CACHE_VERSION)){
        pp_rows[[length(pp_rows)+1L]] <- ce$pp_batch
        message(sprintf("  polsby d=%02d batch %s: cached (%d plans)",
                        ndists, tools::file_path_sans_ext(basename(bf)), ce$nplans))
        next
      }
    }

    e <- new.env()
    tryCatch(load(bf, envir = e),
             error = function(err){
               message("    skipping (load error): ", basename(bf)); NULL })
    plans_obj <- if(exists("plans_obj", envir = e)) e$plans_obj
                 else if(exists("pa_plans_bg", envir = e)) e$pa_plans_bg
                 else NULL
    if(is.null(plans_obj)) next
    pm <- .plans_matrix(plans_obj)
    pm <- pm[, -1, drop = FALSE]
    nplans <- ncol(pm)

    t0 <- Sys.time()
    message(sprintf("  polsby d=%02d batch %s: %d plans ...",
                    ndists, tools::file_path_sans_ext(basename(bf)), nplans))
    pp_raw <- redistmetrics::comp_polsby(plans   = pm,
                                          shp     = bgs.full,
                                          perim_df = perim_df)
    if(is.data.frame(pp_raw)){
      pp_mat <- matrix(pp_raw[["polsby"]], nrow = ndists, ncol = nplans)
    }else{
      pp_mat <- matrix(as.numeric(pp_raw), nrow = ndists, ncol = nplans)
    }

    pp_batch <- data.frame(
      mean_polsby_popper = colMeans(pp_mat, na.rm = TRUE),
      min_polsby_popper  = apply(pp_mat, 2, min, na.rm = TRUE),
      max_polsby_popper  = apply(pp_mat, 2, max, na.rm = TRUE)
    )

    if(!is.null(cache_f)){
      cache_version <- CACHE_VERSION
      tmp_f <- paste0(cache_f, ".tmp")
      save(pp_batch, nplans, cache_version, file = tmp_f, compress = "xz")
      file.rename(tmp_f, cache_f)
    }
    pp_rows[[length(pp_rows)+1L]] <- pp_batch

    wall <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("    wall=%.2f min, median PP=%.3f",
                    wall, median(pp_batch$mean_polsby_popper, na.rm = TRUE)))
  }

  do.call(rbind, pp_rows)
}

#Block-group-level diagnostics for one enacted-plan assignment column
#(e.g. "RoughCD116" or "RoughCD118"). Returns a one-row summary list with
#district count, county/municipality splits, mean Polsby-Popper compactness,
#and single-member Droop Democratic seat share. Used by the enacted-plan
#comparison in Supplementary Materials Section 3; results are precomputed in
#BBW - Data Prep BG.R and saved to ./Output Data/enacted_metrics.Rdata.
enacted_metrics <- function(col, bgs.full, perim_df = NULL){
  if(is.null(perim_df))
    perim_df <- redistmetrics::prep_perims(shp = bgs.full)
  plan <- as.integer(bgs.full[[col]])
  ok   <- !is.na(plan)
  cnty <- bgs.full$county_fips
  muni <- bgs.full$cousub_id
  okc  <- ok & !is.na(cnty)
  okm  <- ok & !is.na(muni)
  #a county / municipality is "split" if its block groups span >1 district
  n_cnty <- sum(rowSums(table(cnty[okc], plan[okc]) > 0) > 1)
  n_muni <- sum(rowSums(table(muni[okm], plan[okm]) > 0) > 1)
  #single-member Droop seat share (toWin = half the two-party vote)
  vd <- tapply(bgs.full$sum_biden[ok], plan[ok], sum, na.rm = TRUE)
  vr <- tapply(bgs.full$sum_trump[ok], plan[ok], sum, na.rm = TRUE)
  toWin <- (vd + vr) / 2
  dem <- floor(vd / toWin); rep <- floor(vr / toWin)
  dem_share <- sum(dem) / sum(dem + rep)
  #Polsby-Popper needs a complete assignment; the single unpopulated NA
  #block group is folded into district 1 (negligible effect on any district).
  pv <- plan; pv[is.na(pv)] <- 1L
  pp <- as.numeric(redistmetrics::comp_polsby(plans = matrix(pv, ncol = 1),
                                              shp = bgs.full, perim_df = perim_df))
  list(ndists = length(unique(plan[ok])), n_county = n_cnty, n_muni = n_muni,
       pp_mean = mean(pp, na.rm = TRUE), dem_share = dem_share)
}

# ============================================================================
# Penalty-sweep tools (Supplementary Materials Section 3 robustness check)
# ----------------------------------------------------------------------------
# A "variant" is a directory suffix encoding a penalty cell, e.g. "cc_c8_m12"
# (combined county+cousub, county_strength 8, cousub_strength 12).
# ensure_smc_one_combined(..., dir_suffix = variant) writes PID-tagged SMC
# batches to smc_plans/d{NN}_<variant>/; generate_metrics_bg()/
# generate_polsby_bg() with admin_col = variant score those plans and cache
# under metrics/om/d{NN}_<variant>/ and metrics/polsby/d{NN}_<variant>/.
# The "county_cousub" main ensemble is itself a variant by this convention.
# Readers ignore *.Rdata.tmp and tryCatch() loads so a OneDrive mid-sync file
# is skipped rather than fatal. These helpers resolve bg_root / bg_smc_dir by
# lexical scope from the calling script, same as the other ensemble helpers.
# ============================================================================

#Configs in scope for the sweep (d=50 dropped). Headline metrics tracked.
VARIANT_CONFIGS <- c(3, 6, 9, 17, 18)
VARIANT_METRICS <- c("ShareDemSeats", "BlackSeats", "BlackSeats37",
                     "WastedVotes", "mean_competitiveness_q",
                     "prob_seat_change_sigma3", "n_marginal_5pt",
                     "n_county_splits", "n_muni_splits")

#Tail-vs-shortburst mapping: ensemble metric <- sb_overlay target (same column).
.TAIL_MAP <- list(
  list(metric = "ShareDemSeats",           target = "dem_seats_max"),
  list(metric = "BlackSeats",              target = "black_seats_max"),
  list(metric = "BlackSeats37",            target = "black_seats37_max"),
  list(metric = "n_marginal_5pt",          target = "marg5pt_max"),
  list(metric = "prob_seat_change_sigma3", target = "flipprob_max")
)

#List per-batch cache files for a variant/config under metrics/<root_sub>/,
#ordered by mtime (≈ landing order). Ignores *.tmp.
.list_variant_files <- function(root_sub, variant, ndists){
  d_dir <- file.path(bg_root, root_sub, sprintf("d%02d_%s", ndists, variant))
  if(!dir.exists(d_dir)) return(character(0))
  fs <- list.files(d_dir, pattern = "\\.Rdata$", full.names = TRUE)
  fs <- fs[!grepl("\\.tmp$", fs)]
  if(!length(fs)) return(fs)
  fs[order(file.info(fs)$mtime)]
}

#Per-batch plan-level metric tables (pl_row), ordered named list.
read_variant_metric_batches <- function(variant, ndists){
  fs <- .list_variant_files(file.path("metrics", "om"), variant, ndists)
  if(!length(fs)) return(list())
  out <- lapply(fs, function(f){
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
    if(!ok || !exists("pl_row", envir = e)) return(NULL)
    e$pl_row
  })
  names(out) <- tools::file_path_sans_ext(basename(fs))
  out[!vapply(out, is.null, logical(1))]
}

#All plan-level metrics for a variant/config, row-bound.
read_variant_metrics <- function(variant, ndists){
  b <- read_variant_metric_batches(variant, ndists)
  if(!length(b)) return(NULL)
  dplyr::bind_rows(b)
}

#All Polsby-Popper rows (mean/min/max per plan) for a variant/config.
read_variant_polsby <- function(variant, ndists){
  fs <- .list_variant_files(file.path("metrics", "polsby"), variant, ndists)
  if(!length(fs)) return(NULL)
  rows <- lapply(fs, function(f){
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
    if(!ok || !exists("pp_batch", envir = e)) return(NULL)
    e$pp_batch
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if(!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

#SMC batch markers (used strengths / pop_temper / wall time) for a
#variant/config. NOTE: load()s the full SMC batch files, so I/O-heavy;
#`max_files` caps how many of the most-recent batches are inspected.
read_variant_smc_markers <- function(variant, ndists, max_files = Inf){
  d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, variant))
  if(!dir.exists(d_dir)) return(NULL)
  fs <- list.files(d_dir, pattern = "^batch_.*\\.Rdata$", full.names = TRUE)
  fs <- fs[!grepl("\\.tmp$", fs)]
  if(!length(fs)) return(NULL)
  fs <- fs[order(file.info(fs)$mtime)]
  if(is.finite(max_files) && length(fs) > max_files) fs <- tail(fs, max_files)
  rows <- lapply(fs, function(f){
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
    if(!ok) return(NULL)
    getn <- function(n, d = NA_real_) if(exists(n, envir = e)) get(n, envir = e) else d
    data.frame(file = basename(f), this_nsims = getn("this_nsims"),
               used_pop_temper = getn("used_pop_temper"),
               used_cnty_str = getn("used_county_strength"),
               used_cousub_str = getn("used_cousub_strength"),
               pop_tol = getn("pop_tol"), wall_min = getn("wall_min"),
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if(!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

#Enacted reference (m116 = 18 districts, m118 = 17 districts). NULL when no
#enacted plan exists at this magnitude or the cache is missing.
.enacted_for <- function(ndists){
  f <- "./Output Data/enacted_metrics.Rdata"
  if(!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  if(ndists == 18 && exists("m116", envir = e)) return(e$m116)
  if(ndists == 17 && exists("m118", envir = e)) return(e$m118)
  NULL
}

#Build (cache-gated) one penalty cell: SMC into the variant dir, then metrics +
#Polsby-Popper scored off those plans. No-ops where caches already exist.
#pop_tol is fixed at 0.01; the fallback ladder climbs pop_temper only (no
#strength-drop step), so a 0.03 failure is a real "cannot solve at 1%" signal.
ensure_penalty_sweep_cell <- function(county_strength, cousub_strength, variant,
                                      configs = c(17, 18),
                                      target_nsims = 2000L, batch_size = 2000L,
                                      metrics, metric_vals, bgs, perim_df = NULL,
                                      ncores = max(1L, parallel::detectCores() - 1L)){
  if(is.null(perim_df)) perim_df <- redistmetrics::prep_perims(shp = bgs)
  for(nd in configs){
    ensure_smc_one_combined(
      ndists = nd, target_nsims = target_nsims, batch_size = batch_size,
      county_strength = county_strength, cousub_strength = cousub_strength,
      pop_tol = 0.01, pop_temper = 0,
      fallback_ladder = list(list(pop_temper = 0.01),
                             list(pop_temper = 0.02),
                             list(pop_temper = 0.03)),
      dir_suffix = variant, ncores = ncores)
    generate_metrics_bg(metrics, metric_vals, ndists = nd, bgs, admin_col = variant)
    generate_polsby_bg(ndists = nd, bgs.full = bgs, admin_col = variant,
                       perim_df = perim_df)
  }
  invisible()
}

#Assign the per-config seed plan (same RoughCD116/RoughCD118/RoughSLDU nesting
#used by ensure_smc_one_combined), returned as an integer vector. Factored out
#so the county-structural variant below can reuse it.
.assign_bg_seed <- function(bgs, ndists){
  seed <- rep(1L, nrow(bgs))
  if(ndists == 3){
    seed[bgs$RoughCD116 > 6]  <- 2L; seed[bgs$RoughCD116 > 12] <- 3L
  }else if(ndists == 6){
    seed[bgs$RoughCD116 > 3]  <- 2L; seed[bgs$RoughCD116 > 6]  <- 3L
    seed[bgs$RoughCD116 > 9]  <- 4L; seed[bgs$RoughCD116 > 12] <- 5L
    seed[bgs$RoughCD116 > 15] <- 6L
  }else if(ndists == 9){
    seed[bgs$RoughCD116 > 2]  <- 2L; seed[bgs$RoughCD116 > 4]  <- 3L
    seed[bgs$RoughCD116 > 6]  <- 4L; seed[bgs$RoughCD116 > 8]  <- 5L
    seed[bgs$RoughCD116 > 10] <- 6L; seed[bgs$RoughCD116 > 12] <- 7L
    seed[bgs$RoughCD116 > 14] <- 8L; seed[bgs$RoughCD116 > 16] <- 9L
  }else if(ndists == 17){ seed <- as.integer(bgs$RoughCD118)
  }else if(ndists == 18){ seed <- as.integer(bgs$RoughCD116)
  }else if(ndists == 50){ seed <- as.integer(bgs$RoughSLDU)
  }else stop("Unsupported ndists: ", ndists)
  seed[is.na(seed)] <- 1L
  seed[seed < 1 | seed > ndists] <- 1L
  seed
}

#County-structural SMC variant (Supplementary Section 3 follow-up): instead of
#a soft county penalty, pass the county column to redist_smc's HARD `counties=`
#constraint (the sampler then produces only plans that split <= ndists-1
#counties), and keep a SOFT add_constr_splits penalty on municipalities only.
#This is the redist-intended mechanism for admin preservation; see
#PLAN_tighter_constraints.md / the Section 3 research note. Same isolated
#d{NN}_<dir_suffix>/ layout, PID-tagged writes, and pop_temper-only fallback
#ladder as ensure_smc_one_combined; pop_tol fixed by the caller (0.01).
ensure_smc_one_county_struct <- function(ndists,
                                         cousub_strength = 5,
                                         target_nsims = 2000L,
                                         batch_size   = 2000L,
                                         pop_tol      = 0.01,
                                         pop_temper   = 0,
                                         fallback_ladder = list(
                                           list(pop_temper = 0.01),
                                           list(pop_temper = 0.02),
                                           list(pop_temper = 0.03)),
                                         dir_suffix = "ccstruct",
                                         ncores = max(1L, parallel::detectCores() - 1L)){
  d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", ndists, dir_suffix))
  dir.create(d_dir, showWarnings = FALSE, recursive = TRUE)
  worker_pid <- Sys.getpid()
  list_batches <- function(){
    fs <- sort(list.files(d_dir, pattern = "^batch_.*\\.Rdata$", full.names = TRUE))
    fs[!grepl("\\.tmp$", fs)]
  }
  current <- length(list_batches()) * batch_size
  if(current >= target_nsims) return(invisible())
  if(!file.exists(dataset_cache)) stop("dataset cache missing — re-run stage 1.")
  load(dataset_cache)   # bgs
  bgs$seed <- .assign_bg_seed(bgs, ndists)
  pa_map <- build_pa_map_bg(bgs, ndists, seed_col = "seed", pop_tol = pop_tol)

  initial_attempt <- list(pop_temper = pop_temper, cousub_strength = cousub_strength)
  fill_defaults <- function(step){
    if(is.null(step$pop_temper))      step$pop_temper      <- initial_attempt$pop_temper
    if(is.null(step$cousub_strength)) step$cousub_strength <- initial_attempt$cousub_strength
    step
  }
  attempts <- c(list(initial_attempt), lapply(fallback_ladder, fill_defaults))

  worker_idx <- 0L
  while(current < target_nsims){
    worker_idx <- worker_idx + 1L
    this_nsims <- as.integer(min(batch_size, target_nsims - current))
    out_f      <- file.path(d_dir, sprintf("batch_pid%d_%03d.Rdata", worker_pid, worker_idx))
    batch_seed <- as.integer((as.numeric(Sys.time()) + worker_idx * 7919 +
                                worker_pid * 104729) %% .Machine$integer.max)
    plans_obj <- NULL; used_pop_temper <- NA_real_
    used_cousub_strength <- NA_real_; wall_min <- NA_real_
    for(k in seq_along(attempts)){
      a <- attempts[[k]]
      constr <- redist::redist_constr(pa_map)
      constr <- redist::add_constr_splits(constr, strength = a$cousub_strength,
                                          admin = cousub_id)
      message(sprintf("d=%d ccstruct pid=%d batch %03d step %d/%d (counties=county_fips[HARD] cousub_str=%g pop_temper=%g pop_tol=%g): redist_smc nsims=%d ncores=%d  (have %d/%d, seed=%d) ...",
                      ndists, worker_pid, worker_idx, k, length(attempts),
                      a$cousub_strength, a$pop_temper, pop_tol, this_nsims, ncores,
                      current, target_nsims, batch_seed))
      t0 <- Sys.time(); set.seed(batch_seed)
      result <- tryCatch(
        redist::redist_smc(pa_map, nsims = this_nsims, counties = county_fips,
                           compactness = 1, constraints = constr, ncores = ncores,
                           pop_temper = a$pop_temper),
        error = function(e){
          message("    SMC failed at step ", k, ": ", conditionMessage(e)); NULL })
      if(!is.null(result)){
        plans_obj <- result; used_pop_temper <- a$pop_temper
        used_cousub_strength <- a$cousub_strength
        wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        message(sprintf("  pid=%d batch %03d done in %.1f min (step %d/%d)",
                        worker_pid, worker_idx, wall_min, k, length(attempts)))
        break
      }
    }
    if(is.null(plans_obj))
      stop(sprintf("ensure_smc_one_county_struct: exhausted ladder for d=%d pid=%d batch %03d",
                   ndists, worker_pid, worker_idx))
    #County is a HARD structural constraint (counties=), so it has no strength.
    used_county_strength <- NA_real_
    tmp_f <- paste0(out_f, ".tmp")
    save(plans_obj, batch_seed, this_nsims, wall_min, used_pop_temper,
         used_cousub_strength, used_county_strength, pop_tol, worker_pid, worker_idx,
         file = tmp_f, compress = "xz")
    file.rename(tmp_f, out_f)
    message("  wrote ", out_f)
    rm(plans_obj); gc(verbose = FALSE)
    current <- length(list_batches()) * batch_size
  }
  invisible()
}

#Build (cache-gated) one county-structural cell across configs: SMC (hard
#county + soft cousub) into the variant dir, then metrics + Polsby-Popper.
ensure_county_struct_cell <- function(cousub_strength, variant,
                                      configs = c(3, 6, 9, 17, 18),
                                      target_nsims = 2000L, batch_size = 2000L,
                                      metrics, metric_vals, bgs, perim_df = NULL,
                                      ncores = max(1L, parallel::detectCores() - 1L)){
  if(is.null(perim_df)) perim_df <- redistmetrics::prep_perims(shp = bgs)
  for(nd in configs){
    ensure_smc_one_county_struct(
      ndists = nd, cousub_strength = cousub_strength,
      target_nsims = target_nsims, batch_size = batch_size,
      pop_tol = 0.01, pop_temper = 0, dir_suffix = variant, ncores = ncores)
    generate_metrics_bg(metrics, metric_vals, ndists = nd, bgs, admin_col = variant)
    generate_polsby_bg(ndists = nd, bgs.full = bgs, admin_col = variant,
                       perim_df = perim_df)
  }
  invisible()
}

#Per-district Dem-share distribution + seat-outcome class distribution across an
#SMC ensemble, feeding the d=9-anomaly figures (Supplementary Materials Section
#4: fig-supp-rank-shares). Reads the stored plan matrices directly from the
#variant's SMC batch dirs (bg_smc_dir/d{NN}_{admin_col}), so it is
#ensemble-agnostic: pass admin_col = "ccstruct_m5" (hard-county MAIN) or
#"county_cousub" (soft-penalty contrast). Returns
#list(rank_summary, outcome_summary, per_plan_outcomes). Lifted verbatim from
#the former standalone verify_anomaly_numbers.R so the consumed cache is now
#built inside the Data Prep pipeline like every other object the qmds load.
#`configs` is a list of list(nd=, M=); d=50 is intentionally excluded (the rank
#figure only covers 3/6/9/17/18).
build_anomaly_bg <- function(bgs, admin_col,
                             configs = list(list(nd =  3, M = 6),
                                            list(nd =  6, M = 3),
                                            list(nd =  9, M = 2),
                                            list(nd = 17, M = 1),
                                            list(nd = 18, M = 1))){
  biden <- as.numeric(bgs$sum_biden); biden[is.na(biden)] <- 0
  trump <- as.numeric(bgs$sum_trump); trump[is.na(trump)] <- 0
  twop  <- biden + trump
  group_pct_c <- get("group_pct", envir = asNamespace("redist"))

  list_batches <- function(nd){
    d_dir <- file.path(bg_smc_dir, sprintf("d%02d_%s", nd, admin_col))
    fs <- sort(list.files(d_dir, pattern = "^batch_.*\\.Rdata$", full.names = TRUE))
    fs[!grepl("\\.tmp$", fs)]
  }

  analyze_one <- function(nd, M){
    fs <- list_batches(nd)
    if(length(fs) == 0L){
      message(sprintf("build_anomaly_bg: d=%d (%s): no batches found", nd, admin_col))
      return(NULL)
    }
    message(sprintf("build_anomaly_bg: d=%d (M=%d, %s): %d batches", nd, M, admin_col, length(fs)))
    sorted_pieces  <- list()   # nd x nplans, sorted descending by column
    outcome_pieces <- list()   # (M+1) x nplans, row k = districts with (k-1) Dem seats
    total_plans    <- 0L
    for(i in seq_along(fs)){
      f <- fs[i]
      e <- new.env()
      ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
      if(!ok || !exists("plans_obj", envir = e)){ rm(e); next }
      pm <- redist::get_plans_matrix(e$plans_obj)
      pm <- pm[, -1, drop = FALSE]
      if(ncol(pm) == 0L){ rm(e); next }
      pct_dem <- group_pct_c(pm, biden, twop, as.integer(nd))
      # Sorted per-district shares (rank 1 = most Democratic in each plan).
      sorted_pct <- apply(pct_dem, 2, sort, decreasing = TRUE)
      if(!is.matrix(sorted_pct)) sorted_pct <- matrix(sorted_pct, nrow = nd)
      # Droop seat allocation: Dem seats per district = floor(pct_dem * (M+1)).
      dem_seats <- floor(pct_dem * (M + 1))
      outcome_counts <- apply(dem_seats, 2, function(v) tabulate(v + 1L, nbins = M + 1L))
      if(!is.matrix(outcome_counts)) outcome_counts <- matrix(outcome_counts, nrow = M + 1L)
      sorted_pieces[[length(sorted_pieces)+1L]]   <- sorted_pct
      outcome_pieces[[length(outcome_pieces)+1L]] <- outcome_counts
      total_plans <- total_plans + ncol(pm)
      rm(e, pm, pct_dem, sorted_pct, dem_seats, outcome_counts); gc(verbose = FALSE)
    }
    if(total_plans == 0L) return(NULL)
    sorted_all  <- do.call(cbind, sorted_pieces)
    outcome_all <- do.call(cbind, outcome_pieces)
    rank_summary <- tibble::tibble(
      nd = nd, M = M, rank = seq_len(nd),
      median_pct_dem = apply(sorted_all, 1, median),
      mean_pct_dem   = apply(sorted_all, 1, mean),
      q05_pct_dem    = apply(sorted_all, 1, quantile, probs = 0.05),
      q25_pct_dem    = apply(sorted_all, 1, quantile, probs = 0.25),
      q75_pct_dem    = apply(sorted_all, 1, quantile, probs = 0.75),
      q95_pct_dem    = apply(sorted_all, 1, quantile, probs = 0.95))
    d_seats_vec <- seq(0L, M); r_seats_vec <- seq(M, 0L, by = -1L)
    outcome_summary <- tibble::tibble(
      nd = nd, members = M, d_seats = d_seats_vec, r_seats = r_seats_vec,
      class  = sprintf("%dD-%dR", d_seats_vec, r_seats_vec),
      median = apply(outcome_all, 1, median),
      mean   = apply(outcome_all, 1, mean),
      q05    = apply(outcome_all, 1, quantile, probs = 0.05),
      q95    = apply(outcome_all, 1, quantile, probs = 0.95),
      share_of_plans_any = rowMeans(outcome_all > 0))
    per_plan <- tibble::as_tibble(t(outcome_all), .name_repair = "minimal")
    colnames(per_plan) <- sprintf("n_%dD%dR", seq(0, M), seq(M, 0, by = -1))
    per_plan$nd <- nd; per_plan$M <- M
    list(rank_summary = rank_summary, outcome_summary = outcome_summary,
         per_plan = per_plan, total_plans = total_plans)
  }

  results <- lapply(configs, function(cfg) analyze_one(cfg$nd, cfg$M))
  list(
    rank_summary      = dplyr::bind_rows(lapply(results, `[[`, "rank_summary")),
    outcome_summary   = dplyr::bind_rows(lapply(results, `[[`, "outcome_summary")),
    per_plan_outcomes = dplyr::bind_rows(lapply(results, `[[`, "per_plan"))
  )
}

#Assemble a tidy summary across penalty cells × configs for the Supplementary
#Section 3 figure/table. `cells` is a list of
#  list(county=, cousub=, variant=, label=)
#Returns one row per (cell, config, metric) with ensemble q05/median/q95/min/max
#plus the matching enacted value (NA where no enacted plan exists).
assemble_penalty_sweep <- function(cells, configs = c(17, 18)){
  metric_specs <- list(
    list(name = "County splits",       col = "n_county_splits", src = "met", enacted = "n_county"),
    list(name = "Municipality splits", col = "n_muni_splits",   src = "met", enacted = "n_muni"),
    list(name = "Mean Polsby-Popper",  col = "mean_polsby_popper", src = "pp", enacted = "pp_mean"),
    list(name = "Dem seat share",      col = "ShareDemSeats",   src = "met", enacted = "dem_share")
  )
  rows <- list()
  for(cell in cells){
    for(nd in configs){
      met <- read_variant_metrics(cell$variant, nd)
      pp  <- read_variant_polsby (cell$variant, nd)
      if(is.null(met)) next
      enacted <- .enacted_for(nd)
      for(ms in metric_specs){
        vec <- if(ms$src == "pp"){ if(!is.null(pp)) pp[[ms$col]] else NULL } else met[[ms$col]]
        if(is.null(vec) || !length(vec)) next
        v <- vec[is.finite(vec)]
        q <- quantile(v, c(0.05, 0.5, 0.95), names = FALSE, na.rm = TRUE)
        rows[[length(rows)+1L]] <- data.frame(
          variant = cell$variant, label = cell$label,
          county_strength = cell$county, cousub_strength = cell$cousub,
          ndists = nd, n = length(v), metric = ms$name,
          q05 = q[1], median = q[2], q95 = q[3], min = min(v), max = max(v),
          enacted = if(!is.null(enacted)) enacted[[ms$enacted]] else NA_real_,
          stringsAsFactors = FALSE)
      }
    }
  }
  dplyr::bind_rows(rows)
}

#Console diagnostic: where does each enacted value sit within the variant
#ensemble (gate = within central `central` on county splits, muni splits, PP)?
penalty_match_report <- function(variant, ndists, central = 0.90, markers = TRUE){
  met <- read_variant_metrics(variant, ndists)
  pp  <- read_variant_polsby (variant, ndists)
  if(is.null(met)){
    cat(sprintf("[match] d=%d variant=%s: no metric caches found.\n", ndists, variant))
    return(invisible(NULL))
  }
  enacted <- .enacted_for(ndists)
  lo <- (1 - central) / 2; hi <- 1 - lo
  specs <- list(
    list(name = "n_county_splits",    vec = met$n_county_splits,
         enacted = if(!is.null(enacted)) enacted$n_county else NA_real_, gate = TRUE),
    list(name = "n_muni_splits",      vec = met$n_muni_splits,
         enacted = if(!is.null(enacted)) enacted$n_muni   else NA_real_, gate = TRUE),
    list(name = "mean_polsby_popper", vec = if(!is.null(pp)) pp$mean_polsby_popper else NULL,
         enacted = if(!is.null(enacted)) enacted$pp_mean  else NA_real_, gate = TRUE),
    list(name = "ShareDemSeats",      vec = met$ShareDemSeats,
         enacted = if(!is.null(enacted)) enacted$dem_share else NA_real_, gate = FALSE)
  )
  rows <- list(); gate_pass <- logical(0)
  for(s in specs){
    if(is.null(s$vec) || !length(s$vec)) next
    v <- s$vec[is.finite(s$vec)]
    q <- quantile(v, c(lo, 0.5, hi), names = FALSE, na.rm = TRUE)
    inside <- if(is.na(s$enacted)) NA else (s$enacted >= q[1] & s$enacted <= q[3])
    if(isTRUE(s$gate) && !is.na(inside)) gate_pass <- c(gate_pass, inside)
    rows[[length(rows)+1L]] <- data.frame(metric = s$name,
      ens_q_lo = q[1], ens_median = q[2], ens_q_hi = q[3], enacted = s$enacted,
      enacted_pctl = if(is.na(s$enacted)) NA_real_ else mean(v <= s$enacted) * 100,
      inside = inside, gated = s$gate, stringsAsFactors = FALSE)
  }
  tab <- dplyr::bind_rows(rows)
  cat(sprintf("\n=== Penalty match: d=%d  variant=%s  (n=%d plans, central=%.0f%%) ===\n",
              ndists, variant, length(met$ShareDemSeats), 100 * central))
  if(is.null(enacted)) cat("  (no enacted plan at this magnitude — distribution only, no gate)\n")
  cat(sprintf("%-20s %10s %10s %10s %10s %8s %7s\n", "metric",
              sprintf("q%02.0f", 100*lo), "median", sprintf("q%02.0f", 100*hi),
              "enacted", "pctl", "inside"))
  cat(strrep("-", 80), "\n", sep = "")
  for(i in seq_len(nrow(tab)))
    cat(sprintf("%-20s %10.4g %10.4g %10.4g %10s %7s %7s\n",
                tab$metric[i], tab$ens_q_lo[i], tab$ens_median[i], tab$ens_q_hi[i],
                ifelse(is.na(tab$enacted[i]), "—", formatC(tab$enacted[i], format="g", digits=4)),
                ifelse(is.na(tab$enacted_pctl[i]), "—", sprintf("%.0f", tab$enacted_pctl[i])),
                ifelse(is.na(tab$inside[i]), "—", ifelse(tab$inside[i], "YES", "no"))))
  if(length(gate_pass))
    cat(sprintf("\n  GATE (county splits & muni splits & PP within central %.0f%%): %s\n",
                100*central, if(all(gate_pass)) "PASS" else "FAIL"))
  if(isTRUE(markers)){
    mk <- tryCatch(read_variant_smc_markers(variant, ndists), error = function(e) NULL)
    if(!is.null(mk) && nrow(mk)){
      cat(sprintf("\n  capacity-to-solve: %d batch(es), %s sims, wall %.1f min total (%.1f/batch)\n",
                  nrow(mk), formatC(sum(mk$this_nsims, na.rm=TRUE), big.mark=",", format="d"),
                  sum(mk$wall_min, na.rm=TRUE), mean(mk$wall_min, na.rm=TRUE)))
      cat(sprintf("    pop_temper used: max=%.3g (any > 0 means the fallback ladder fired)\n",
                  max(mk$used_pop_temper, na.rm = TRUE)))
      cat(sprintf("    strengths used : county=%s  cousub=%s\n",
                  paste(unique(mk$used_cnty_str), collapse="/"),
                  paste(unique(mk$used_cousub_str), collapse="/")))
    }
  }
  invisible(tab)
}

.ks_d <- function(x, y){
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if(length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(as.numeric(stats::ks.test(x, y)$statistic))
}

#Console diagnostic: per-batch convergence (KS_D of each new batch vs the
#cumulative-prior distribution + normalized mean shift) for a variant/config,
#plus an optional tail-vs-shortburst gap. Appends a per-(variant,config) CSV.
variant_convergence_report <- function(variant, ndists, sb_overlay = NULL,
                                       ks_thresh = 0.03, delta_thresh = 0.05){
  batches <- read_variant_metric_batches(variant, ndists)
  if(!length(batches)){
    cat(sprintf("[conv] d=%d variant=%s: no metric caches found.\n", ndists, variant))
    return(invisible(NULL))
  }
  sizes <- vapply(batches, nrow, integer(1))
  cat(sprintf("\n=== Convergence: d=%d  variant=%s  %d batch(es), %d plans ===\n",
              ndists, variant, length(batches), sum(sizes)))
  for(i in seq_along(batches))
    cat(sprintf("  %-26s n=%6d  cum_n=%7d\n", names(batches)[i], sizes[i], sum(sizes[1:i])))
  m_use <- intersect(VARIANT_METRICS, colnames(batches[[1]]))
  if(length(batches) < 2) cat("\n  Need >=2 batches for convergence shifts. (Reporting tail only.)\n")
  rows <- list()
  if(length(batches) >= 2){
    for(i in 2:length(batches)){
      cum_prev <- dplyr::bind_rows(batches[1:(i-1)]); new_b <- batches[[i]]
      for(m in m_use){
        xp <- cum_prev[[m]]; xn <- new_b[[m]]
        if(anyNA(xp) || anyNA(xn)) next
        sdp <- sd(xp); if(!is.finite(sdp) || sdp == 0) next
        rows[[length(rows)+1L]] <- data.frame(batch_idx = i, batch_tag = names(batches)[i],
          cum_n_prev = nrow(cum_prev), cum_n_now = nrow(cum_prev) + nrow(new_b),
          metric = m, cum_mean_prev = mean(xp), batch_mean = mean(xn),
          delta_batch_norm = (mean(xn) - mean(xp)) / sdp,
          ks_d_batch = .ks_d(xp, xn), stringsAsFactors = FALSE)
      }
    }
  }
  shifts <- dplyr::bind_rows(rows)
  if(nrow(shifts)){
    last <- shifts[shifts$batch_idx == max(shifts$batch_idx), ]
    cat(sprintf("\n  --- latest batch (%s) shifts ---\n", last$batch_tag[1]))
    cat(sprintf("  %-26s %12s %14s %12s\n", "metric", "batch_mean", "delta_batch/sd", "KS_D(batch)"))
    cat(strrep("-", 68), "\n", sep = "")
    for(i in seq_len(nrow(last)))
      cat(sprintf("  %-26s %12.4g %+14.4f %12.4f\n", last$metric[i], last$batch_mean[i],
                  last$delta_batch_norm[i], last$ks_d_batch[i]))
    converged <- all(abs(last$delta_batch_norm) <= delta_thresh, na.rm = TRUE) &&
                 all(last$ks_d_batch <= ks_thresh, na.rm = TRUE)
    cat(sprintf("\n  VERDICT: %s (KS_D <= %.3g and |delta/sd| <= %.3g on all metrics)\n",
                if(converged) "CONVERGED — more sims unlikely to move the bulk"
                else          "NOT YET — distribution still shifting; more sims may help",
                ks_thresh, delta_thresh))
  }
  if(!is.null(sb_overlay)){
    cum_all <- dplyr::bind_rows(batches)
    sb_cfg  <- sb_overlay[sb_overlay$ensemble == "county_cousub" &
                          sb_overlay$ndists == ndists, , drop = FALSE]
    tl <- list()
    for(t in .TAIL_MAP){
      if(!t$metric %in% colnames(cum_all)) next
      sb_rows <- sb_cfg[sb_cfg$target == t$target, , drop = FALSE]
      if(!nrow(sb_rows) || !t$metric %in% colnames(sb_rows)) next
      tl[[length(tl)+1L]] <- data.frame(metric = t$metric,
        ensemble_max = max(cum_all[[t$metric]], na.rm = TRUE),
        shortburst_opt = max(sb_rows[[t$metric]], na.rm = TRUE), stringsAsFactors = FALSE)
    }
    tl <- dplyr::bind_rows(tl)
    if(nrow(tl)){
      tl$gap <- tl$shortburst_opt - tl$ensemble_max
      cat("\n  --- tail vs shortburst (cumulative ensemble) ---\n")
      cat("  (shortburst optima are from the MAIN constraints — an approximate ceiling)\n")
      cat(sprintf("  %-26s %12s %14s %10s\n", "metric", "ensemble_max", "shortburst_opt", "gap"))
      cat(strrep("-", 66), "\n", sep = "")
      for(i in seq_len(nrow(tl)))
        cat(sprintf("  %-26s %12.4g %14.4g %10.4g\n",
                    tl$metric[i], tl$ensemble_max[i], tl$shortburst_opt[i], tl$gap[i]))
    }
  }
  if(nrow(shifts)){
    shifts$variant <- variant; shifts$ndists <- ndists
    shifts$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_path <- file.path(bg_root, sprintf("convergence_%s_d%02d.csv", variant, ndists))
    write.table(shifts, file = log_path, sep = ",", row.names = FALSE,
                col.names = !file.exists(log_path), append = file.exists(log_path), quote = TRUE)
    cat(sprintf("\n  appended %d rows to %s\n", nrow(shifts), log_path))
  }
  invisible(shifts)
}

#Per-BG competitiveness
#  Marginal_{NN}: fraction of plans where the BG's containing district is
#                 "marginal" (dist_d < threshold).#
#  FlipProb_{NN}: expected probability that the BG's containing district
#                 flips a seat under a single inter-election swing modeled
#                 as Normal(0, sigma). 
build_marginal_bg <- function(bgs.full,
                              configs   = list(
                                list(nd =  3, M = 6),
                                list(nd =  6, M = 3),
                                list(nd =  9, M = 2),
                                list(nd = 17, M = 1),
                                list(nd = 18, M = 1),
                                list(nd = 50, M = 1)),
                              admin_col = "county_cousub",
                              threshold = 0.05,
                              sigma     = 0.03){
  group_pct_c <- get("group_pct", envir = asNamespace("redist"))
  biden <- as.numeric(bgs.full$sum_biden); biden[is.na(biden)] <- 0
  trump <- as.numeric(bgs.full$sum_trump); trump[is.na(trump)] <- 0
  twop  <- biden + trump

  suffix <- if(is.null(admin_col)) "" else paste0("_", admin_col)

  compute_one <- function(nd, M){
    d_dir <- file.path(bg_smc_dir, sprintf("d%02d%s", nd, suffix))
    files <- sort(list.files(d_dir, pattern = "^batch_.*\\.Rdata$",
                              full.names = TRUE))
    files <- files[!grepl("\\.tmp$", files)]
    if(length(files) == 0L) return(NULL)
    q       <- 1 / (M + 1)
    bg_marg <- rep(0, nrow(bgs.full))
    bg_flip <- rep(0, nrow(bgs.full))
    n_plans <- 0L
    for(f in files){
      e <- new.env()
      ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
      if(!ok || !exists("plans_obj", envir = e)){ rm(e); next }
      pm <- redist::get_plans_matrix(e$plans_obj)
      pm <- pm[, -1, drop = FALSE]
      nplans <- ncol(pm)
      pct_dem <- group_pct_c(pm, biden, twop, as.integer(nd))
      v_mod   <- pct_dem %% q
      dist_d  <- pmin(v_mod, q - v_mod)
      marg_d  <- (dist_d < threshold)
      flip_d  <- flip_prob_q(dist_d, q, sigma)
      for(p in seq_len(nplans)){
        bg_marg <- bg_marg + as.numeric(marg_d[pm[, p], p])
        bg_flip <- bg_flip + flip_d[pm[, p], p]
      }
      n_plans <- n_plans + nplans
      rm(e, pm, pct_dem, v_mod, dist_d, marg_d, flip_d); gc(verbose = FALSE)
    }
    list(frac    = bg_marg / n_plans,
         flip    = bg_flip / n_plans,
         n_plans = n_plans)
  }

  marginal_bg <- bgs.full[, "GEOID"]
  for(cfg in configs){
    message(sprintf("build_marginal_bg: d=%02d (M=%d) ...", cfg$nd, cfg$M))
    r <- compute_one(cfg$nd, cfg$M)
    if(!is.null(r)){
      marginal_bg[[sprintf("Marginal_%d", cfg$nd)]] <- r$frac
      marginal_bg[[sprintf("FlipProb_%d", cfg$nd)]] <- r$flip
      message(sprintf("  Marginal median=%.3f q95=%.3f | FlipProb median=%.3f q95=%.3f | n_plans=%d",
                      median(r$frac, na.rm = TRUE),
                      quantile(r$frac, 0.95, names = FALSE, na.rm = TRUE),
                      median(r$flip, na.rm = TRUE),
                      quantile(r$flip, 0.95, names = FALSE, na.rm = TRUE),
                      r$n_plans))
    }
  }
  marginal_bg
}

#Shortburst overlay table
#Feeds the shortburst-point overlay in fig-eg-wv and fig-bw-bs.
score_shortburst_optima <- function(bgs){
  rows <- list()
  for(spec in district_specs){
    nd <- spec$ndists
    nm_label <- c("3" = "Three Districts", "6" = "Six Districts",
                  "9" = "Nine Districts",  "17"= "Seventeen Districts",
                  "18"= "Eighteen Districts","50" = "Fifty Districts")[as.character(nd)]
    #county_cousub = OLD main (retained for Supplementary); ccstruct_m5 = NEW
    #main (Beyond_the_Border_Wars.qmd overlays). Both are scored when their shortburst cache
    #exists; configs/ensembles with no cache are silently skipped.
    for(ensemble in c("unconstrained", "county_cousub", "ccstruct_m5")){
      suffix <- if(identical(ensemble, "unconstrained")) "" else paste0("_", ensemble)
      cache_f <- file.path(bg_sb_dir, sprintf("shortburst_bg_d%02d%s.Rdata",
                                              nd, suffix))
      if(!file.exists(cache_f)) next
      load(cache_f)   # populates `out`

      for(run_name in names(out)){
        sb <- out[[run_name]]
        #Parse run_name. Two formats coexist:
        #  "<target>_k<N>"  e.g. "dem_max_k1"          -- k-indexed sweep
        #  "<target>"       e.g. "black_seats_max"     -- scalar objective
        #Detect by trailing _k<digits>.
        if(grepl("_k[0-9]+$", run_name)){
          target <- sub("_k[0-9]+$", "", run_name)
          k_val  <- as.integer(sub(".*_k", "", run_name))
        }else{
          target <- run_name
          k_val  <- NA_integer_
        }

        #Pick the optimum plan. New list format vs legacy redist_plans.
        plan <- if(is.list(sb) && "plan" %in% names(sb)){
          sb$plan
        }else{
          pm <- .plans_matrix(sb); pm[, ncol(pm)]
        }

        df <- st_drop_geometry(bgs)
        df$District <- plan
        df$Voters   <- df$sum_biden + df$sum_trump

        num_members <- 18/nd; if(num_members < 2) num_members <- 1

        #Wasted votes / Dem seats
        sb_d  <- sapply(seq_len(nd), \(d) sum(df$sum_biden[df$District==d], na.rm=TRUE))
        vt_d  <- sapply(seq_len(nd), \(d) sum(df$Voters   [df$District==d], na.rm=TRUE))
        toWin <- vt_d / (num_members + 1)
        dem_w <- floor(sb_d / toWin);  rep_w <- floor((vt_d - sb_d) / toWin)
        DemSeats      <- sum(dem_w);   RepSeats <- sum(rep_w)
        ShareDemSeats <- DemSeats / max(DemSeats + RepSeats, 1)
        dem_wasted <- sb_d - dem_w * toWin
        rep_wasted <- (vt_d - sb_d) - rep_w * toWin
        WastedVotes <- (sum(dem_wasted) - sum(rep_wasted)) / sum(vt_d)

        #Black metrics
        pb_d   <- sapply(seq_len(nd), \(d) sum(df$sum_BVAP[df$District==d], na.rm=TRUE))
        pvap_d <- sapply(seq_len(nd), \(d) sum(df$sum_VAP [df$District==d], na.rm=TRUE))
        bk_toWin  <- pvap_d / (num_members + 1)
        bk_toWin2 <- pvap_d / (num_members + 2)
        black_w   <- floor(pb_d / bk_toWin)
        black_w37 <- floor(pb_d / bk_toWin2)
        other_w   <- floor((pvap_d - pb_d) / bk_toWin)
        BlackSeats     <- sum(black_w);   BlackSeats37 <- sum(black_w37)
        bk_wasted      <- pb_d - black_w * bk_toWin
        ot_wasted      <- (pvap_d - pb_d) - other_w * bk_toWin
        BlackWasted    <- (sum(bk_wasted) - sum(ot_wasted)) / sum(pvap_d)
        BlackPctWasted <- (sum(bk_wasted) / sum(pb_d)) -
          (sum(ot_wasted) / sum(pvap_d - pb_d))

        ShareRepSeats <- RepSeats / max(DemSeats + RepSeats, 1)

        #Competitiveness summaries — mirror n_marginal_5pt /
        #prob_seat_change_sigma3 from generate_metrics_bg so the fig-cd-box
        #overlay reads from the same algebra as the boxplots.
        dp_d                    <- sb_d / vt_d
        q_drp                   <- 1 / (num_members + 1)
        v_modq_d                <- dp_d %% q_drp
        delta_d                 <- pmin(v_modq_d, q_drp - v_modq_d)
        n_marginal_5pt          <- sum(delta_d < 0.05)
        share_marginal          <- n_marginal_5pt / nd
        prob_seat_change_sigma3 <- mean(flip_prob_q(delta_d, q_drp, 0.03))

        rows[[length(rows)+1L]] <- data.frame(
          Type            = nm_label,
          ndists          = nd,
          ensemble        = ensemble,
          target          = target,
          k               = k_val,
          DemSeats        = DemSeats,
          RepSeats        = RepSeats,
          ShareDemSeats   = ShareDemSeats,
          ShareRepSeats   = ShareRepSeats,
          WastedVotes     = WastedVotes,
          BlackSeats      = BlackSeats,
          BlackSeats37    = BlackSeats37,
          BlackWasted     = BlackWasted,
          BlackPctWasted  = BlackPctWasted,
          n_marginal_5pt          = n_marginal_5pt,
          share_marginal          = share_marginal,
          prob_seat_change_sigma3 = prob_seat_change_sigma3
        )
      }
    }
  }

  out_df <- bind_rows(rows)
  if(nrow(out_df) == 0L) return(out_df)
  out_df$Type <- factor(out_df$Type,
                  levels = c("Three Districts","Six Districts","Nine Districts",
                  "Seventeen Districts","Eighteen Districts","Fifty Districts"))
  out_df
}

#quick function to replace space with newline in labels
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
#Strip the embedded newline that addline_format() 
strip_type <- function(d){ d$Type <- factor(gsub("\n", " ", as.character(d$Type)),levels = c("Three Districts","Six Districts", "Nine Districts","Eighteen Districts", "Fifty Districts")); d }
#Display a cached figure: during knit
show_cached_fig <- function(path, preview_max_px = 800){
  if(isTRUE(getOption("knitr.in.progress"))){
    knitr::include_graphics(path)
  }else{
    img <- png::readPNG(path)
    s   <- max(1L, as.integer(ceiling(max(dim(img)[1:2]) / preview_max_px)))
    if(s > 1L)
      img <- img[seq(1L, nrow(img), by = s),
                 seq(1L, ncol(img), by = s), , drop = FALSE]
    grid::grid.newpage()
    grid::grid.raster(img)
  }
}

#Helper for jaccard
ratio   <- function(a, b){ d <- pmax(a, b); ifelse(d == 0, NA_real_, pmin(a, b) / d) }
#calculate jaccard similarity
jaccard <- function(a, b){
  ok <- !is.na(a) & !is.na(b)
  num <- sum(pmin(a[ok], b[ok])); den <- sum(pmax(a[ok], b[ok]))
  if(den == 0) NA_real_ else num / den
}