####################################################
### Data Prep for Beyond the Border Wars — Block Group version
### 5/12/26
####################################################

#Libraries
require(reshape2)
require(redist)
require(redistmetrics)
require(igraph)
require(patchwork)
require(parallel)

source("BBW - Custom Functions.R")

##
#### Setup section.
##

#R Global variable set
overwrite = FALSE
#Pipeline verbosity. Set FALSE *before* sourcing this script (e.g. in
#Beyond_the_Border_Wars.qmd's setup chunk: `verbose <- FALSE` then source(...)) to
#silence the per-batch progress messages from the SMC, shortburst, and
#metric helpers. Errors and warnings always print regardless.
if(!exists("verbose")) verbose <- TRUE
maybe_silent <- function(expr){
  if(isTRUE(verbose)) expr else suppressMessages(expr)
}
options(tigris_use_cache = TRUE)
pa_crs <- "+proj=lcc +lat_1=40.88333333333333 +lat_2=41.95 +lat_0=40.16666666666666 +lon_0=-77.75 +x_0=600000.0000000001 +y_0=0 +datum=NAD83 +units=us-ft +no_defs"

#Build directory structure to support simulation and shortburst outputs
bg_root     <- "./Intermediate Data/bgs"
bg_smc_dir  <- file.path(bg_root, "smc_plans")
bg_sb_dir   <- file.path(bg_root, "shortburst")
for (d in c(bg_root, bg_smc_dir, bg_sb_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

#simulation input files
dataset_cache <- file.path(bg_root, "bgs.dataset.Rdata")


#District configurations 
district_specs <- list(
  list(ndists =  3, num_members = 6),
  list(ndists =  6, num_members = 3),
  list(ndists =  9, num_members = 2),
  list(ndists = 17, num_members = 1),
  list(ndists = 18, num_members = 1),
  list(ndists = 50, num_members = 1)
)

# Block-group dataset assembly
# End result is a set of 2020 PA block groups with the
# sum_TotalPop, sum_biden, sum_trump, sum_VAP, sum_BVAP
# plus ACS demographics (NHBlack, TotalPop_ACS, POP65Plus, Households, HHwithChld,
# InPoverty, PovDeterm, EmpTotal, EmpManufac) and existing-plan identifiers
# (RoughCD118, RoughSLDU, RoughCD116).
if(!file.exists(dataset_cache) | overwrite==TRUE){

  #Per-stage caches live under bg_root. Get deleted if dataset is successfully created
  blocks_raw_cache   <- file.path(bg_root, "stage_blocks_census.Rdata")
  blocks_votes_cache <- file.path(bg_root, "stage_blocks_votes.Rdata")
  bgs_acs_cache      <- file.path(bg_root, "stage_bgs_acs.Rdata")
  muni_cache         <- file.path(bg_root, "pa_cousub_2020.Rdata")
  per_county_dir     <- file.path(bg_root, "blocks_by_county")
  dir.create(per_county_dir, showWarnings = FALSE, recursive = TRUE)

  #Stage 1: download 2020 census blocks via tidycensus 
  #Per-county download with caching so a partial run resumes cleanly.
    if(!file.exists(blocks_raw_cache)){
      pa_county_fips <- fips_codes %>% filter(state == "PA") %>% pull(county_code)
      block_list <- lapply(pa_county_fips, function(cty){
        f <- file.path(per_county_dir, sprintf("blocks_42%s.rds", cty))
        if(file.exists(f)){
          readRDS(f)
        }else{
          out <- get_decennial(
            geography = "block", sumfile = "pl", output = "wide",
            variables = c(TotalPop = "P1_001N",
                          VAP      = "P3_001N",
                          BVAP     = "P3_004N"), #single race, no ethnicity distinction. Could also use P4_006N which is the Not Hispanic equivalent 
            year = 2020, state = "PA", county = cty,
            geometry = TRUE, keep_geo_vars = TRUE)
          saveRDS(out, f)
          out
        }
      })
      blocks_raw <- do.call(rbind, block_list)
      blocks_raw <- blocks_raw[, c("GEOID","TRACTCE20","TotalPop","VAP","BVAP","geometry")]
      blocks_raw <- st_transform(blocks_raw, pa_crs) %>% st_make_valid()
      save(blocks_raw, file = blocks_raw_cache, compress = "xz")
    }else{
      load(blocks_raw_cache)
    }

    #Stage 2: load VTDs and apportion votes to blocks via tcw_blocks()
    if(!file.exists(blocks_votes_cache)){
      vtd <- read_sf(dsn = "./Input Data", layer = "pa_vtds20") %>%
        st_transform(pa_crs) %>%
        st_make_valid()
      vtd <- vtd[, c("GEOID20","trump","biden")]

      blocks_sf <- tcw_blocks(vtd, blocks_raw)
      save(blocks_sf, file = blocks_votes_cache, compress = "xz")
    }else{
      load(blocks_votes_cache)
    }

    #Stage 3: aggregate block-level attributes to block-group level
    bg_sum_cols <- c("sum_biden","sum_trump","sum_TotalPop","sum_VAP","sum_BVAP")
    bgs_attrs <- st_drop_geometry(blocks_sf) %>%
      mutate(GEOID = substr(GEOID, 1, 12)) %>%
      group_by(GEOID) %>%
      summarise(across(all_of(bg_sum_cols), \(x) sum(x, na.rm = TRUE)),
                .groups = "drop")

    #Stage 4: download BG-level ACS data with official BG geometry; join the
    #stage-3 attributes onto it.
    if(!file.exists(bgs_acs_cache)){
      #BG-level ACS. Poverty pulled from C17002 (income-to-poverty ratio)
      #rather than B17001 (poverty by sex by age); B17001 is suppressed at
      #BG geography
      b01001_65plus <- paste0("B01001_0", c(20:25, 44:49))
      bg_vars <- c(PovDeterm    = "C17002_001",
                   InPovBelow50 = "C17002_002",
                   InPov50to99  = "C17002_003",
                   Households   = "B09019_002",
                   HHwithChld   = "B09019_008",
                   NHBlack      = "B03002_004",
                   TotalPop_ACS = "B03002_001",
                   setNames(b01001_65plus, paste0("Age_", b01001_65plus)))

      bgs_acs <- get_acs(geography="block group", variables=bg_vars,
                         state="PA", geometry=TRUE, year=2020, output="wide") %>%
        st_transform(pa_crs) %>%
        st_make_valid() %>%
        transmute(GEOID,
                  PovDeterm    = PovDetermE,
                  InPoverty    = InPovBelow50E + InPov50to99E,  # below 1.00 of poverty line
                  Households   = HouseholdsE,
                  HHwithChld   = HHwithChldE,
                  NHBlack      = NHBlackE,
                  TotalPop_ACS = TotalPop_ACSE,
                  POP65Plus    = rowSums(across(starts_with("Age_") &
                                                ends_with("E")), na.rm=TRUE))

      #Employment / manufacturing (C24030) only available at tract level
      tracts_emp <- get_acs(geography="tract",
                            variables=c("C24030_001","C24030_007","C24030_034"),
                            state="PA", geometry=FALSE, year=2020)
      tracts_emp$moe <- NULL
      tracts_emp$NAME <- NULL
      tracts_emp <- spread(tracts_emp, key="variable", value="estimate", fill=NA, convert=FALSE)
      colnames(tracts_emp) <- c("GEOID11","EmpTotal","EmpMale","EmpFemale")
      tracts_emp$EmpManufac <- tracts_emp$EmpMale + tracts_emp$EmpFemale
      tracts_emp <- tracts_emp[, c("GEOID11","EmpTotal","EmpManufac")]

      save(bgs_acs, tracts_emp, file=bgs_acs_cache, compress="xz")
    }else{
      load(bgs_acs_cache)
    }

    #Join stage-3 BG attributes onto the BG geometry/ACS sf
    bgs <- bgs_acs %>%
      left_join(bgs_attrs, by="GEOID") %>%
      mutate(GEOID11 = substr(GEOID, 1, 11))

    #Allocate tract-level employment counts to BGs by population share
    bg_pop_in_tract <- st_drop_geometry(bgs) %>%
      group_by(GEOID11) %>%
      summarise(tract_bg_pop = sum(sum_TotalPop, na.rm=TRUE), .groups="drop")
    bgs <- bgs %>%
      left_join(bg_pop_in_tract, by="GEOID11") %>%
      left_join(tracts_emp,      by="GEOID11") %>%
      mutate(pop_share  = ifelse(tract_bg_pop == 0, 0, sum_TotalPop / tract_bg_pop),
             EmpTotal   = EmpTotal   * pop_share,
             EmpManufac = EmpManufac * pop_share) %>%
      select(-tract_bg_pop, -pop_share, -GEOID11)

    #Stage 5: admin assignments at BG level
    bgs$county_fips <- substr(bgs$GEOID, 1, 5)
    if(!file.exists(muni_cache)){
      url  <- "https://www2.census.gov/geo/tiger/TIGER2020/COUSUB/tl_2020_42_cousub.zip"
      dest <- file.path(tempdir(), "tl_2020_42_cousub.zip")
      download.file(url, dest, mode = "wb")
      unzip(dest, exdir = file.path(tempdir(), "cousub"))
      cousub <- read_sf(file.path(tempdir(), "cousub", "tl_2020_42_cousub.shp")) %>%
        st_transform(pa_crs) %>%
        transmute(cousub_id = GEOID)
      save(cousub, file = muni_cache, compress = "xz")
    }else{
      load(muni_cache)
    }
    #Centroid-in-polygon assignment of BGs to cousubs. Falls back to
    #nearest cousub when a BG's centroid lies outside every polygon —
    #this happens for water-dominated BGs whose centroid sits in a river
    #or lake.
    bgs_centroid <- st_centroid(bgs)
    cousub_join  <- st_intersects(bgs_centroid, cousub)
    no_match     <- lengths(cousub_join) == 0
    if(any(no_match)){
      nearest_cousub <- st_nearest_feature(bgs_centroid[no_match, ], cousub)
      cousub_join[no_match] <- as.list(nearest_cousub)
      message(sprintf("    fell back to nearest cousub for %d BG(s) with water-area centroids",
                      sum(no_match)))
    }
    bgs$cousub_id <- cousub$cousub_id[vapply(cousub_join, `[`, integer(1), 1L)]

    for(lyr in list(
      list(col = "RoughCD118", url = "https://www2.census.gov/geo/tiger/TIGER_RD18/STATE/42_PENNSYLVANIA/42/tl_rd22_42_cd118.zip"),
      list(col = "RoughSLDU",  url = "https://www2.census.gov/geo/tiger/TIGER_RD18/STATE/42_PENNSYLVANIA/42/tl_rd22_42_sldu.zip"),
      list(col = "RoughCD116", url = "https://www2.census.gov/geo/tiger/TIGER2020PL/STATE/42_PENNSYLVANIA/42/tl_2020_42_cd116.zip")
    )){
      dest <- file.path(tempdir(), basename(lyr$url))
      if(!file.exists(dest)) download.file(lyr$url, dest, mode = "wb")
      layer_dir <- file.path(tempdir(), tools::file_path_sans_ext(basename(lyr$url)))
      unzip(dest, exdir = layer_dir)
      shp_path <- list.files(layer_dir, pattern = "\\.shp$", full.names = TRUE)[1]
      shp <- read_sf(shp_path) %>% st_transform(pa_crs)
      hits <- st_intersects(bgs_centroid, shp)
      bgs[[lyr$col]] <- vapply(hits,
        \(x) if(length(x)) as.integer(x[1]) else NA_integer_,
        integer(1))
      message("    joined ", lyr$col)
    }

    #Zero-fill any NA cells in numeric columns picked up from the joins.
    fill_cols <- c("sum_biden","sum_trump",
                   "TotalPop_ACS","NHBlack","PovDeterm","InPoverty",
                   "Households","HHwithChld","POP65Plus",
                   "EmpTotal","EmpManufac")
    fill_cols <- intersect(fill_cols, colnames(bgs))
    bgs <- bgs %>% mutate(across(all_of(fill_cols), \(x) replace_na(x, 0)))

  #Drop unpopulated BGs and any isolated singletons (PA has 1 disconnected
  #populated BG; downstream redist_smc needs a single connected component).
  #To avoid NA holes in BG-resolution maps, absorb each dropped BG's
  #geometry into its spatially nearest kept neighbor before removing the
  #row. The absorbing BG retains its own attribute values (dropped BGs
  #are unpopulated or singleton and contribute no analytic content).
  pop_idx     <- which(bgs$sum_TotalPop > 0)
  adj         <- redist::redist.adjacency(bgs[pop_idx, ])
  g           <- igraph::graph_from_adj_list(lapply(adj, `+`, 1L), mode = "all")
  comp        <- igraph::components(g)
  in_largest  <- comp$membership == which.max(comp$csize)
  keep_geoids <- bgs$GEOID[pop_idx][in_largest]
  keep_rows   <- which(bgs$GEOID %in% keep_geoids)
  drop_rows   <- setdiff(seq_len(nrow(bgs)), keep_rows)
  if(length(drop_rows) > 0){
    nearest_kept <- st_nearest_feature(bgs[drop_rows, ], bgs[keep_rows, ])
    absorber     <- keep_rows[nearest_kept]
    geom         <- st_geometry(bgs)
    for(i in seq_along(drop_rows))
      geom[absorber[i]] <- st_union(c(geom[absorber[i]], geom[drop_rows[i]]))
    st_geometry(bgs) <- geom
  }
  bgs <- bgs[keep_rows, ] %>% st_make_valid()
  message(sprintf("Filtered bgs to %d connected populated BGs; absorbed geometry of %d dropped BG(s) into spatial neighbors.",
                  nrow(bgs), length(drop_rows)))

  save(bgs, file = dataset_cache, compress = "xz")
}else{
  load(dataset_cache)
}

### Valid map generation + shortburst 
#There was a lot of complexity involved in getting sims run for different districting
#combinations using different machines and generating cloud storage conflicts. These
#two functions examine how many simulations are completed already and call the appropriate
#functions to run more simulations if the target (250K sims per config) isn't yet met.
#Unconstrained SMC ensemble -- not used in the main paper, but required for the
#unconstrained comparison in Supplementary Materials Section 2 (heterogeneity_bg /
#eg_bg). Cache-gated: no-ops once each config holds 250k plans.
maybe_silent(ensure_smc_all(target_nsims = 250000L, batch_size = 10000L))

#Admin-preserving plans: counties -- DROPPED. The county-only ensemble was
#removed from the analysis: it was degenerate (only a handful of distinct plans
#per config) and Supplementary Materials Section 2 no longer uses it.
#maybe_silent(ensure_smc_all_admin(admin_col = "county_fips", target_nsims = 10000L, batch_size = 10000L))
#Cousub-preserving (municipality-only) — DROPPED. Like the county-only run
#above, this variant is no longer used: nothing downstream builds metrics from
#it (there is no "_cousub_id" entry in .metric_variants), so it was orphaned
#compute. Re-enable if a supplementary section needs the cousub-only ensemble.
#Settings note (if re-enabled): strength = 5 / pop_temper = 0.01 by default;
#smaller configs (d=3..18) sampled cleanly under these settings. d=50 hits a
#bottleneck on the final few splits (the unassigned region runs out of
#cousub-respecting + pop-balanced cuts), so we ease BOTH the constraint
#weight (strength 5 → 3) AND the population-balance enforcement
#(pop_temper 0.01 → 0.03) for that config alone.
#maybe_silent(ensure_smc_all_admin(admin_col = "cousub_id",   target_nsims = 10000L, batch_size = 10000L,
#                                   strength   = c(default = 5,    "50" = 3),
#                                   pop_temper = c(default = 0.01, "50" = 0.03)))

#Combined county+cousub ensemble — SOFT penalties on BOTH county and
#municipality. This is the OLD main, retained for the Supplementary Materials
#(the NEW main is the hard-county ccstruct_m5 ensemble below). Generation logs
#correlated to the on-disk batch PIDs show the production plans were ALL built at
#SOFT strengths county = 1 / cousub = 3 AND pop_temper = 0.01, UNIFORMLY across
#every config including d=50 (d=3 also shows a few abandoned early 3/5 /
#pop_temper=0 attempts, but its production batches are 1/3 at 0.01). The ensemble
#will NOT complete at 1% pop_tol under stronger penalties (see the Supplementary
#Section 3 penalty sweep). Everything is PINNED here because the function's own
#defaults (county/cousub = 3/5, pop_temper = 0) do NOT match production -- leaving
#them unset (as the call originally did) would make a from-scratch rebuild
#diverge from, and fail to reproduce, the analysis ensemble. The fallback ladder
#climbs pop_temper (0.02 -> 0.03) only if 0.01 fails; no strength-raising rung.
maybe_silent(ensure_smc_all_combined(
  target_nsims    = 250000L, batch_size = 10000L,
  county_strength = 1, cousub_strength = 3,
  pop_temper      = 0.01,
  fallback_ladder = list(list(pop_temper = 0.02),
                         list(pop_temper = 0.03))))

#NEW MAIN ensemble — ccstruct_m5 (HARD county via redist_smc(counties=) + SOFT
#municipality penalty at cousub strength 5). PLAN generation folded in from the
#former external Run_CountyStruct.R / Run_CountyStructWorker.R (retained only as
#optional multi-machine accelerators that write the same d{NN}_ccstruct_m5/
#caches). Plans only here, mirroring ensure_smc_all_combined above; the metric
#and Polsby-Popper caches for ccstruct_m5 are produced by the variant loops
#below (generate_metrics_bg / generate_polsby_bg run on these plans). Cache-
#gated per config (ensure_smc_one_county_struct no-ops once a config hits the
#cap). The on-disk ccstruct plans were ALL generated at pop_temper = 0 --
#verified from the worker logs correlated to the batch PIDs, including d=50: the
#HARD county constraint keeps even the fifty-district run feasible at 1% without
#tempering (unlike the SOFT combined run above, which needed 0.02 at d=50). The
#default fallback ladder still climbs pop_temper if any batch fails.
for(spec in district_specs){
  nd <- spec$ndists
  maybe_silent(ensure_smc_one_county_struct(
    ndists       = nd, cousub_strength = 5,
    target_nsims = 250000L, batch_size = 10000L,
    pop_temper   = 0,
    dir_suffix   = "ccstruct_m5"))
}

#Shortburst across both ensembles. The unconstrained pass adds no admin
#constraints; the combined pass stacks county+cousub constraints
#matching the SMC settings (strengths 1, 3). Each writes its own cache
#(shortburst_bg_d{NN}.Rdata vs shortburst_bg_d{NN}_county_cousub.Rdata).
#Unconstrained pass is not surfaced in any current figure (Beyond_the_Border_Wars.qmd
#filters sb_overlay to ensemble == "county_cousub"), so don't pay the
#compute cost. Re-enable if a supplementary section needs the unconstrained
#shortburst envelope.
#Targets searched per (config x ensemble): black_seats_max, black_seats37_max,
#dem_seats_max, rep_seats_max, marg5pt_max, flipprob_max, plus dem_max_k /
#dem_min_k for k = 1..ndists. New competitiveness targets (marg5pt_max,
#flipprob_max) are picked up incrementally on rerun: existing cache entries
#stay, only the two new keys are computed and atomically saved.
#ensure_shortburst_all()
#ccstruct_m5 = NEW main (hard county + soft municipality), overlaid by Beyond_the_Border_Wars.qmd;
#county_cousub = OLD main, retained for the Supplementary Materials. Both caches
#already exist on disk, so these calls are cache-gated no-ops.
maybe_silent(ensure_shortburst_all(admin_col = "ccstruct_m5"))
maybe_silent(ensure_shortburst_all(admin_col = "county_cousub"))

#Motivated multi-seed BlackSeats / BlackSeats37 sweep, applied uniformly to
#every configuration (5000 bursts x 20, 5 diverse seeds, county+cousub
#constraints). Each config gates on its own marker file
#(./Output Data/extra_shortburst_d{NN}_black.Rdata) and merges any improvement
#back into its main county_cousub cache so score_shortburst_optima picks up the
#improved entries automatically.
maybe_silent(ensure_extra_shortburst_all_black(admin_col = "ccstruct_m5", overwrite = overwrite))
maybe_silent(ensure_extra_shortburst_all_black(admin_col = "county_cousub", overwrite = overwrite))

#Generate metrics used in figures from plans.
#DemSeats, ShareDemSeats, WastedVotes, CompetitiveSeats,
#CompetitiveDistricts, BlackSeats, BlackSeats37, BlackWasted,
#BlackPctWasted, plus six PctXxx heterogeneity SDs and Total
# added new metrics based on review:
#n_marginal_2pt, n_marginal_5pt, mean_competitiveness_q,
#prob_seat_change_sigma3, n_county_splits, n_muni_splits
#Metrics scalars + lookup vector. Built once and cached; needed to drive
#every generate_metrics_bg call below.
if(!file.exists("./Output Data/metrics_bg.Rdata") | overwrite==TRUE){
  TotalPop         <- sum(bgs$sum_TotalPop)
  NHBlackPct       <- sum(bgs$NHBlack)/sum(bgs$TotalPop_ACS)
  PctPoverty       <- sum(bgs$InPoverty,na.rm=TRUE)/sum(bgs$PovDeterm,na.rm=TRUE)
  Pct65Plus        <- sum(bgs$POP65Plus,na.rm=TRUE)/TotalPop
  PctHHwChild      <- sum(bgs$HHwithChld)/sum(bgs$Households)
  PctManufacturing <- sum(bgs$EmpManufac)/sum(bgs$EmpTotal)
  Dem20      <- sum(bgs$sum_biden,na.rm=TRUE)
  Rep20      <- sum(bgs$sum_trump,na.rm=TRUE)
  VotesCast  <- Dem20+Rep20
  ShareDem   <- Dem20/VotesCast
  metrics     <- c("PctBlack","PctPoverty","Pct65Plus","PctHHwChild","PctManufacturing","ShareDem")
  metric_vals <- c(NHBlackPct,PctPoverty,Pct65Plus,PctHHwChild,PctManufacturing,ShareDem)
  save(list = c("metrics","metric_vals"), file = "./Output Data/metrics_bg.Rdata")
}else{
  load("./Output Data/metrics_bg.Rdata")
}

#Per-config metric tables across ensemble variants. Each cache is per-file
#gated so a partial run resumes cleanly; per-batch caches under
#./Intermediate Data/bgs/metrics/om/ keep recomputation cheap when only
#some batches are new.
#  "_county_cousub"  — combined county+cousub ensemble (MAIN analytical run)
#The county-only ("_county_fips") and unconstrained (admin_col = NULL,
#suffix = "") variants are commented out -- the main paper uses _county_cousub
#and the supplementary materials no longer build a county-only overlay
#(it was degenerate). Re-enable an entry below to reinstate it.
#All built objects stay in scope -- the heterogeneity/eg/t_comp variant
#builders below look them up by name via .om_for().
.metric_variants <- list(
  list(admin_col = NULL,            suffix = ""),               # unconstrained BG (Supplementary Section 2)
  #list(admin_col = "county_fips",   suffix = "_county_fips"),
  list(admin_col = "ccstruct_m5",   suffix = "_ccstruct_m5"),  # NEW main (Beyond_the_Border_Wars.qmd)
  list(admin_col = "county_cousub", suffix = "_county_cousub") # OLD main (Supplementary)
)
for(v in .metric_variants){
  for(spec in district_specs){
    nd     <- spec$ndists
    out_f  <- sprintf("./Output Data/om%02d_bg%s.Rdata", nd, v$suffix)
    var_nm <- sprintf("output_metrics_%d%s", nd, v$suffix)
    if(file.exists(out_f) & !overwrite){
      load(out_f)
      next
    }
    om <- maybe_silent(generate_metrics_bg(metrics, metric_vals,
                                            ndists    = nd, bgs,
                                            admin_col = v$admin_col,
                                            overwrite = overwrite))
    assign(var_nm, om)
    save(list = var_nm, file = out_f)
  }
}
rm(.metric_variants)

#Per-variant helpers for the heterogeneity, eg, and t_comp tables.
#  variant_suffix = "" -> unconstrained ensemble (object names like
#                        output_metrics_17), saved without a suffix
#                        (heterogeneity_bg.Rdata etc.)
#  variant_suffix = "_county_cousub" -> constrained ensemble, saved as
#                        heterogeneity_bg_county_cousub.Rdata etc.
.om_for <- function(nd, variant_suffix){
  nm <- sprintf("output_metrics_%d%s", nd, variant_suffix)
  if(!exists(nm, inherits = TRUE)) return(NULL)
  get(nm, inherits = TRUE)
}

.build_heterogeneity_bg <- function(variant_suffix, dst_file){
  if(file.exists(dst_file)) return(invisible(NULL))
  type_labels <- c("3"  = "Three Districts",   "6"  = "Six Districts",
                   "9"  = "Nine Districts",    "17" = "Seventeen Districts",
                   "18" = "Eighteen Districts","50" = "Fifty Districts")
  #Single-member comparator is d=18 (evenly divides multi-member plans
  #into 2-, 3-, and 6-member districts)
  pieces <- list()
  for(nd in c(3, 6, 9, 18, 50)){
    om <- .om_for(nd, variant_suffix)
    if(is.null(om)) next
    d <- reshape2::melt(data = om[[1]], id.vars = "Plan",
                        measure.vars = c("Pct65Plus","PctHHwChild","PctManufacturing",
                                          "PctPoverty","PctBlack","ShareDem","Total"),
                        variable.name = "Metric",
                        value.name    = "DistrictSquaredError")
    d$Type <- type_labels[as.character(nd)]
    pieces[[length(pieces)+1L]] <- d
  }
  if(length(pieces) == 0L) return(invisible(NULL))
  heterogeneity <- do.call(rbind, pieces)
  colnames(heterogeneity)[colnames(heterogeneity) == "DistrictSquaredError"] <- "Value"
  heterogeneity <- heterogeneity %>%
    mutate(Metric = fct_relevel(Metric, "Pct65Plus","PctHHwChild","PctManufacturing",
                                "PctPoverty","PctBlack","ShareDem"))
  heterogeneity$Type <- factor(heterogeneity$Type,
                                levels = c("Three Districts","Six Districts","Nine Districts",
                                           "Eighteen Districts","Fifty Districts"))
  levels(heterogeneity$Type) <- addline_format(levels(heterogeneity$Type))
  heterogeneity <- heterogeneity[heterogeneity$Metric != "Total", ]
  save(heterogeneity, file = dst_file)
  invisible(heterogeneity)
}

.build_eg_bg <- function(variant_suffix, dst_file){
  if(file.exists(dst_file)) return(invisible(NULL))
  base_cols <- c("ShareDemSeats","CompetitiveSeats","CompetitiveDistricts",
                 "WastedVotes","BlackSeats","BlackWasted","BlackPctWasted","BlackSeats37")
  new_cols  <- c("n_marginal_2pt","n_marginal_5pt","mean_competitiveness_q",
                 "prob_seat_change_sigma3","n_county_splits","n_muni_splits")
  pieces <- list()
  for(nd in c(3, 6, 9, 17, 18, 50)){
    om <- .om_for(nd, variant_suffix)
    if(is.null(om)) next
    data <- om[[1]]
    keep <- intersect(c(base_cols, new_cols), colnames(data))
    sub  <- data[, keep, drop = FALSE]
    sub$Type <- as.character(nd)
    pieces[[length(pieces)+1L]] <- sub
  }
  if(length(pieces) == 0L) return(invisible(NULL))
  eg <- do.call(rbind, pieces)
  eg$Type <- factor(eg$Type, levels = c(3,6,9,17,18,50),
                    labels = c("Three Districts","Six Districts","Nine Districts",
                               "Seventeen Districts","Eighteen Districts","Fifty Districts"))
  save(eg, file = dst_file)
  invisible(eg)
}

.build_t_comp_bg <- function(variant_suffix, dst_file){
  if(file.exists(dst_file)) return(invisible(NULL))
  t_comp  <- bgs
  joined  <- 0L
  for(nd in c(3, 6, 9, 17, 18, 50)){
    om <- .om_for(nd, variant_suffix)
    if(is.null(om)) next
    data <- om[[2]]
    colnames(data)[-1] <- paste0(colnames(data)[-1], "_", nd)
    t_comp <- left_join(t_comp, data, by = "GEOID")
    joined <- joined + 1L
  }
  if(joined == 0L) return(invisible(NULL))   # nothing to write
  save(t_comp, file = dst_file)
  invisible(t_comp)
}

#Build (or load) the ensemble variants. Only the MAIN county+cousub variant
#is built now.
.variants <- list(
  list(suffix = "",               dst_tag = "",               label = "unconstrained BG (Supplementary Section 2)"),
  list(suffix = "_ccstruct_m5",   dst_tag = "_ccstruct_m5",   label = "ccstruct_m5 hard-county+soft-muni (NEW MAIN, GA)"),
  list(suffix = "_county_cousub", dst_tag = "_county_cousub", label = "county+cousub (OLD main, Supplementary)")
  #county-only variant dropped (degenerate)
  #list(suffix = "_county_fips",   dst_tag = "_county_fips",   label = "county-preserving (supplementary)")
)

for(v in .variants){
  het_f <- sprintf("./Output Data/heterogeneity_bg%s.Rdata", v$dst_tag)
  eg_f  <- sprintf("./Output Data/eg_bg%s.Rdata",            v$dst_tag)
  tc_f  <- sprintf("./Output Data/t_comp_bg%s.Rdata",        v$dst_tag)
  .build_heterogeneity_bg(v$suffix, het_f)
  .build_eg_bg           (v$suffix, eg_f)
  .build_t_comp_bg       (v$suffix, tc_f)
}

#Polsby-Popper compactness across the ensemble. Standalone from
#generate_metrics_bg so the geometry-dependent computation (which is
#unchanged across the seat/wasted-vote metric edits) caches
#independently. Saves a per-config plan-level table under
#./Output Data/pp{NN}_bg{suffix}.Rdata, then concatenates them into a
#single `polsby` data frame parallel to `eg`, saved at
#./Output Data/polsby_bg{suffix}.Rdata. The aggregate has one row per
#plan with columns Type, mean_polsby_popper, min_polsby_popper,
#max_polsby_popper.
.pp_variants <- list(
  list(admin_col = "ccstruct_m5",   suffix = "_ccstruct_m5"),  # NEW main (Beyond_the_Border_Wars.qmd)
  list(admin_col = "county_cousub", suffix = "_county_cousub") # OLD main (Supplementary)
  #county-only variant dropped (see .metric_variants above)
  #list(admin_col = "county_fips",   suffix = "_county_fips")
)
.pp_perim <- NULL   #precomputed once on first use and reused across configs
for(v in .pp_variants){
  agg_f <- sprintf("./Output Data/polsby_bg%s.Rdata", v$suffix)
  if(file.exists(agg_f) & !overwrite) next
  type_labels <- c("3"  = "Three Districts",   "6"  = "Six Districts",
                   "9"  = "Nine Districts",    "17" = "Seventeen Districts",
                   "18" = "Eighteen Districts","50" = "Fifty Districts")
  pp_pieces <- list()
  for(spec in district_specs){
    nd    <- spec$ndists
    out_f <- sprintf("./Output Data/pp%02d_bg%s.Rdata", nd, v$suffix)
    if(file.exists(out_f) & !overwrite){
      load(out_f)   #populates pp_NN
    }else{
      if(is.null(.pp_perim)){
        message("polsby: precomputing perim_df ...")
        .pp_perim <- redistmetrics::prep_perims(shp = bgs)
      }
      pp_df <- maybe_silent(generate_polsby_bg(ndists    = nd,
                                                bgs.full  = bgs,
                                                admin_col = v$admin_col,
                                                perim_df  = .pp_perim,
                                                overwrite = overwrite))
      assign(sprintf("pp_%d", nd), pp_df)
      save(list = sprintf("pp_%d", nd), file = out_f)
    }
    pp_df <- get(sprintf("pp_%d", nd))
    if(nrow(pp_df) == 0L){ rm(list = sprintf("pp_%d", nd)); next }
    pp_df$Type <- type_labels[as.character(nd)]
    pp_pieces[[length(pp_pieces)+1L]] <- pp_df
    rm(list = sprintf("pp_%d", nd))
  }
  if(length(pp_pieces) > 0L){
    polsby <- do.call(rbind, pp_pieces)
    polsby$Type <- factor(polsby$Type,
                          levels = c("Three Districts","Six Districts","Nine Districts",
                                     "Seventeen Districts","Eighteen Districts","Fifty Districts"))
    polsby <- polsby[, c("Type","mean_polsby_popper","min_polsby_popper","max_polsby_popper")]
    save(polsby, file = agg_f)
    rm(polsby)
  }
}
rm(.pp_variants, .pp_perim)

#Per-district Dem-share distribution + seat-outcome classes for the d=9 anomaly
#discussion (Supplementary Materials Section 4: fig-supp-rank-shares). Like the
#Polsby-Popper block above this reads the stored SMC plan matrices directly, so
#it is built once per ensemble variant: the hard-county MAIN (ccstruct_m5) that
#the figure now uses, plus the soft county+cousub ensemble retained for the
#Supplementary contrast. Formerly the standalone verify_anomaly_numbers.R;
#folded in here so the consumed cache is built with everything else and honors
#`overwrite`. Output: anomaly_diagnostics{suffix}.Rdata holding rank_summary /
#outcome_summary / per_plan_outcomes (d = 3,6,9,17,18; d=50 excluded by design).
.anom_variants <- list(
  list(admin_col = "ccstruct_m5",   suffix = "_ccstruct_m5"),  # NEW main (GA + SI Section 4)
  list(admin_col = "county_cousub", suffix = "_county_cousub") # OLD main (SI contrast)
)
for(v in .anom_variants){
  anom_f <- sprintf("./Output Data/anomaly_diagnostics%s.Rdata", v$suffix)
  if(file.exists(anom_f) & !overwrite) next
  .anom <- maybe_silent(build_anomaly_bg(bgs, admin_col = v$admin_col))
  rank_summary      <- .anom$rank_summary
  outcome_summary   <- .anom$outcome_summary
  per_plan_outcomes <- .anom$per_plan_outcomes
  save(rank_summary, outcome_summary, per_plan_outcomes, file = anom_f)
  rm(.anom, rank_summary, outcome_summary, per_plan_outcomes)
}
rm(.anom_variants)

#Tract-resolution comparison ensemble (Supplementary Materials Section 2).
#The tract-level SMC ensemble is built at pop_tol = 0.01 with NO
#administrative-preservation constraints, matching the block-group
#UNCONSTRAINED ensemble on tolerance and constraints so the contrast isolates
#sub-unit resolution (census tract vs block group). Generation now lives in
#this pipeline via ensure_tract_smc_all() (the helpers were lifted out of the
#former standalone Run_Tract1Pct.R); it is cache-gated per config, so missing
#tract plans are regenerated and existing ones no-op.
#Here we then recompute the FULL headline metric set on those tract plans —
#the same heterogeneity SDs, seat share, wasted votes, Black-preferred
#seats, competitiveness, and admin splits that generate_metrics_bg()
#produces for block groups — and assemble tidy heterogeneity / eg / polsby
#objects parallel to the block-group builders above. Gated on output-file
#existence, so this heavy scoring runs once and then no-ops.
.tract_out <- c("./Output Data/heterogeneity_tract1pct.Rdata",
                "./Output Data/eg_tract1pct.Rdata",
                "./Output Data/polsby_tract1pct.Rdata")
if(!all(file.exists(.tract_out)) | overwrite){
  message("tract1pct: recomputing full metric set on the tract-resolution ensemble ...")
  #Ensure the tract SMC plans exist (generate any missing config; cache-gated).
  maybe_silent(ensure_tract_smc_all(nsims = 50000L))
  .ti  <- local({ load("./Intermediate Data/tracts.interpolated.Rdata"); tracts.interpolated })
  .trc <- local({ load("./Intermediate Data/tract_1pct/new_metrics/tract_admin_join.Rdata"); trc })
  .u <- sf::st_drop_geometry(.ti); .u <- .u[.u$sum_TotalPop > 0, ]
  stopifnot(identical(as.character(.u$GEOID), as.character(.trc$GEOID)))
  .u$county_fips <- .trc$county_fips; .u$cousub_id <- .trc$cousub_id
  .u$Voters <- .u$sum_biden + .u$sum_trump

  .tmetrics     <- c("PctBlack","PctPoverty","Pct65Plus","PctHHwChild","PctManufacturing","ShareDem")
  .tmetric_vals <- c(
    sum(.u$NHBlack)/sum(.u$TotalPop_ACS),
    sum(.u$InPoverty, na.rm=TRUE)/sum(.u$PovDeterm, na.rm=TRUE),
    sum(.u$POP65Plus, na.rm=TRUE)/sum(.u$sum_TotalPop),
    sum(.u$HHwithChld)/sum(.u$Households),
    sum(.u$EmpManufac)/sum(.u$EmpTotal),
    sum(.u$sum_biden)/(sum(.u$sum_biden)+sum(.u$sum_trump)))
  names(.tmetric_vals) <- .tmetrics

  .tcols <- c("sum_biden","Voters","sum_BVAP","sum_VAP","NHBlack","TotalPop_ACS",
              "InPoverty","PovDeterm","POP65Plus","sum_TotalPop","HHwithChld",
              "Households","EmpManufac","EmpTotal")
  .X <- as.matrix(.u[, .tcols]); .X[is.na(.X)] <- 0
  .cidx <- setNames(seq_along(.tcols), .tcols)
  .cnty <- .u$county_fips; .okc <- !is.na(.cnty)
  .muni <- .u$cousub_id;   .okm <- !is.na(.muni)

  #Score one block of plan columns; mirrors generate_metrics_bg()'s math.
  .score_tract_block <- function(M, nd, nm){
    b <- ncol(M); q <- 1/(nm+1)
    PD <- array(0, dim = c(nd, length(.tcols), b))
    for(p in seq_len(b)){ rs <- rowsum(.X, M[,p], reorder=TRUE); PD[as.integer(rownames(rs)),,p] <- rs }
    asm <- function(z) if(is.matrix(z)) z else matrix(z, nrow=nd)
    g <- function(n) asm(PD[, .cidx[n], ])
    sb<-g("sum_biden"); vt<-g("Voters"); pb<-g("sum_BVAP"); pvap<-g("sum_VAP")
    nhb<-g("NHBlack"); tpacs<-g("TotalPop_ACS"); pov<-g("InPoverty"); povd<-g("PovDeterm")
    p65<-g("POP65Plus"); tp<-g("sum_TotalPop"); hhc<-g("HHwithChld"); hh<-g("Households")
    empm<-g("EmpManufac"); empt<-g("EmpTotal")
    dp<-sb/vt
    toWin<-vt/(nm+1); dem_w<-floor(sb/toWin); rep_w<-floor((vt-sb)/toWin)
    DemSeats<-colSums(dem_w); RepSeats<-colSums(rep_w)
    ShareDemSeats<-DemSeats/pmax(DemSeats+RepSeats,1)
    WastedVotes<-(colSums(sb-dem_w*toWin)-colSums((vt-sb)-rep_w*toWin))/colSums(vt)
    CompetitiveDistricts<-colSums((dp>0.45)&(dp<0.55))
    CompetitiveSeats<-if(nm==1) CompetitiveDistricts else colSums((dem_w!=0)&(rep_w!=0))
    bkw<-pvap/(nm+1); bkw2<-pvap/(nm+2)
    bw<-floor(pb/bkw); bw37<-floor(pb/bkw2); ow<-floor((pvap-pb)/bkw)
    BlackSeats<-colSums(bw); BlackSeats37<-colSums(bw37)
    BlackWasted<-(colSums(pb-bw*bkw)-colSums((pvap-pb)-ow*bkw))/colSums(pvap)
    BlackPctWasted<-(colSums(pb-bw*bkw)/colSums(pb))-(colSums((pvap-pb)-ow*bkw)/colSums(pvap-pb))
    vmod<-dp %% q; delta<-pmin(vmod,q-vmod)
    out <- data.frame(
      PctBlack=sd_metric(nhb,tpacs,.tmetric_vals["PctBlack"]),
      PctPoverty=sd_metric(pov,povd,.tmetric_vals["PctPoverty"]),
      Pct65Plus=sd_metric(p65,tp,.tmetric_vals["Pct65Plus"]),
      PctHHwChild=sd_metric(hhc,hh,.tmetric_vals["PctHHwChild"]),
      PctManufacturing=sd_metric(empm,empt,.tmetric_vals["PctManufacturing"]),
      ShareDem=sd_metric(sb,vt,.tmetric_vals["ShareDem"]),
      DemSeats=DemSeats, ShareDemSeats=ShareDemSeats,
      CompetitiveSeats=CompetitiveSeats, CompetitiveDistricts=CompetitiveDistricts,
      WastedVotes=WastedVotes, BlackSeats=BlackSeats, BlackSeats37=BlackSeats37,
      BlackWasted=BlackWasted, BlackPctWasted=BlackPctWasted,
      n_marginal_2pt=colSums(delta<0.02), n_marginal_5pt=colSums(delta<0.05),
      mean_competitiveness_q=colMeans(1-2*(nm+1)*delta),
      prob_seat_change_sigma3=colMeans(flip_prob_q(delta,q,0.03)),
      row.names=NULL)
    out$n_county_splits <- vapply(seq_len(b), function(p)
      if(any(.okc)) sum(rowSums(table(.cnty[.okc],M[.okc,p])>0)>1) else NA_integer_, integer(1))
    out$n_muni_splits <- vapply(seq_len(b), function(p)
      if(any(.okm)) sum(rowSums(table(.muni[.okm],M[.okm,p])>0)>1) else NA_integer_, integer(1))
    out
  }

  .tspecs <- c(3,6,9,17,18,50)
  .tplan <- list(); .tpp <- list()
  for(nd in .tspecs){
    nm <- 18/nd; if(nm < 2) nm <- 1
    .pe <- new.env()
    load(sprintf("./Intermediate Data/tract_1pct/smc_plans/pa_plans_tract1pct_d%02d.Rdata", nd), envir=.pe)
    Mall <- attr(.pe[[ls(.pe)[1]]], "plans")[, -1, drop=FALSE]
    nplans <- ncol(Mall); pieces <- list(); done <- 0L
    while(done < nplans){
      idx <- (done+1L):min(done+10000L, nplans)
      pieces[[length(pieces)+1L]] <- .score_tract_block(Mall[, idx, drop=FALSE], nd, nm)
      done <- done + length(idx)
    }
    tab <- do.call(rbind, pieces); tab$Plan <- seq_len(nrow(tab)); tab$ndists <- nd
    .tplan[[as.character(nd)]] <- tab
    .ppf <- sprintf("./Intermediate Data/tract_1pct/new_metrics/pp_d%02d.Rdata", nd)
    if(file.exists(.ppf)){
      .ppe <- new.env(); load(.ppf, envir=.ppe); ppv <- .ppe[[ls(.ppe)[1]]]
      ppm <- matrix(ppv, nrow=nd)[, -1, drop=FALSE]
      .tpp[[as.character(nd)]] <- data.frame(ndists=nd,
        mean_polsby_popper=colMeans(ppm, na.rm=TRUE),
        min_polsby_popper=apply(ppm,2,min,na.rm=TRUE),
        max_polsby_popper=apply(ppm,2,max,na.rm=TRUE))
    }
    rm(.pe, Mall); gc(verbose=FALSE)
  }

  .tlab <- c("3"="Three Districts","6"="Six Districts","9"="Nine Districts",
             "17"="Seventeen Districts","18"="Eighteen Districts","50"="Fifty Districts")
  .tlev <- c("Three Districts","Six Districts","Nine Districts",
             "Seventeen Districts","Eighteen Districts","Fifty Districts")

  het_pieces <- list()
  for(nd in c(3,6,9,18,50)){
    d <- .tplan[[as.character(nd)]]
    for(m in c("Pct65Plus","PctHHwChild","PctManufacturing","PctPoverty","PctBlack","ShareDem"))
      het_pieces[[length(het_pieces)+1L]] <- data.frame(Plan=d$Plan, Metric=m,
        Value=d[[m]], Type=.tlab[as.character(nd)])
  }
  heterogeneity <- do.call(rbind, het_pieces)
  heterogeneity$Metric <- factor(heterogeneity$Metric,
    levels=c("Pct65Plus","PctHHwChild","PctManufacturing","PctPoverty","PctBlack","ShareDem"))
  heterogeneity$Type <- factor(heterogeneity$Type, levels=.tlev)
  save(heterogeneity, file="./Output Data/heterogeneity_tract1pct.Rdata")

  eg_cols <- c("ShareDemSeats","CompetitiveSeats","CompetitiveDistricts","WastedVotes",
               "BlackSeats","BlackWasted","BlackPctWasted","BlackSeats37","n_marginal_2pt",
               "n_marginal_5pt","mean_competitiveness_q","prob_seat_change_sigma3",
               "n_county_splits","n_muni_splits")
  eg <- do.call(rbind, lapply(.tspecs, function(nd){
    s <- .tplan[[as.character(nd)]][, eg_cols]; s$Type <- .tlab[as.character(nd)]; s }))
  eg$Type <- factor(eg$Type, levels=.tlev)
  save(eg, file="./Output Data/eg_tract1pct.Rdata")

  if(length(.tpp)){
    polsby <- do.call(rbind, lapply(.tspecs, function(nd){
      p <- .tpp[[as.character(nd)]]; if(is.null(p)) return(NULL)
      p$Type <- .tlab[as.character(nd)]; p }))
    polsby$Type <- factor(polsby$Type, levels=.tlev)
    polsby <- polsby[, c("Type","mean_polsby_popper","min_polsby_popper","max_polsby_popper")]
    save(polsby, file="./Output Data/polsby_tract1pct.Rdata")
  }
  rm(.ti,.trc,.u,.X,.tplan,.tpp,.score_tract_block,het_pieces); gc(verbose=FALSE)
}
rm(.tract_out)

#Bring the NEW main ccstruct_m5 variant (hard county + soft municipality)
#into scope under the unsuffixed names that Beyond_the_Border_Wars.qmd expects. The
#OLD county_cousub variant stays on disk (*_county_cousub.Rdata) for the
#Supplementary Materials, which load it explicitly where needed.
load("./Output Data/heterogeneity_bg_ccstruct_m5.Rdata")  # heterogeneity
load("./Output Data/eg_bg_ccstruct_m5.Rdata")             # eg
load("./Output Data/t_comp_bg_ccstruct_m5.Rdata")         # t_comp
load("./Output Data/polsby_bg_ccstruct_m5.Rdata")         # polsby

rm(.om_for, .build_heterogeneity_bg, .build_eg_bg, .build_t_comp_bg, .variants)

#Per-BG Droop-marginal fractions (one Marginal_{NN} column per config).
#Computed from the NEW main ccstruct_m5 ensemble by build_marginal_bg
marginal_bg_f <- file.path(bg_root, "marginal_bg.Rdata")
if(!file.exists(marginal_bg_f) | overwrite == TRUE){
  marginal_bg <- build_marginal_bg(bgs, admin_col = "ccstruct_m5")
  tmp_f <- paste0(marginal_bg_f, ".tmp")
  save(marginal_bg, file = tmp_f, compress = "xz")
  file.rename(tmp_f, marginal_bg_f)
}else{
  load(marginal_bg_f)   # populates marginal_bg
}
if(inherits(marginal_bg, "sf")) marginal_bg <- sf::st_drop_geometry(marginal_bg)
t_comp <- left_join(t_comp, marginal_bg, by = "GEOID")


#Pull out maximizing plans from shortburst runs
if(!file.exists("./Output Data/sb_overlay.Rdata") | overwrite==TRUE){
  sb_overlay <- score_shortburst_optima(bgs)
  save(sb_overlay, file = "./Output Data/sb_overlay.Rdata")
}else{
  load("./Output Data/sb_overlay.Rdata")
}

#Congressional Black Caucus membership vs district BVAP for U.S. House districts
#(Supplementary Materials BVAP-threshold figure).
if(!file.exists("./Output Data/BlackCBC.Rdata") | overwrite==TRUE){
  BlackCBC <- read.csv("./Input Data/PresidentialResultsAndBVAP.csv",
                       stringsAsFactors = FALSE)
  BlackCBC$bvap_pct  <- BlackCBC$bvap_pct * 100
  BlackCBC$is_in_cbc <- factor(BlackCBC$is_in_cbc, levels = c(TRUE, FALSE),
                               labels = c("Yes", "No"))
  save(BlackCBC, file = "./Output Data/BlackCBC.Rdata")
}else{
  load("./Output Data/BlackCBC.Rdata")
}

#Block-group-level diagnostics for Pennsylvania's two most recent enacted
#congressional plans (Supplementary Materials Section 3). Precompute here so
#the qmd just loads the result rather than rerunning prep_perims/comp_polsby.
if(!file.exists("./Output Data/enacted_metrics.Rdata") | overwrite==TRUE){
  .enacted_perim <- redistmetrics::prep_perims(shp = bgs)
  m116 <- enacted_metrics("RoughCD116", bgs, .enacted_perim)  # 116th Congress: 18 districts
  m118 <- enacted_metrics("RoughCD118", bgs, .enacted_perim)  # current 118th Congress: 17 districts
  save(m116, m118, file = "./Output Data/enacted_metrics.Rdata")
  rm(.enacted_perim)
}else{
  load("./Output Data/enacted_metrics.Rdata")
}

#Soft-penalty sweep GENERATION (folded in from the former standalone
#Run_PenaltySweep.R, which is retained only as an optional multi-machine
#accelerator writing the same caches). For each stronger SOFT-penalty cell, run
#the MAIN combined pipeline (SMC -> metrics -> Polsby-Popper) into an ISOLATED
#variant dir. Cache-gated per cell so it no-ops once populated (and skips the
#prep_perims cost entirely on the cache-hit path). 2 000-plan feasibility runs:
#a single SMC batch is near-deterministic on split counts, so the medians are
#stable. The 1/3 "main" reference cell is the county_cousub ensemble generated
#above (line ~305), so it is NOT regenerated here.
.penalty_gen_cells <- list(
  list(variant = "cc_c5_m8_pilot",   county = 5,  cousub = 8),
  list(variant = "cc_c8_m12_pilot",  county = 8,  cousub = 12),
  list(variant = "cc_c12_m20_pilot", county = 12, cousub = 20)
)
.cell_populated <- function(variant) all(vapply(c(17L, 18L), function(nd)
  length(.list_variant_files(file.path("metrics", "om"),     variant, nd)) > 0L &&
  length(.list_variant_files(file.path("metrics", "polsby"), variant, nd)) > 0L,
  logical(1)))
.penalty_perim <- NULL
for(cell in .penalty_gen_cells){
  if(.cell_populated(cell$variant) & !overwrite) next
  if(is.null(.penalty_perim)) .penalty_perim <- redistmetrics::prep_perims(shp = bgs)
  maybe_silent(ensure_penalty_sweep_cell(
    county_strength = cell$county, cousub_strength = cell$cousub,
    variant = cell$variant, configs = c(17, 18),
    target_nsims = 2000L, batch_size = 2000L,
    metrics = metrics, metric_vals = metric_vals, bgs = bgs,
    perim_df = .penalty_perim))
}
rm(.penalty_gen_cells, .cell_populated, .penalty_perim)

#Soft-penalty sweep table (Supplementary Materials Section 3: "Can stronger
#preservation penalties close the gap?"). Ensemble q05/median/q95 for county
#splits, municipality splits, and Polsby-Popper across increasing SOFT
#county/cousub penalties, against the enacted values. Assembled from the
#isolated penalty-cell variant dirs generated by Run_PenaltySweep.R (which
#stays the generator); needs enacted_metrics (above) for the reference values
#via .enacted_for(). Formerly assembled by hand in an interactive session --
#folded in here so the consumed cache is built by the pipeline and honors
#`overwrite`. The cells are the SOFT-penalty family on purpose: this table is
#the evidence that turning up the soft penalty does NOT close the gap, which is
#what motivated the hard-county ccstruct_m5 main ensemble.
if(!file.exists("./Output Data/penalty_sweep.Rdata") | overwrite==TRUE){
  .penalty_cells <- list(
    list(variant = "county_cousub",    label = "1 / 3 (main)", county = 1,  cousub = 3),
    list(variant = "cc_c5_m8_pilot",   label = "5 / 8",        county = 5,  cousub = 8),
    list(variant = "cc_c8_m12_pilot",  label = "8 / 12",       county = 8,  cousub = 12),
    list(variant = "cc_c12_m20_pilot", label = "12 / 20",      county = 12, cousub = 20)
  )
  penalty_sweep <- assemble_penalty_sweep(.penalty_cells, configs = c(17, 18))
  save(penalty_sweep, file = "./Output Data/penalty_sweep.Rdata")
  rm(.penalty_cells, penalty_sweep)
}

#Clean up the workspace for the manuscript. Beyond_the_Border_Wars.qmd consumes only
#bgs, heterogeneity, t_comp, eg, and sb_overlay -- everything else is
#scaffolding from the build pipeline.
rm(
   #control variables and path/constant strings used only during this script
   overwrite, maybe_silent,
   pa_crs, bg_root, bg_smc_dir, bg_sb_dir, dataset_cache, marginal_bg_f,
   #functions sourced from BBW - Custom Functions.R, used only during this script
   tcw_blocks, .plans_matrix, list_smc_batches,
   smc_nsims_accumulated, get_one_smc_init_plan, build_pa_map_bg,
   ensure_smc_one, ensure_smc_all,
   ensure_smc_one_admin, ensure_smc_all_admin,
   ensure_smc_one_combined, ensure_smc_all_combined,
   build_tract_shp, tract_district_specs,
   ensure_tract_smc_one, ensure_tract_smc_all,
   ensure_shortburst_one, ensure_shortburst_all,
   ensure_extra_shortburst_black, ensure_extra_shortburst_all_black,
   .shortburst_plan_matrix, .score_plans_topk,
   make_seats_scorer, make_black_seats_scorer,
   make_competitiveness_scorer, .score_plans_competitiveness, flip_prob_q,
   .score_plans_seats, .score_plans_black_seats, find_best_seeds,
   sd_metric, build_marginal_bg,
   generate_metrics_bg, generate_polsby_bg, enacted_metrics,
   score_shortburst_optima, addline_format,
   #lookup tables only needed to drive generate_metrics_bg
   metrics, metric_vals,
   #data objects rolled up into heterogeneity/eg/t_comp or unused by qmd.
   #BlackCBC is kept in scope -- supplementary materials will use it.
   marginal_bg, district_specs)

#Scalar denominators only enter scope when the metrics_bg.Rdata cache is
#built fresh; on the cache-hit path they are never assigned. Use
#intersect() so the rm doesn't warn about names that don't exist.
.scalar_names <- c("TotalPop","NHBlackPct","PctPoverty","Pct65Plus",
                   "PctHHwChild","PctManufacturing",
                   "Dem20","Rep20","VotesCast","ShareDem")
rm(list = intersect(.scalar_names, ls()))
rm(.scalar_names)

#Per-config metric tables (output_metrics_*) have been folded into the
#heterogeneity/eg/t_comp summaries. Drop all variants -- whichever combination
#ended up in scope (now just _county_cousub).
rm(list = ls(pattern = "^output_metrics_"))

#Loop iteration variables left behind by the metric-build, dir-creation,
#and variant-build loops. Use intersect() because not all of them exist on
#every code path (e.g. eg_f / het_f / tc_f are only assigned when the
#variant loop ran a fresh build).
.loop_leftovers <- c("v", "spec", "d", "nd", "out_f", "var_nm", "om",
                     "het_f", "eg_f", "tc_f",
                     "agg_f", "type_labels", "pp_pieces", "pp_df")
rm(list = intersect(.loop_leftovers, ls()))
rm(.loop_leftovers)
