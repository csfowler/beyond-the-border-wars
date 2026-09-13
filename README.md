# Beyond the Border Wars: Representation, Scale, and Multi-Member Districts

Replication code for:

> Fowler, Christopher S., and Linda L. Fowler. 2026. "Beyond the Border Wars:
> Representation, Scale, and Multi-Member Districts." *Geographical Analysis*.
> https://doi.org/10.1111/gean.70057

This repository contains the version of the code and manuscript source accepted for
publication. Volume, issue, and page numbers will be added once the article is assigned
to an issue.

## Citation

If you use this code, please cite the article:

```bibtex
@article{fowler2026border,
  title   = {Beyond the Border Wars: Representation, Scale, and Multi-Member Districts},
  author  = {Fowler, Christopher S. and Fowler, Linda L.},
  journal = {Geographical Analysis},
  year    = {2026},
  doi     = {10.1111/gean.70057}
}
```

## Repository contents

| File | Purpose |
| --- | --- |
| `BBW - Data Prep BG.R` | Full pipeline: builds the block-group dataset, generates or resumes the SMC ensembles and short-burst optima, and computes every metric object used by the manuscript and supplement. |
| `BBW - Custom Functions.R` | Helper functions sourced by the pipeline (vote apportionment, SMC/short-burst drivers, metric calculation, caching). |
| `Beyond_the_Border_Wars.qmd` | The accepted text. Renders the full paper. |
| `Supplementary_Materials.qmd` | Supplementary Materials source. |
| `Input Data/PresidentialResultsAndBVAP.csv` | U.S. House district presidential results, BVAP, and Congressional Black Caucus membership; used for the BVAP-threshold comparison in the Supplementary Materials. |
| `references.bib`, `custom-reference-doc.docx`, `elsevier-harvard.csl`, `_extensions/` | Bibliography, Word template, citation style, and Quarto/Elsevier formatting support. |
## Data

Most inputs are downloaded automatically the first time the pipeline runs:

- **2020 Decennial Census** (PL 94-171) block population, voting-age population, and Black
  voting-age population, and **2016–2020 ACS 5-year** block-group and tract estimates,
  via [`tidycensus`](https://walker-data.com/tidycensus/). A Census API key is required
  (`tidycensus::census_api_key()`).
- **TIGER/Line 2020** county subdivisions, plus 116th/118th congressional district and
  state senate boundaries, downloaded directly from the Census Bureau.

Two inputs are **not** included in this repository. Download them and place them in `Input Data/`:

- `pa_vtds20.*`: 2020 presidential vote totals by voting district for Pennsylvania,
  compiled by Ruth Buck from Census voting districts and OpenElections returns.
  The shapefile and its documentation are available at
  [github.com/RKBuck1/pa_vtds_2020](https://github.com/RKBuck1/pa_vtds_2020).
  Section 1 of the Supplementary Materials compares this file with the VEST precinct data.
- `pa_2020/pa_2020.*`: the Voting and Election Science Team (VEST) 2020 Pennsylvania
  precinct file ([Harvard Dataverse, doi:10.7910/DVN/K7760H](https://doi.org/10.7910/DVN/K7760H)),
  used only for the data-validation comparison in the Supplementary Materials.

## Software
Analyses were run with R 4.6.0 and Quarto 1.7. Key package versions:
| Package | Version |
| --- | --- |
| redist | 4.3.2 |
| redistmetrics | 1.0.11 |
| sf | 1.1.0 |
| tidyverse | 2.0.0 |
| tidycensus | 1.7.5 |
| tigris | 2.2.1 |
| igraph | 2.3.0 |
| ggpubr | 0.6.3 |
| gridExtra | 2.3 |
| patchwork | 1.3.2 |
| RColorBrewer | 1.1.3 |
| reshape2, knitr, scales | current CRAN |

## Reproducing the results

Open `Beyond the Border Wars.Rproj`, or set the working directory to the repository root.

1. **Run the pipeline.**

   ```r
   source("BBW - Data Prep BG.R")
   ```
2. **Render the manuscript and supplement.**

   ```sh
   quarto render Beyond_the_Border_Wars.qmd --to docx
   quarto render Supplementary_Materials.qmd
   ```

   `Beyond_the_Border_Wars.qmd` sources `BBW - Data Prep BG.R` itself, so rendering it on its
   own runs the entire pipeline and produces the complete paper with Table 1 and all figures.
   Running step 1 first simply separates the long computation from the render.

Simulation uses random seeds, so a from-scratch rebuild will produce a statistically
equivalent ensemble rather than identical plans. Distributional summaries and the
figures' substantive patterns should match the published results.