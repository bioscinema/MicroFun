# MicroFun <img src="man/figures/logo.png" align="right" width="160" alt="MicroFun logo" />

MicroFun is an R package for differential analysis and visualization of taxon-stratified microbial functional profiles predicted by PICRUSt2. It integrates functional contribution tables with `phyloseq` data to explore predicted functions within taxa (taxon-to-function) and taxonomic contributors to functions (function-to-taxon).

[Documentation and tutorials](https://bioscinema.github.io/MicroFun_website_tutorials/) · [Study data](https://bioscinema.github.io/MicroFun_website_tutorials/articles/study-data.html) · [GitHub](https://github.com/bioscinema/MicroFun_0.1.0)

## Graphical abstract

![MicroFun graphical abstract](man/figures/microfun_abstract.png)

## Features

- Taxon-centered and function-centered differential analyses from stratified PICRUSt2 contributions.
- KO, EC, MetaCyc and KEGG pathway annotation and aggregation.
- Class-to-Genus and other supported taxonomy-rank analyses, depending on the input taxonomy.
- Sankey and chord diagrams, boxplots, association dendrograms and effect-size visualizations.
- Community-level PCA, heatmaps and abundance summary plots.

## Installation

```r
install.packages("remotes")
remotes::install_github("bioscinema/MicroFun_0.1.0", dependencies = TRUE)
library(MicroFun)
```

For Bioconductor dependencies and reproducible versions, see the [installation guide](https://bioscinema.github.io/MicroFun_website_tutorials/articles/installation.html).

## Example: Riquelme 2019 demonstration

The complete filtered Riquelme et al. (2019) demonstration is bundled with the package: 68 samples, 2,487 ASVs, four Standard profiles and four full stratified contribution tables. Load only the datasets needed for your analysis. No separate download or KO-part assembly is required.

```r
library(MicroFun)
library(phyloseq)

# Core phyloseq object, Standard tables, selections and provenance.
data("microfun_demo", package = "MicroFun")
ps <- microfun_demo$ps
ko_table <- microfun_demo$standard$KO

# Each full contribution table can be loaded independently.
data("microfun_demo_KO", package = "MicroFun")
# data("microfun_demo_EC", package = "MicroFun")
# data("microfun_demo_KEGG", package = "MicroFun")
# data("microfun_demo_MetaCyc", package = "MicroFun")

ps <- prune_samples(sample_data(ps)$Cohort == "MDA", ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Select project KOs before fitting the function-centered analysis.
selected <- sub("^ko:", "", microfun_demo_KO[["function"]]) %in%
  microfun_demo$selections$KO
contrib <- microfun_demo_KO[selected &
  microfun_demo_KO$sample %in% sample_names(ps), , drop = FALSE]
result_fun <- pathway_sdaa(
  stratified_data = contrib,
  physeq = ps,
  taxon_level = "Genus",
  group = "group",
  direction = "function",
  pathway = "KO",
  daa_method = "Maaslin2",
  reference = "STS",
  p.adjust.method = "BH",
  threshold = 0.05
)
```

Use `?microfun_demo` and `?microfun_demo_KO` for the data structure and provenance. `data()` loads bundled datasets; `load()` reads a standalone R data file. `pathway_sdaa()` accepts exactly one of `stratified_data` (a loaded table) or `stratified_path` (a file path). Full contribution tables can require substantial memory after loading; use `rm(microfun_demo_KO); gc()` when finished.

The source-study data terms and third-party annotation terms must be checked before redistribution. The software MIT license does not establish rights to all derived study data or reference annotations. The original PICRUSt2 version and reference database were not recorded in the supplied files. The dated KEGG table is an alternative snapshot and differs from the package's KO-to-KEGG mapping used by the tutorials.

The [Standard demonstration](https://bioscinema.github.io/MicroFun_website_tutorials/articles/demonstration-standard.html) covers community-level KO and KEGG analyses; the [Stratified demonstration](https://bioscinema.github.io/MicroFun_website_tutorials/articles/demonstration-stratified.html) covers both taxon-centered and function-centered analyses using selected KOs and several visualization methods.

## Authors

Ruitao Liu, Chenyu Liu, Yaoran Yu, Zilin Li and Liangliang Zhang.

Corresponding author: [Liangliang Zhang](mailto:lxz716@case.edu).

Maintainers: [Chenyu Liu](mailto:cxl1299@case.edu), [Yaoran Yu](mailto:yxy1625@case.edu) and [Zilin Li](mailto:zxl1715@case.edu).

## License

MIT. See [LICENSE](LICENSE).
