## Introduction

MicroFun is a lightweight but powerful analysis tool for statistical inference and visualizing PICRUSt2 functional outputs. By leveraging stratified PICRUSt2 outputs, MicroFun fills a critical gap between taxonomic and functional analyses, explicitly linking differential taxa to their predicted functional contributions. This extends the capabilities of ggpicrust2, which is unable to determine the taxonomic sources of predicted functions. The tool is implemented as an R package.

<img width="9006" height="3856" alt="microfun_abstract_1022" src="https://github.com/user-attachments/assets/59d1298c-ecd6-447f-bf89-0e0c54b58cf5" />

## Installation

```r
# install devtools if needed
install.packages("devtools")

# install from GitHub
devtools::install_github("bioscinema/MicroFun")
```

## Functions overview

- `tax2fun_sankey()`: Sankey diagram from taxonomy to function

- `fun2tax_sankey()`: Sankey diagram from function to taxonomy

- `tax2fun_box()`: Box plot of function contributions across group for every taxonomy

- `fun2tax_box()`: Box plot of taxonomy abundance across group for every function

- `tax2fun_dendro()`: Dendrogram for the relationship between taxonomy and functions

- `phy2fun_dendro()`: Tree plot connecting phylogenetic tree with associated functions

## Website

We provide a detailed tutorial on [https://bioscinema.github.io/MicroFun_website/]. If you have any issues, please post your questions in the Issue section and we will answer the questions as soon as possible!
