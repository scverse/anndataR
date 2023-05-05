# anndataR

<!-- README.md is generated from README.qmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/anndataR.png)](https://CRAN.R-project.org/package=anndataR)
<!-- badges: end -->

`{anndataR}` is an R package that brings the power and flexibility of
AnnData to the R ecosystem, allowing you to effortlessly manipulate and
analyze your single-cell data. This package lets you work with backed
h5ad and zarr files, directly access various slots (e.g. X, obs, var,
obsm, obsp), or convert the data into SingleCellExperiment and Seurat
objects.

## Design

This package was initially created at the [scverse 2023-04
hackathon](https://scverse.org/events/2023_04_hackathon/) in Heidelberg.

When fully implemented, it will be a complete replacement for
[theislab/zellkonverter](https://github.com/theislab/zellkonverter),
[mtmorgan/h5ad](github.com/mtmorgan/h5ad/) and
[dynverse/anndata](https://github.com/dynverse/anndata).

## Installation

You can install the development version of `{anndataR}` like so:

``` r
devtools::install_github("scverse/anndataR")
```

## Example

Here’s a quick example of how to use `{anndataR}`:

``` r
library(anndataR)

download.file("https://github.com/chanzuckerberg/cellxgene/raw/main/example-dataset/pbmc3k.h5ad", "pbmc3k.h5ad")

# Read an h5ad file
adata <- read_h5ad("pbmc3k.h5ad")

# Access AnnData slots
adata$X
adata$obs
adata$var

# Convert the AnnData object to a SingleCellExperiment object
sce <- adata$to_SingleCellExperiment()

# Convert the AnnData object to a Seurat object
sce <- adata$to_Seurat()
```
