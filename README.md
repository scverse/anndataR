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

Here’s a quick example of how to use `{anndataR}`. First, we download an
h5ad file.

``` r
library(anndataR)

h5ad_path <- system.file("extdata", "example.h5ad", package = "anndataR")
```

Read an h5ad file:

``` r
adata <- read_h5ad(h5ad_path, to = "InMemoryAnnData")
```

View structure:

``` r
adata
#> class: InMemoryAnnData
#> dim: 50 obs x 100 var
#> X: dgRMatrix
#> layers: counts csc_counts dense_X dense_counts
#> obs: Float FloatNA Int IntNA Bool BoolNA n_genes_by_counts
#>   log1p_n_genes_by_counts total_counts log1p_total_counts leiden
#> var: String n_cells_by_counts mean_counts log1p_mean_counts
#>   pct_dropout_by_counts total_counts log1p_total_counts highly_variable
#>   means dispersions dispersions_norm
```

Access AnnData slots:

``` r
dim(adata$X)
#> [1]  50 100
adata$obs[1:5, 1:6]
#>   Float FloatNA Int IntNA  Bool BoolNA
#> 1 42.42     NaN   0    NA FALSE  FALSE
#> 2 42.42   42.42   1    42  TRUE     NA
#> 3 42.42   42.42   2    42  TRUE   TRUE
#> 4 42.42   42.42   3    42  TRUE   TRUE
#> 5 42.42   42.42   4    42  TRUE   TRUE
adata$var[1:5, 1:6]
#>    String n_cells_by_counts mean_counts log1p_mean_counts pct_dropout_by_counts
#> 1 String0                44        1.94          1.078410                    12
#> 2 String1                42        2.04          1.111858                    16
#> 3 String2                43        2.12          1.137833                    14
#> 4 String3                41        1.72          1.000632                    18
#> 5 String4                42        2.06          1.118415                    16
#>   total_counts
#> 1           97
#> 2          102
#> 3          106
#> 4           86
#> 5          103
```

## Interoperability

Convert the AnnData object to a SingleCellExperiment object:

``` r
sce <- adata$to_SingleCellExperiment()
sce
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(0):
#> assays(5): X counts csc_counts dense_X dense_counts
#> rownames(100): Gene000 Gene001 ... Gene098 Gene099
#> rowData names(11): String n_cells_by_counts ... dispersions
#>   dispersions_norm
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(11): Float FloatNA ... log1p_total_counts leiden
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Convert the AnnData object to a Seurat object:

``` r
obj <- adata$to_Seurat()
#> Loading required namespace: SeuratObject
#> The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
#> which was just loaded, will retire in October 2023.
#> Please refer to R-spatial evolution reports for details, especially
#> https://r-spatial.org/r/2023/05/15/evolution4.html.
#> It may be desirable to make the sf package available;
#> package maintainers should consider adding sf to Suggests:.
#> The sp package is now running under evolution status 2
#>      (status 2 uses the sf package in place of rgdal)
#> Warning: Keys should be one or more alphanumeric characters followed by an
#> underscore, setting key from rna to rna_
#> Warning: Keys should be one or more alphanumeric characters followed by an
#> underscore, setting key from csc_counts_ to csccounts_
#> Warning: Keys should be one or more alphanumeric characters followed by an
#> underscore, setting key from dense_x_ to densex_
#> Warning: Keys should be one or more alphanumeric characters followed by an
#> underscore, setting key from dense_counts_ to densecounts_
obj
#> An object of class Seurat 
#> 500 features across 50 samples within 5 assays 
#> Active assay: RNA (100 features, 0 variable features)
#>  4 other assays present: counts, csc_counts, dense_X, dense_counts
```
