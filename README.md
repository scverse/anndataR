# anndataR


<!-- README.md is generated from README.qmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/anndataR.png)](https://CRAN.R-project.org/package=anndataR)
<!-- badges: end -->

`{anndataR}` aims to make the AnnData format a first-class citizen in
the R ecosystem, and to make it easy to work with AnnData files in R,
either directly or by converting it to a SingleCellExperiment or Seurat
object.

Feature list:

- Provide an `R6` class to work with AnnData objects in R (either
  in-memory or on-disk).
- Read/write `*.h5ad` files natively
- Convert to/from `SingleCellExperiment` objects
- Convert to/from `Seurat` objects

> [!WARNING]
>
> This package is still in the experimental stage, and may not work as
> expected. You can find the status of development of anndataR on the
> [feature tracking
> page](https://anndatar.data-intuitive.com/articles/design.html#feature-tracking)
> of the website.Please report any issues you encounter.

## Installation

You can install the development version of `{anndataR}` like so:

``` r
devtools::install_github("scverse/anndataR")
```

You might need to install suggested dependencies manually, depending on
the task you want to perform.

- To read/write \*.h5ad files, you need to install
  [hdf5r](https://cran.r-project.org/package=hdf5r):  
  `BiocManager::install("hdf5r")`
- To read/write \*.zarr files, you need to install
  [zarr](https://github.com/keller-mark/pizzarr):  
  `devtools::install_github("keller-mark/pizzarr")`
- To convert to/from `SingleCellExperiment` objects, you need to install
  [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html):  
  `BiocManager::install("SingleCellExperiment")`
- To convert to/from `Seurat` objects, you need to install
  [SeuratObject](https://cran.r-project.org/package=SeuratObject):  
  `install.packages("SeuratObject")`

If you’re feeling adventurous, you can install all suggested
dependencies at once:

``` r
devtools::install_github("scverse/anndataR", dependencies = TRUE)
```

## Example

Here’s a quick example of how to use `{anndataR}`. First, we download an
h5ad file.

``` r
library(anndataR)

h5ad_path <- system.file("extdata", "example.h5ad", package = "anndataR")
```

Read an h5ad file in memory:

``` r
adata <- read_h5ad(h5ad_path)
```

Read an h5ad file on disk:

``` r
adata <- read_h5ad(h5ad_path, to = "HDF5AnnData")
```

View structure:

``` r
adata
#> AnnData object with n_obs × n_vars = 50 × 100
#>     obs: 'Float', 'FloatNA', 'Int', 'IntNA', 'Bool', 'BoolNA', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'leiden'
#>     var: 'String', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'PCs'
#>     layers: 'counts', 'csc_counts', 'dense_X', 'dense_counts'
#>     obsp: 'connectivities', 'distances'
```

Access AnnData slots:

``` r
dim(adata$X)
#> [1]  50 100
adata$obs[1:5, 1:6]
#>         Float FloatNA Int IntNA  Bool BoolNA
#> Cell000 42.42     NaN   0    NA FALSE  FALSE
#> Cell001 42.42   42.42   1    42  TRUE     NA
#> Cell002 42.42   42.42   2    42  TRUE   TRUE
#> Cell003 42.42   42.42   3    42  TRUE   TRUE
#> Cell004 42.42   42.42   4    42  TRUE   TRUE
adata$var[1:5, 1:6]
#>          String n_cells_by_counts mean_counts log1p_mean_counts pct_dropout_by_counts total_counts
#> Gene000 String0                44        1.94          1.078410                    12           97
#> Gene001 String1                42        2.04          1.111858                    16          102
#> Gene002 String2                43        2.12          1.137833                    14          106
#> Gene003 String3                41        1.72          1.000632                    18           86
#> Gene004 String4                42        2.06          1.118415                    16          103
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
#> rowData names(11): String n_cells_by_counts ... dispersions dispersions_norm
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(11): Float FloatNA ... log1p_total_counts leiden
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Convert the AnnData object to a Seurat object:

``` r
obj <- adata$to_Seurat()
obj
#> An object of class Seurat 
#> 500 features across 50 samples within 5 assays 
#> Active assay: RNA (100 features, 0 variable features)
#>  2 layers present: counts, data
#>  4 other assays present: counts, csc_counts, dense_X, dense_counts
```

## Manually create an object

``` r
adata <- AnnData(
  X = matrix(rnorm(100), nrow = 10),
  obs = data.frame(
    cell_type = factor(rep(c("A", "B"), each = 5))
  ),
  var = data.frame(
    gene_name = paste0("gene_", 1:10)
  )
)

adata
#> AnnData object with n_obs × n_vars = 10 × 10
#>     obs: 'cell_type'
#>     var: 'gene_name'
```
