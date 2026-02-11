# Python Integration with anndataR

## Introduction

*[anndataR](https://bioconductor.org/packages/3.22/anndataR)* works with
Python `AnnData` objects through
*[reticulate](https://CRAN.R-project.org/package=reticulate)*. You can
load Python objects, apply Python functions to them, and convert them to
`Seurat` or `SingleCellExperiment` objects.

``` r
message(
  "Python packages scanpy and mudata are required to run this vignette. Code chunks will not be evaluated."
)
#> Python packages scanpy and mudata are required to run this vignette. Code chunks will not be evaluated.
```

### Prerequisites

This vignette requires Python with the *scanpy* and *mudata* packages
installed. If these are not available, the code chunks will not be
evaluated but the examples remain visible.

## Basic Integration with *scanpy*

Install required Python packages if needed:

``` r
reticulate::py_require("scanpy")
```

``` r
library(anndataR)
library(reticulate)
sc <- import("scanpy")
```

Load a dataset directly from *scanpy*:

``` r
adata <- sc$datasets$pbmc3k_processed()
adata
```

Apply *scanpy* functions directly:

``` r
sc$pp$filter_cells(adata, min_genes = 200L)
sc$pp$normalize_total(adata, target_sum = 1e4)
sc$pp$log1p(adata)
```

## Conversion to R objects

Convert to `SingleCellExperiment` (see
[`vignette("usage_singlecellexperiment")`](https://anndataR.data-intuitive.com/articles/usage_singlecellexperiment.md)):

``` r
sce_obj <- adata$as_SingleCellExperiment()
sce_obj
```

Convert to `Seurat` (see
[`vignette("usage_seurat")`](https://anndataR.data-intuitive.com/articles/usage_seurat.md)):

``` r
seurat_obj <- adata$as_Seurat()
seurat_obj
```

## Multi-modal data with *mudata*

Install required Python packages if needed:

``` r
reticulate::py_install("mudata")
```

``` r
md <- import("mudata")
```

Load a `MuData` object from file:

``` r
cache <- BiocFileCache::BiocFileCache(ask = FALSE)
h5mu_file <- BiocFileCache::bfcrpath(
  cache,
  "https://github.com/gtca/h5xx-datasets/raw/b1177ac8877c89d8bb355b072164384b4e9cc81d/datasets/minipbcite.h5mu"
)

mdata <- md$read_h5mu(h5mu_file)
```

Access individual modalities and convert them:

``` r
rna_mod <- mdata$mod[["rna"]]

rna_seurat <- rna_mod$as_Seurat()
print(rna_seurat)

rna_sce <- rna_mod$as_SingleCellExperiment()
print(rna_sce)
```

## Session info

### R

``` r
sessionInfo()
```

### Python

``` r
reticulate::py_config()

reticulate::py_list_packages()
```
