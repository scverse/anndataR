# Read/write SingleCellExperiment objects using anndataR

## Introduction

This vignette demonstrates how to read and write `SingleCellExperiment`
objects using the
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* package,
leveraging the interoperability between `SingleCellExperiment` and the
`AnnData` format.

`SingleCellExperiment` is a widely used class for storing single-cell
data in R, especially within the
[Bioconductor](https://bioconductor.org/) ecosystem.
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* enables
conversion between `SingleCellExperiment` objects and `AnnData` objects,
allowing you to leverage the strengths of both the
[scverse](https://scverse.org/) and
[Bioconductor](https://bioconductor.org/) ecosystems.

### Prerequisites

This vignette requires
*[SingleCellExperiment](https://bioconductor.org/packages/3.24/SingleCellExperiment)*
in addition to
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)*. You can
install them using the following code:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("SingleCellExperiment")
```

## Reading H5AD files and Zarr stores to a `SingleCellExperiment` object

Using an example `.h5ad` file included in the package, we will
demonstrate how to read an `.h5ad` file and convert it to a
`SingleCellExperiment` object.

``` r

library(anndataR)
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from 'package:utils':
#> 
#>     data
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians

h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
```

Read the `.h5ad` file as a `SingleCellExperiment` object:

``` r

sce_obj <- read_h5ad(h5ad_file, as = "SingleCellExperiment")
sce_obj
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(18): Bool BoolNA ... rank_genes_groups umap
#> assays(5): counts csc_counts dense_X dense_counts X
#> rownames(100): Gene000 Gene001 ... Gene098 Gene099
#> rowData names(11): String n_cells_by_counts ... dispersions
#>   dispersions_norm
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(11): Float FloatNA ... log1p_total_counts leiden
#> reducedDimNames(2): X_pca X_umap
#> mainExpName: NULL
#> altExpNames(0):
```

This is equivalent to reading in the `.h5ad` file as an `AnnData` and
explicitly converting:

``` r

adata <- read_h5ad(h5ad_file)
sce <- adata$as_SingleCellExperiment()
sce
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(18): Bool BoolNA ... rank_genes_groups umap
#> assays(5): counts csc_counts dense_X dense_counts X
#> rownames(100): Gene000 Gene001 ... Gene098 Gene099
#> rowData names(11): String n_cells_by_counts ... dispersions
#>   dispersions_norm
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(11): Float FloatNA ... log1p_total_counts leiden
#> reducedDimNames(2): X_pca X_umap
#> mainExpName: NULL
#> altExpNames(0):
```

Similarly, we can read from a Zarr store which we also demonstrate with
an example `.zarr` store:

``` r

# Please use "example_v3.zarr.zip" for AnnData stored as Zarr version 3
zarr_path <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(zarr_path, exdir = td)
zarr_path <- file.path(td, "example_v2.zarr")

sce_zarr <- read_zarr(zarr_path, as = "SingleCellExperiment")
sce_zarr
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(18): Bool BoolNA ... StringScalar umap
#> assays(5): counts csc_counts dense_counts dense_X X
#> rownames(100): Gene000 Gene001 ... Gene098 Gene099
#> rowData names(11): String n_cells_by_counts ... dispersions
#>   dispersions_norm
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(11): Float FloatNA ... log1p_total_counts leiden
#> reducedDimNames(2): X_pca X_umap
#> mainExpName: NULL
#> altExpNames(0):
```

or

``` r

adata <- read_zarr(zarr_path)
sce_zarr <- adata$as_SingleCellExperiment()
sce_zarr
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(18): Bool BoolNA ... StringScalar umap
#> assays(5): counts csc_counts dense_counts dense_X X
#> rownames(100): Gene000 Gene001 ... Gene098 Gene099
#> rowData names(11): String n_cells_by_counts ... dispersions
#>   dispersions_norm
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(11): Float FloatNA ... log1p_total_counts leiden
#> reducedDimNames(2): X_pca X_umap
#> mainExpName: NULL
#> altExpNames(0):
```

## Mapping between `AnnData` and `SingleCellExperiment`

Figure @ref(fig:mapping) shows the structures of the `AnnData` and
`SingleCellExperiment` objects and how
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* maps
between them. It is important to note that matrices in the two objects
are transposed relative to each other.

![Mapping between \`AnnData\` and \`SingleCellExperiment\`
objects](../reference/figures/AnnData-SCE.png)

Mapping between `AnnData` and `SingleCellExperiment` objects

By default, all items in most slots are converted using the same names.
Items in the `varm` slot are only converted when they are specified in a
mapping argument. See
[`?as_SingleCellExperiment`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)
for more details on the default mapping.

## Customizing the conversion

You can customize the conversion process by providing specific mappings
for each slot in the `SingleCellExperiment` object.

Each of the mapping arguments can be provided with one of the following:

- `TRUE`: all items in the slot will be copied using the default mapping
- `FALSE`: the slot will not be copied
- A (named) character vector: the names are the names of the slot in the
  `SingleCellExperiment` object, the values are the names of the slot in
  the `AnnData` object.

See
[`?as_SingleCellExperiment`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)
for more details on how to customize the conversion process. For
instance:

``` r

adata$as_SingleCellExperiment(
  x_mapping = "counts",
  assays_mapping = c("csc_counts"),
  colData_mapping = c("Int", "IntNA"),
  rowData_mapping = c(rowdata1 = "String", rowdata2 = "total_counts"),
  reducedDims_mapping = list(
    "pca" = c(sampleFactors = "X_pca", featureLoadings = "PCs"),
    "umap" = c(sampleFactors = "X_umap")
  ),
  colPairs_mapping = TRUE,
  rowPairs_mapping = FALSE,
  metadata_mapping = c(value1 = "Bool", value2 = "IntScalar")
)
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(2): value1 value2
#> assays(2): counts csc_counts
#> rownames(100): Gene000 Gene001 ... Gene098 Gene099
#> rowData names(2): rowdata1 rowdata2
#> colnames(50): Cell000 Cell001 ... Cell048 Cell049
#> colData names(2): Int IntNA
#> reducedDimNames(2): pca umap
#> mainExpName: NULL
#> altExpNames(0):
```

The mapping arguments can also be passed directly to
[`read_h5ad()`](https://anndataR.scverse.org/reference/read_h5ad.md).

## Writing a `SingleCellExperiment` object to a H5AD file or Zarr store

The reverse conversion is also possible, allowing you to convert a
`SingleCellExperiment` object back to an `AnnData` object, or to just
write out the `SingleCellExperiment` object as an `.h5ad` file or
`.zarr` store.

``` r

write_h5ad(sce_obj, tempfile(fileext = ".h5ad"))
write_zarr(sce_obj, tempfile(fileext = ".zarr"))
```

This is equivalent to converting the `SingleCellExperiment` object to an
`AnnData` object and then writing it out:

``` r

adata <- as_AnnData(sce_obj)
adata$write_h5ad(tempfile(fileext = ".h5ad"))
adata$write_zarr(tempfile(fileext = ".zarr"))
```

You can again customize the conversion process by providing specific
mappings for each slot in the `AnnData` object. For more details, see
[`?as_AnnData`](https://anndataR.scverse.org/reference/as_AnnData.md).

Here’s an example:

``` r

as_AnnData(
  sce_obj,
  x_mapping = "counts",
  layers_mapping = c("csc_counts"),
  obs_mapping = c(metadata1 = "Int", metadata2 = "IntNA"),
  var_mapping = FALSE,
  obsm_mapping = list(X_pca = "X_pca", X_umap = "X_umap"),
  obsp_mapping = TRUE,
  uns_mapping = c("Bool", "IntScalar")
)
#> InMemoryAnnData object with n_obs × n_vars = 50 × 100
#>     obs: 'metadata1', 'metadata2'
#>     uns: 'Bool', 'IntScalar'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'X_pca'
#>     layers: 'csc_counts'
#>     obsp: 'connectivities', 'distances'
#>     varp: 'test_varp'
```

The mapping arguments can also be passed directly to
[`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)
or
[`write_zarr()`](https://anndataR.scverse.org/reference/write_zarr.md).

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SingleCellExperiment_1.35.1 SummarizedExperiment_1.43.0
#>  [3] Biobase_2.73.1              GenomicRanges_1.65.0       
#>  [5] Seqinfo_1.3.0               IRanges_2.47.2             
#>  [7] S4Vectors_0.51.3            BiocGenerics_0.59.7        
#>  [9] generics_0.1.4              MatrixGenerics_1.25.0      
#> [11] matrixStats_1.5.0           anndataR_1.3.0             
#> [13] BiocStyle_2.41.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] xfun_0.59           bslib_0.11.0        httr2_1.2.2        
#>  [4] htmlwidgets_1.6.4   rhdf5_2.57.1        lattice_0.22-9     
#>  [7] rhdf5filters_1.25.0 vctrs_0.7.3         tools_4.6.0        
#> [10] curl_7.1.0          paws.common_0.8.10  R.oo_1.27.1        
#> [13] Matrix_1.7-5        desc_1.4.3          lifecycle_1.0.5    
#> [16] compiler_4.6.0      textshaping_1.0.5   htmltools_0.5.9    
#> [19] sass_0.4.10         yaml_2.3.12         pkgdown_2.2.0      
#> [22] crayon_1.5.3        jquerylib_0.1.4     R.utils_2.13.0     
#> [25] grumpy_0.1.1        DelayedArray_0.39.3 cachem_1.1.0       
#> [28] abind_1.4-8         digest_0.6.39       paws.storage_0.10.0
#> [31] purrr_1.2.2         bookdown_0.47       fastmap_1.2.0      
#> [34] grid_4.6.0          cli_3.6.6           SparseArray_1.13.2 
#> [37] Rarr_2.1.18         magrittr_2.0.5      S4Arrays_1.13.0    
#> [40] withr_3.0.3         rappdirs_0.3.4      rmarkdown_2.31     
#> [43] XVector_0.53.0      otel_0.2.0          reticulate_1.46.0  
#> [46] ragg_1.5.2          png_0.1-9           R.methodsS3_1.8.2  
#> [49] evaluate_1.0.5      knitr_1.51          rlang_1.2.0        
#> [52] Rcpp_1.1.1-1.1      glue_1.8.1          BiocManager_1.30.27
#> [55] jsonlite_2.0.0      R6_2.6.1            Rhdf5lib_2.1.0     
#> [58] systemfonts_1.3.2   fs_2.1.0
```
