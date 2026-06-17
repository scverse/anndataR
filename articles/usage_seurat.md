# Read/write Seurat objects using anndataR

## Introduction

This vignette demonstrates how to read and write `Seurat` objects using
the *[anndataR](https://bioconductor.org/packages/3.24/anndataR)*
package, leveraging the interoperability between `Seurat` and the
`AnnData` format.

*[Seurat](https://CRAN.R-project.org/package=Seurat)* is a widely used
toolkit for single-cell analysis in R.
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* enables
conversion between `Seurat` objects and `AnnData` objects, allowing you
to leverage the strengths of both the [scverse](https://scverse.org/)
and [Seurat](https://satijalab.org/seurat/) ecosystems.

### Prerequisites

This vignette requires the
*[Seurat](https://CRAN.R-project.org/package=Seurat)* package in
addition to
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)*. You can
install them using the following code:

``` r

install.packages("Seurat")
```

## Reading H5AD files and Zarr stores to a `Seurat` Object

Using an example `.h5ad` file included in the package, we will
demonstrate how to read an `.h5ad` file and convert it to a `Seurat`
object.

``` r

library(anndataR)
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t

h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
```

Read the `.h5ad` file as a `Seurat` object:

``` r

seurat_obj <- read_h5ad(h5ad_file, as = "Seurat")
seurat_obj
#> An object of class Seurat 
#> 100 features across 50 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  5 layers present: counts, csc_counts, dense_X, dense_counts, X
#>  2 dimensional reductions calculated: X_pca, X_umap
```

This is equivalent to reading in the `.h5ad` file as an `AnnData` and
explicitly converting:

``` r

adata <- read_h5ad(h5ad_file)
seurat_obj <- adata$as_Seurat()
seurat_obj
#> An object of class Seurat 
#> 100 features across 50 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  5 layers present: counts, csc_counts, dense_X, dense_counts, X
#>  2 dimensional reductions calculated: X_pca, X_umap
```

Similarly, we can read from a Zarr store which we also demonstrate with
an example `.zarr` store:

``` r

# Please use "example_v3.zarr.zip" for AnnData stored as Zarr version 3
zarr_path <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(zarr_path, exdir = td)
zarr_path <- file.path(td, "example_v2.zarr")

seurat_obj_zarr <- read_zarr(zarr_path, as = "Seurat")
seurat_obj_zarr
#> An object of class Seurat 
#> 100 features across 50 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  5 layers present: counts, csc_counts, dense_counts, dense_X, X
#>  2 dimensional reductions calculated: X_pca, X_umap
```

or

``` r

adata <- read_zarr(zarr_path)
seurat_obj_zarr <- adata$as_Seurat()
seurat_obj_zarr
#> An object of class Seurat 
#> 100 features across 50 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  5 layers present: counts, csc_counts, dense_counts, dense_X, X
#>  2 dimensional reductions calculated: X_pca, X_umap
```

## Mapping between `AnnData` and `Seurat`

Figure @ref(fig:mapping) shows the structures of the `AnnData` and
`Seurat` objects and how
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* maps
between them. It is important to note that matrices in the two objects
are transposed relative to each other.

![Mapping between \`AnnData\` and \`Seurat\`
objects](../reference/figures/AnnData-Seurat.png)

Mapping between `AnnData` and `Seurat` objects

By default, all items in most slots are converted using the same names.
An exception is the `varp` slot which doesn’t have a corresponding slot
in `Seurat`. Items in the `varm` slot are only converted when they are
specified in a mapping argument. The `Neighbors` and `Images` slots are
not mapped when converting from `Seurat`. See
[`?as_Seurat`](https://anndataR.scverse.org/reference/as_Seurat.md) for
more details on the default mapping.

## Customizing the conversion

You can customize the conversion process by providing specific mappings
for each slot in the `Seurat` object.

Each of the mapping arguments can be provided with one of the following:

- `TRUE`: all items in the slot will be copied using the default mapping
- `FALSE`: the slot will not be copied
- A (named) character vector: the names are the names of the slot in the
  `Seurat` object, the values are the names of the slot in the `AnnData`
  object.

See [`?as_Seurat`](https://anndataR.scverse.org/reference/as_Seurat.md)
for more details on how to customize the conversion process. For
instance:

``` r

seurat_obj <- adata$as_Seurat(
  layers_mapping = c("counts", "dense_counts"),
  object_metadata_mapping = c(metadata1 = "Int", metadata2 = "Float"),
  assay_metadata_mapping = FALSE,
  reduction_mapping = list(
    pca = c(key = "PC_", embeddings = "X_pca", loadings = "PCs"),
    umap = c(key = "UMAP_", embeddings = "X_umap")
  ),
  graph_mapping = TRUE,
  misc_mapping = c(misc1 = "Bool", misc2 = "IntScalar")
)
seurat_obj
#> An object of class Seurat 
#> 100 features across 50 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  3 layers present: counts, dense_counts, X
#>  2 dimensional reductions calculated: pca, umap
```

The mapping arguments can also be passed directly to
[`read_h5ad()`](https://anndataR.scverse.org/reference/read_h5ad.md).

## Writing a `Seurat` object to a H5AD file or Zarr store

The reverse conversion is also possible, allowing you to convert the
`Seurat` object back to an `AnnData` object, or to just write out the
`Seurat` object as an `.h5ad` file or `.zarr` store.

``` r

write_h5ad(seurat_obj, tempfile(fileext = ".h5ad"))
write_zarr(seurat_obj, tempfile(fileext = ".zarr"))
```

This is equivalent to converting the `Seurat` object to an `AnnData`
object and then writing it out:

``` r

adata <- as_AnnData(seurat_obj)
adata$write_h5ad(tempfile(fileext = ".h5ad"))
adata$write_zarr(tempfile(fileext = ".zarr"))
```

You can again customize the conversion process by providing specific
mappings for each slot in the `AnnData` object. For more details, see
[`?as_AnnData`](https://anndataR.scverse.org/reference/as_AnnData.md).

Here’s an example:

``` r

adata <- as_AnnData(
  seurat_obj,
  assay_name = "RNA",
  x_mapping = "counts",
  layers_mapping = c("dense_counts"),
  obs_mapping = c(RNA_count = "nCount_RNA", metadata1 = "metadata1"),
  var_mapping = FALSE,
  obsm_mapping = list(X_pca = "pca", X_umap = "umap"),
  obsp_mapping = TRUE,
  uns_mapping = c("misc1", "misc2")
)
adata
#> InMemoryAnnData object with n_obs × n_vars = 50 × 100
#>     obs: 'RNA_count', 'metadata1'
#>     uns: 'misc1', 'misc2'
#>     obsm: 'X_pca', 'X_umap'
#>     layers: 'dense_counts'
#>     obsp: 'connectivities', 'distances'
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
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Seurat_5.5.0       SeuratObject_5.4.0 sp_2.2-1           anndataR_1.3.0    
#> [5] BiocStyle_2.41.0  
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.5        
#>   [4] spatstat.utils_3.2-3   farver_2.1.2           rmarkdown_2.31        
#>   [7] fs_2.1.0               ragg_1.5.2             vctrs_0.7.3           
#>  [10] ROCR_1.0-12            spatstat.explore_3.8-1 htmltools_0.5.9       
#>  [13] curl_7.1.0             Rhdf5lib_2.1.0         rhdf5_2.57.1          
#>  [16] sass_0.4.10            sctransform_0.4.3      parallelly_1.47.0     
#>  [19] KernSmooth_2.23-26     bslib_0.11.0           htmlwidgets_1.6.4     
#>  [22] desc_1.4.3             ica_1.0-3              httr2_1.2.2           
#>  [25] plyr_1.8.9             plotly_4.12.0          zoo_1.8-15            
#>  [28] cachem_1.1.0           igraph_2.3.2           mime_0.13             
#>  [31] lifecycle_1.0.5        pkgconfig_2.0.3        Matrix_1.7-5          
#>  [34] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6    
#>  [37] future_1.70.0          shiny_1.13.0           digest_0.6.39         
#>  [40] paws.storage_0.10.0    patchwork_1.3.2        tensor_1.5.1          
#>  [43] RSpectra_0.16-2        irlba_2.3.7            textshaping_1.0.5     
#>  [46] progressr_0.19.0       Rarr_2.1.18            spatstat.sparse_3.2-0 
#>  [49] httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
#>  [52] compiler_4.6.0         S7_0.2.2               fastDummies_1.7.6     
#>  [55] grumpy_0.1.1           R.utils_2.13.0         MASS_7.3-65           
#>  [58] rappdirs_0.3.4         tools_4.6.0            lmtest_0.9-40         
#>  [61] otel_0.2.0             httpuv_1.6.17          future.apply_1.20.2   
#>  [64] goftest_1.2-3          R.oo_1.27.1            glue_1.8.1            
#>  [67] nlme_3.1-169           rhdf5filters_1.25.0    promises_1.5.0        
#>  [70] grid_4.6.0             Rtsne_0.17             cluster_2.1.8.2       
#>  [73] reshape2_1.4.5         generics_0.1.4         gtable_0.3.6          
#>  [76] spatstat.data_3.1-9    R.methodsS3_1.8.2      tidyr_1.3.2           
#>  [79] data.table_1.18.4      spatstat.geom_3.8-1    RcppAnnoy_0.0.23      
#>  [82] ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1         
#>  [85] stringr_1.6.0          spam_2.11-4            RcppHNSW_0.7.0        
#>  [88] later_1.4.8            splines_4.6.0          dplyr_1.2.1           
#>  [91] lattice_0.22-9         survival_3.8-6         deldir_2.0-4          
#>  [94] paws.common_0.8.9      tidyselect_1.2.1       miniUI_0.1.2          
#>  [97] pbapply_1.7-4          knitr_1.51             gridExtra_2.3         
#> [100] bookdown_0.47          scattermore_1.2        xfun_0.58             
#> [103] matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.3        
#> [106] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
#> [109] tibble_3.3.1           BiocManager_1.30.27    cli_3.6.6             
#> [112] uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0     
#> [115] systemfonts_1.3.2      jquerylib_0.1.4        Rcpp_1.1.1-1.1        
#> [118] globals_0.19.1         spatstat.random_3.5-0  png_0.1-9             
#> [121] spatstat.univar_3.2-0  parallel_4.6.0         pkgdown_2.2.0         
#> [124] ggplot2_4.0.3          dotCall64_1.2          listenv_0.10.1        
#> [127] viridisLite_0.4.3      scales_1.4.0           ggridges_0.5.7        
#> [130] crayon_1.5.3           purrr_1.2.2            rlang_1.2.0           
#> [133] cowplot_1.2.0
```
