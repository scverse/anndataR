# Python Integration with anndataR

## Introduction

*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* works with
Python `AnnData` objects through
*[reticulate](https://CRAN.R-project.org/package=reticulate)*. You can
load Python objects, apply Python functions to them, and convert them to
`Seurat` or `SingleCellExperiment` objects.

``` r

message(
  "Python packages scanpy and mudata are required to run this vignette. Code chunks will not be evaluated."
)
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
#> ReticulateAnnData object with n_obs × n_vars = 2638 × 1838
#>     obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain'
#>     var: 'n_cells'
#>     uns: 'draw_graph', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
#>     obsm: 'X_pca', 'X_tsne', 'X_umap', 'X_draw_graph_fr'
#>     varm: 'PCs'
#>     obsp: 'distances', 'connectivities'
```

Apply *scanpy* functions directly:

``` r

sc$pp$filter_cells(adata, min_genes = 200L)
sc$pp$normalize_total(adata, target_sum = 1e4)
sc$pp$log1p(adata)
```

## Conversion to R objects

Convert to `SingleCellExperiment` (see
[`vignette("usage_singlecellexperiment")`](https://anndataR.scverse.org/articles/usage_singlecellexperiment.md)):

``` r

sce_obj <- adata$as_SingleCellExperiment()
sce_obj
#> class: SingleCellExperiment 
#> dim: 1838 538 
#> metadata(7): draw_graph louvain ... rank_genes_groups log1p
#> assays(1): X
#> rownames(1838): TNFRSF4 CPSF3L ... S100B PRMT2
#> rowData names(1): n_cells
#> colnames(538): AAACATTGAGCTAC-1 AAACATTGATCAGC-1 ... TTTCGAACTCTCAT-1
#>   TTTCTACTGAGGCA-1
#> colData names(4): n_genes percent_mito n_counts louvain
#> reducedDimNames(4): X_pca X_tsne X_umap X_draw_graph_fr
#> mainExpName: NULL
#> altExpNames(0):
```

Convert to `Seurat` (see
[`vignette("usage_seurat")`](https://anndataR.scverse.org/articles/usage_seurat.md)):

``` r

seurat_obj <- adata$as_Seurat()
#> Warning: No "counts" or "data" layer found in `names(layers_mapping)`, this may lead to
#> unexpected results when using the resulting <Seurat> object.
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
seurat_obj
#> An object of class Seurat 
#> 1838 features across 538 samples within 1 assay 
#> Active assay: RNA (1838 features, 0 variable features)
#>  1 layer present: X
#>  4 dimensional reductions calculated: X_pca, X_tsne, X_umap, X_draw_graph_fr
```

## Multi-modal data with *mudata*

Install required Python packages if needed:

``` r

reticulate::py_install("mudata<0.4")
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
#> adding rname 'https://github.com/gtca/h5xx-datasets/raw/b1177ac8877c89d8bb355b072164384b4e9cc81d/datasets/minipbcite.h5mu'

mdata <- md$read_h5mu(h5mu_file)
```

Access individual modalities and convert them:

``` r

rna_mod <- mdata$mod[["rna"]]

rna_seurat <- rna_mod$as_Seurat()
#> Warning: No "counts" or "data" layer found in `names(layers_mapping)`, this may lead to
#> unexpected results when using the resulting <Seurat> object.
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
print(rna_seurat)
#> An object of class Seurat 
#> 27 features across 411 samples within 1 assay 
#> Active assay: RNA (27 features, 0 variable features)
#>  1 layer present: X
#>  2 dimensional reductions calculated: X_pca, X_umap

rna_sce <- rna_mod$as_SingleCellExperiment()
print(rna_sce)
#> class: SingleCellExperiment 
#> dim: 27 411 
#> metadata(8): celltype_colors hvg ... rank_genes_groups umap
#> assays(1): X
#> rownames(27): NKG7 KLRC2 ... MS4A1 KLF4
#> rowData names(13): gene_ids feature_types ... mean std
#> colnames(411): CAGCCAGGTCTCGACG-1 TTCTTCCTCTCGGTAA-1 ...
#>   GACTCTCCAGCTCTGG-1 GAAATGACAAGCACCC-1
#> colData names(6): n_genes_by_counts total_counts ... leiden celltype
#> reducedDimNames(2): X_pca X_umap
#> mainExpName: NULL
#> altExpNames(0):
```

## Session info

### R

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
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
#> [1] reticulate_1.47.0 anndataR_1.3.4    BiocStyle_2.41.0 
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.3.0                   pbapply_1.7-5              
#>   [3] httr2_1.3.0                 gridExtra_2.3.1            
#>   [5] rlang_1.3.0                 magrittr_2.0.5             
#>   [7] RcppAnnoy_0.0.23            otel_0.3.0                 
#>   [9] RSQLite_3.53.3              matrixStats_1.5.0          
#>  [11] ggridges_0.5.7              compiler_4.6.2             
#>  [13] reshape2_1.4.5              png_0.1-9                  
#>  [15] systemfonts_1.3.2           vctrs_0.7.3                
#>  [17] stringr_1.6.0               pkgconfig_2.0.3            
#>  [19] fastmap_1.2.0               dbplyr_2.6.0               
#>  [21] XVector_0.53.0              rmarkdown_2.32             
#>  [23] ragg_1.5.2                  bit_4.6.0                  
#>  [25] purrr_1.2.2                 xfun_0.62                  
#>  [27] cachem_1.1.0                jsonlite_2.0.0             
#>  [29] blob_1.3.0                  DelayedArray_0.39.8        
#>  [31] irlba_2.4.1                 parallel_4.6.2             
#>  [33] cluster_2.1.8.3             R6_2.6.1                   
#>  [35] ica_1.0-3                   stringi_1.8.9              
#>  [37] bslib_0.12.0                RColorBrewer_1.1-3         
#>  [39] parallelly_1.50.0           GenomicRanges_1.65.4       
#>  [41] jquerylib_0.1.4             scattermore_1.2            
#>  [43] Rcpp_1.1.2                  Seqinfo_1.3.2              
#>  [45] bookdown_0.48               SummarizedExperiment_1.43.0
#>  [47] knitr_1.52                  future.apply_1.21.0        
#>  [49] sctransform_0.4.3           IRanges_2.47.5             
#>  [51] Matrix_1.7-6                igraph_2.3.4               
#>  [53] tidyselect_1.2.1            abind_1.4-8                
#>  [55] yaml_2.3.12                 codetools_0.2-20           
#>  [57] curl_8.1.0                  listenv_1.1.0              
#>  [59] plyr_1.8.9                  lattice_0.23-1             
#>  [61] tibble_3.3.1                withr_3.0.3                
#>  [63] Biobase_2.73.2              S7_0.2.2                   
#>  [65] evaluate_1.0.5              Rtsne_0.17                 
#>  [67] future_1.82.0               fastDummies_1.7.6          
#>  [69] desc_1.4.3                  BiocFileCache_3.3.2        
#>  [71] filelock_1.0.3              pillar_1.11.1              
#>  [73] BiocManager_1.30.27         MatrixGenerics_1.25.0      
#>  [75] Seurat_6.0.0                KernSmooth_2.23-27         
#>  [77] stats4_4.6.2                generics_0.1.4             
#>  [79] RcppHNSW_0.7.0              sp_2.2-4                   
#>  [81] S4Vectors_0.51.10           ggplot2_4.0.3              
#>  [83] scales_1.4.0                globals_0.20.0             
#>  [85] glue_1.8.1                  tools_4.6.2                
#>  [87] RSpectra_0.16-2             RANN_2.6.3                 
#>  [89] fs_2.1.0                    dotCall64_1.2              
#>  [91] cowplot_1.2.0               grid_4.6.2                 
#>  [93] SingleCellExperiment_1.35.2 patchwork_1.3.2            
#>  [95] cli_3.6.6                   textshaping_1.0.5          
#>  [97] spam_2.11-4                 S4Arrays_1.13.2            
#>  [99] dplyr_1.2.1                 uwot_0.2.5                 
#> [101] gtable_0.3.6                pls_2.9-0                  
#> [103] sass_0.4.10                 digest_0.6.39              
#> [105] progressr_1.1.0             BiocGenerics_0.59.12       
#> [107] SparseArray_1.13.4          ggrepel_0.9.8              
#> [109] htmlwidgets_1.6.4           SeuratObject_5.4.0         
#> [111] farver_2.1.2                memoise_2.0.1              
#> [113] htmltools_0.5.9             pkgdown_2.2.1              
#> [115] lifecycle_1.0.5             httr_1.4.9                 
#> [117] bit64_4.8.6                 MASS_7.3-66
```

### Python

``` r

reticulate::py_config()
#> python:         /opt/hostedtoolcache/Python/3.14.10/x64/bin/python3
#> libpython:      /opt/hostedtoolcache/Python/3.14.10/x64/lib/libpython3.14.so
#> pythonhome:     /opt/hostedtoolcache/Python/3.14.10/x64:/opt/hostedtoolcache/Python/3.14.10/x64
#> version:        3.14.10 (main, Oct 26 2026, 14:17:03) [GCC 15.2.0]
#> numpy:          /opt/hostedtoolcache/Python/3.14.10/x64/lib/python3.14/site-packages/numpy
#> numpy_version:  2.5.4
#> scanpy:         /opt/hostedtoolcache/Python/3.14.10/x64/lib/python3.14/site-packages/scanpy
#> 
#> NOTE: Python version was forced by RETICULATE_PYTHON

reticulate::py_list_packages()
#>              package     version               requirement
#> 1            anndata      0.13.4           anndata==0.13.4
#> 2    annotated-types       0.8.0    annotated-types==0.8.0
#> 3   array-api-compat      1.16.0  array-api-compat==1.16.0
#> 4            certifi  2026.10.23       certifi==2026.10.23
#> 5        cloudpickle       3.1.2        cloudpickle==3.1.2
#> 6          contourpy       1.4.0          contourpy==1.4.0
#> 7             cycler      0.12.1            cycler==0.12.1
#> 8             donfig 0.8.1.post1       donfig==0.8.1.post1
#> 9   fast-array-utils       1.5.1   fast-array-utils==1.5.1
#> 10         fonttools      4.68.0         fonttools==4.68.0
#> 11         formulaic       1.2.2          formulaic==1.2.2
#> 12     google-crc32c       1.9.0      google-crc32c==1.9.0
#> 13              h5py      3.16.0              h5py==3.16.0
#> 14    interface_meta       2.0.1     interface_meta==2.0.1
#> 15            joblib       1.6.0             joblib==1.6.0
#> 16        kiwisolver       1.5.1         kiwisolver==1.5.1
#> 17   legacy-api-wrap         1.5      legacy-api-wrap==1.5
#> 18          llvmlite      0.50.0          llvmlite==0.50.0
#> 19        matplotlib      3.11.2        matplotlib==3.11.2
#> 20           msgspec      0.22.0           msgspec==0.22.0
#> 21            mudata      0.3.10            mudata==0.3.10
#> 22          narwhals      2.28.0          narwhals==2.28.0
#> 23           natsort       8.4.0            natsort==8.4.0
#> 24          networkx         3.7             networkx==3.7
#> 25             numba      0.68.1             numba==0.68.1
#> 26         numcodecs      0.17.1         numcodecs==0.17.1
#> 27             numpy       2.5.4              numpy==2.5.4
#> 28         packaging        26.3           packaging==26.3
#> 29            pandas       3.1.1             pandas==3.1.1
#> 30             patsy       1.0.3              patsy==1.0.3
#> 31            pillow      13.0.0            pillow==13.0.0
#> 32          pydantic      2.14.3          pydantic==2.14.3
#> 33 pydantic-settings      2.16.2 pydantic-settings==2.16.2
#> 34     pydantic_core      2.50.4     pydantic_core==2.50.4
#> 35       pynndescent       0.7.0        pynndescent==0.7.0
#> 36         pyparsing       3.3.3          pyparsing==3.3.3
#> 37   python-dateutil       4.0.0    python-dateutil==4.0.0
#> 38     python-dotenv       1.2.4      python-dotenv==1.2.4
#> 39            PyYAML       6.0.3             PyYAML==6.0.3
#> 40            scanpy      1.12.5            scanpy==1.12.5
#> 41      scikit-learn       1.9.1       scikit-learn==1.9.1
#> 42             scipy      1.18.1             scipy==1.18.1
#> 43      scverse-misc       0.1.1       scverse-misc==0.1.1
#> 44           seaborn      0.13.2           seaborn==0.13.2
#> 45     session-info2       0.4.2      session-info2==0.4.2
#> 46       statsmodels      0.15.0       statsmodels==0.15.0
#> 47     threadpoolctl       3.8.0      threadpoolctl==3.8.0
#> 48              tqdm      4.70.1              tqdm==4.70.1
#> 49 typing-inspection       0.4.5  typing-inspection==0.4.5
#> 50 typing_extensions      4.16.0 typing_extensions==4.16.0
#> 51            tzdata      2026.6            tzdata==2026.6
#> 52        umap-learn       0.6.0         umap-learn==0.6.0
#> 53             wrapt  2.6.0.dev1         wrapt==2.6.0.dev1
#> 54              zarr       3.4.2               zarr==3.4.2
```
