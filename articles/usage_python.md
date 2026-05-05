# Python Integration with anndataR

## Introduction

*[anndataR](https://bioconductor.org/packages/3.23/anndataR)* works with
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
#> R Under development (unstable) (2026-05-03 r89994)
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
#> [1] reticulate_1.46.0 anndataR_1.3.0    BiocStyle_2.40.0 
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] magrittr_2.0.5              spatstat.utils_3.2-2       
#>   [5] farver_2.1.2                rmarkdown_2.31             
#>   [7] fs_2.1.0                    ragg_1.5.2                 
#>   [9] vctrs_0.7.3                 ROCR_1.0-12                
#>  [11] memoise_2.0.1               spatstat.explore_3.8-0     
#>  [13] htmltools_0.5.9             S4Arrays_1.12.0            
#>  [15] curl_7.1.0                  SparseArray_1.12.2         
#>  [17] sass_0.4.10                 sctransform_0.4.3          
#>  [19] parallelly_1.47.0           KernSmooth_2.23-26         
#>  [21] bslib_0.10.0                htmlwidgets_1.6.4          
#>  [23] desc_1.4.3                  ica_1.0-3                  
#>  [25] httr2_1.2.2                 plyr_1.8.9                 
#>  [27] plotly_4.12.0               zoo_1.8-15                 
#>  [29] cachem_1.1.0                igraph_2.3.1               
#>  [31] mime_0.13                   lifecycle_1.0.5            
#>  [33] pkgconfig_2.0.3             Matrix_1.7-5               
#>  [35] R6_2.6.1                    fastmap_1.2.0              
#>  [37] MatrixGenerics_1.24.0       fitdistrplus_1.2-6         
#>  [39] future_1.70.0               shiny_1.13.0               
#>  [41] digest_0.6.39               patchwork_1.3.2            
#>  [43] S4Vectors_0.50.0            Seurat_5.5.0               
#>  [45] tensor_1.5.1                RSpectra_0.16-2            
#>  [47] irlba_2.3.7                 RSQLite_2.4.6              
#>  [49] textshaping_1.0.5           GenomicRanges_1.64.0       
#>  [51] filelock_1.0.3              progressr_0.19.0           
#>  [53] spatstat.sparse_3.1-0       polyclip_1.10-7            
#>  [55] httr_1.4.8                  abind_1.4-8                
#>  [57] compiler_4.7.0              withr_3.0.2                
#>  [59] bit64_4.8.0                 S7_0.2.2                   
#>  [61] DBI_1.3.0                   fastDummies_1.7.6          
#>  [63] MASS_7.3-65                 rappdirs_0.3.4             
#>  [65] DelayedArray_0.38.1         tools_4.7.0                
#>  [67] lmtest_0.9-40               otel_0.2.0                 
#>  [69] httpuv_1.6.17               future.apply_1.20.2        
#>  [71] goftest_1.2-3               glue_1.8.1                 
#>  [73] nlme_3.1-169                promises_1.5.0             
#>  [75] grid_4.7.0                  Rtsne_0.17                 
#>  [77] cluster_2.1.8.2             reshape2_1.4.5             
#>  [79] generics_0.1.4              gtable_0.3.6               
#>  [81] spatstat.data_3.1-9         tidyr_1.3.2                
#>  [83] data.table_1.18.2.1         sp_2.2-1                   
#>  [85] XVector_0.52.0              BiocGenerics_0.58.0        
#>  [87] spatstat.geom_3.7-3         RcppAnnoy_0.0.23           
#>  [89] ggrepel_0.9.8               RANN_2.6.2                 
#>  [91] pillar_1.11.1               stringr_1.6.0              
#>  [93] spam_2.11-3                 RcppHNSW_0.6.0             
#>  [95] later_1.4.8                 splines_4.7.0              
#>  [97] dplyr_1.2.1                 BiocFileCache_3.2.0        
#>  [99] lattice_0.22-9              bit_4.6.0                  
#> [101] deldir_2.0-4                survival_3.8-6             
#> [103] tidyselect_1.2.1            SingleCellExperiment_1.34.0
#> [105] miniUI_0.1.2                pbapply_1.7-4              
#> [107] knitr_1.51                  gridExtra_2.3              
#> [109] bookdown_0.46               IRanges_2.46.0             
#> [111] Seqinfo_1.2.0               SummarizedExperiment_1.42.0
#> [113] scattermore_1.2             stats4_4.7.0               
#> [115] xfun_0.57                   Biobase_2.72.0             
#> [117] matrixStats_1.5.0           stringi_1.8.7              
#> [119] lazyeval_0.2.3              yaml_2.3.12                
#> [121] evaluate_1.0.5              codetools_0.2-20           
#> [123] tibble_3.3.1                BiocManager_1.30.27        
#> [125] cli_3.6.6                   uwot_0.2.4                 
#> [127] xtable_1.8-8                systemfonts_1.3.2          
#> [129] jquerylib_0.1.4             Rcpp_1.1.1-1.1             
#> [131] spatstat.random_3.4-5       globals_0.19.1             
#> [133] dbplyr_2.5.2                png_0.1-9                  
#> [135] spatstat.univar_3.1-7       parallel_4.7.0             
#> [137] blob_1.3.0                  pkgdown_2.2.0              
#> [139] ggplot2_4.0.3               dotCall64_1.2              
#> [141] listenv_0.10.1              viridisLite_0.4.3          
#> [143] scales_1.4.0                ggridges_0.5.7             
#> [145] SeuratObject_5.4.0          purrr_1.2.2                
#> [147] rlang_1.2.0                 cowplot_1.2.0
```

### Python

``` r

reticulate::py_config()
#> python:         /opt/hostedtoolcache/Python/3.14.4/x64/bin/python3
#> libpython:      /opt/hostedtoolcache/Python/3.14.4/x64/lib/libpython3.14.so
#> pythonhome:     /opt/hostedtoolcache/Python/3.14.4/x64:/opt/hostedtoolcache/Python/3.14.4/x64
#> version:        3.14.4 (main, Apr  8 2026, 02:27:22) [GCC 13.3.0]
#> numpy:          /opt/hostedtoolcache/Python/3.14.4/x64/lib/python3.14/site-packages/numpy
#> numpy_version:  2.4.4
#> scanpy:         /opt/hostedtoolcache/Python/3.14.4/x64/lib/python3.14/site-packages/scanpy
#> 
#> NOTE: Python version was forced by RETICULATE_PYTHON

reticulate::py_list_packages()
#>              package     version                  requirement
#> 1            anndata     0.12.11             anndata==0.12.11
#> 2   array-api-compat      1.14.0     array-api-compat==1.14.0
#> 3            certifi   2026.4.22           certifi==2026.4.22
#> 4          contourpy       1.3.3             contourpy==1.3.3
#> 5             cycler      0.12.1               cycler==0.12.1
#> 6             donfig 0.8.1.post1          donfig==0.8.1.post1
#> 7   fast-array-utils       1.4.1      fast-array-utils==1.4.1
#> 8          fonttools      4.62.1            fonttools==4.62.1
#> 9      google-crc32c       1.8.0         google-crc32c==1.8.0
#> 10              h5py      3.16.0                 h5py==3.16.0
#> 11            joblib       1.5.3                joblib==1.5.3
#> 12        kiwisolver       1.5.0            kiwisolver==1.5.0
#> 13   legacy-api-wrap         1.5         legacy-api-wrap==1.5
#> 14          llvmlite      0.47.0             llvmlite==0.47.0
#> 15        matplotlib      3.10.9           matplotlib==3.10.9
#> 16            mudata       0.3.7                mudata==0.3.7
#> 17           natsort       8.4.0               natsort==8.4.0
#> 18          networkx       3.6.1              networkx==3.6.1
#> 19             numba      0.65.1                numba==0.65.1
#> 20         numcodecs      0.16.5            numcodecs==0.16.5
#> 21             numpy       2.4.4                 numpy==2.4.4
#> 22         packaging        26.2              packaging==26.2
#> 23            pandas       2.3.3                pandas==2.3.3
#> 24             patsy       1.0.2                 patsy==1.0.2
#> 25            pillow      12.2.0               pillow==12.2.0
#> 26       pynndescent       0.6.0           pynndescent==0.6.0
#> 27         pyparsing       3.3.2             pyparsing==3.3.2
#> 28   python-dateutil 2.9.0.post0 python-dateutil==2.9.0.post0
#> 29              pytz      2026.2                 pytz==2026.2
#> 30            PyYAML       6.0.3                PyYAML==6.0.3
#> 31            scanpy      1.12.1               scanpy==1.12.1
#> 32      scikit-learn       1.8.0          scikit-learn==1.8.0
#> 33             scipy      1.17.1                scipy==1.17.1
#> 34      scverse-misc       0.0.5          scverse-misc==0.0.5
#> 35           seaborn      0.13.2              seaborn==0.13.2
#> 36     session-info2       0.4.1         session-info2==0.4.1
#> 37               six      1.17.0                  six==1.17.0
#> 38       statsmodels      0.14.6          statsmodels==0.14.6
#> 39     threadpoolctl       3.6.0         threadpoolctl==3.6.0
#> 40              tqdm      4.67.3                 tqdm==4.67.3
#> 41 typing_extensions      4.15.0    typing_extensions==4.15.0
#> 42            tzdata      2026.2               tzdata==2026.2
#> 43        umap-learn      0.5.12           umap-learn==0.5.12
#> 44              zarr       3.2.1                  zarr==3.2.1
```
