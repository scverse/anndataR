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

# TEMP(anndata<0.13.0): anndataR does not yet support Python anndata >= 0.13.0
reticulate::py_require("scanpy>=1.10")
reticulate::py_require("anndata<0.13.0")
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
#> R version 4.6.1 (2026-06-24)
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
#> [1] reticulate_1.46.0 anndataR_1.2.1    BiocStyle_2.41.0 
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] magrittr_2.0.5              spatstat.utils_3.2-4       
#>   [5] farver_2.1.2                rmarkdown_2.31             
#>   [7] fs_2.1.0                    ragg_1.5.2                 
#>   [9] vctrs_0.7.3                 ROCR_1.0-12                
#>  [11] memoise_2.0.1               spatstat.explore_3.8-2     
#>  [13] htmltools_0.5.9             S4Arrays_1.13.0            
#>  [15] curl_7.1.0                  SparseArray_1.13.2         
#>  [17] sass_0.4.10                 sctransform_0.4.3          
#>  [19] parallelly_1.48.0           KernSmooth_2.23-26         
#>  [21] bslib_0.12.0                htmlwidgets_1.6.4          
#>  [23] desc_1.4.3                  ica_1.0-3                  
#>  [25] httr2_1.3.0                 plyr_1.8.9                 
#>  [27] plotly_4.12.1               zoo_1.9-0                  
#>  [29] cachem_1.1.0                igraph_2.3.3               
#>  [31] mime_0.13                   lifecycle_1.0.5            
#>  [33] pkgconfig_2.0.3             Matrix_1.7-5               
#>  [35] R6_2.6.1                    fastmap_1.2.0              
#>  [37] MatrixGenerics_1.25.0       fitdistrplus_1.2-6         
#>  [39] future_1.75.0               shiny_1.14.0               
#>  [41] digest_0.6.39               patchwork_1.3.2            
#>  [43] S4Vectors_0.51.6            Seurat_5.5.1               
#>  [45] tensor_1.5.1                RSpectra_0.16-2            
#>  [47] irlba_2.3.7                 RSQLite_3.53.3             
#>  [49] textshaping_1.0.5           GenomicRanges_1.65.1       
#>  [51] filelock_1.0.3              progressr_1.0.0            
#>  [53] spatstat.sparse_3.2-0       httr_1.4.8                 
#>  [55] polyclip_1.10-7             abind_1.4-8                
#>  [57] compiler_4.6.1              withr_3.0.3                
#>  [59] bit64_4.8.2                 S7_0.2.2                   
#>  [61] DBI_1.3.0                   fastDummies_1.7.6          
#>  [63] MASS_7.3-65                 DelayedArray_0.39.5        
#>  [65] tools_4.6.1                 lmtest_0.9-40              
#>  [67] otel_0.2.0                  httpuv_1.6.17              
#>  [69] future.apply_1.20.2         goftest_1.2-3              
#>  [71] glue_1.8.1                  nlme_3.1-169               
#>  [73] promises_1.5.0              grid_4.6.1                 
#>  [75] Rtsne_0.17                  cluster_2.1.8.2            
#>  [77] reshape2_1.4.5              generics_0.1.4             
#>  [79] gtable_0.3.6                spatstat.data_3.1-9        
#>  [81] tidyr_1.3.2                 data.table_1.18.4          
#>  [83] sp_2.2-3                    XVector_0.53.0             
#>  [85] BiocGenerics_0.59.12        spatstat.geom_3.8-2        
#>  [87] RcppAnnoy_0.0.23            ggrepel_0.9.8              
#>  [89] RANN_2.6.2                  pillar_1.11.1              
#>  [91] stringr_1.6.0               spam_2.11-4                
#>  [93] RcppHNSW_0.7.0              later_1.4.8                
#>  [95] splines_4.6.1               dplyr_1.2.1                
#>  [97] BiocFileCache_3.3.0         lattice_0.22-9             
#>  [99] bit_4.6.0                   deldir_2.0-4               
#> [101] survival_3.8-6              tidyselect_1.2.1           
#> [103] SingleCellExperiment_1.35.2 miniUI_0.1.2               
#> [105] pbapply_1.7-4               knitr_1.51                 
#> [107] gridExtra_2.3.1             bookdown_0.47              
#> [109] IRanges_2.47.2              Seqinfo_1.3.0              
#> [111] SummarizedExperiment_1.43.0 scattermore_1.2            
#> [113] stats4_4.6.1                xfun_0.60                  
#> [115] Biobase_2.73.2              matrixStats_1.5.0          
#> [117] stringi_1.8.9               yaml_2.3.12                
#> [119] evaluate_1.0.5              codetools_0.2-20           
#> [121] tibble_3.3.1                BiocManager_1.30.27        
#> [123] cli_3.6.6                   uwot_0.2.4                 
#> [125] xtable_1.8-8                systemfonts_1.3.2          
#> [127] jquerylib_0.1.4             Rcpp_1.1.2                 
#> [129] spatstat.random_3.5-1       globals_0.19.1             
#> [131] dbplyr_2.6.0                png_0.1-9                  
#> [133] spatstat.univar_3.2-0       parallel_4.6.1             
#> [135] blob_1.3.0                  pkgdown_2.2.1              
#> [137] ggplot2_4.0.3               dotCall64_1.2              
#> [139] listenv_1.0.0               viridisLite_0.4.3          
#> [141] scales_1.4.0                ggridges_0.5.7             
#> [143] SeuratObject_5.4.0          purrr_1.2.2                
#> [145] rlang_1.3.0                 cowplot_1.2.0
```

### Python

``` r

reticulate::py_config()
#> python:         /opt/hostedtoolcache/Python/3.14.6/x64/bin/python3
#> libpython:      /opt/hostedtoolcache/Python/3.14.6/x64/lib/libpython3.14.so
#> pythonhome:     /opt/hostedtoolcache/Python/3.14.6/x64:/opt/hostedtoolcache/Python/3.14.6/x64
#> version:        3.14.6 (main, Jun 10 2026, 14:29:35) [GCC 13.3.0]
#> numpy:          /opt/hostedtoolcache/Python/3.14.6/x64/lib/python3.14/site-packages/numpy
#> numpy_version:  2.5.2
#> scanpy:         /opt/hostedtoolcache/Python/3.14.6/x64/lib/python3.14/site-packages/scanpy
#> 
#> NOTE: Python version was forced by RETICULATE_PYTHON

reticulate::py_list_packages()
#>              package      version                  requirement
#> 1            anndata      0.12.19             anndata==0.12.19
#> 2   array-api-compat       1.15.0     array-api-compat==1.15.0
#> 3            certifi    2026.7.22           certifi==2026.7.22
#> 4          contourpy        1.3.3             contourpy==1.3.3
#> 5             cycler       0.12.1               cycler==0.12.1
#> 6             donfig  0.8.1.post1          donfig==0.8.1.post1
#> 7   fast-array-utils          1.5        fast-array-utils==1.5
#> 8          fonttools       4.63.0            fonttools==4.63.0
#> 9      google-crc32c        1.8.0         google-crc32c==1.8.0
#> 10              h5py       3.16.0                 h5py==3.16.0
#> 11            joblib        1.5.3                joblib==1.5.3
#> 12        kiwisolver        1.5.0            kiwisolver==1.5.0
#> 13   legacy-api-wrap          1.5         legacy-api-wrap==1.5
#> 14          llvmlite       0.49.0             llvmlite==0.49.0
#> 15        matplotlib       3.11.1           matplotlib==3.11.1
#> 16            mudata       0.3.10               mudata==0.3.10
#> 17          narwhals       2.24.0             narwhals==2.24.0
#> 18           natsort        8.4.0               natsort==8.4.0
#> 19          networkx        3.6.1              networkx==3.6.1
#> 20             numba       0.67.0                numba==0.67.0
#> 21         numcodecs       0.16.5            numcodecs==0.16.5
#> 22             numpy        2.5.2                 numpy==2.5.2
#> 23         packaging         26.3              packaging==26.3
#> 24            pandas        2.3.3                pandas==2.3.3
#> 25             patsy        1.0.2                 patsy==1.0.2
#> 26            pillow       12.3.0               pillow==12.3.0
#> 27       pynndescent        0.6.0           pynndescent==0.6.0
#> 28         pyparsing        3.3.2             pyparsing==3.3.2
#> 29   python-dateutil  2.9.0.post0 python-dateutil==2.9.0.post0
#> 30              pytz 2026.3.post1           pytz==2026.3.post1
#> 31            PyYAML        6.0.3                PyYAML==6.0.3
#> 32            scanpy       1.12.3               scanpy==1.12.3
#> 33      scikit-learn        1.9.0          scikit-learn==1.9.0
#> 34             scipy       1.18.0                scipy==1.18.0
#> 35      scverse-misc        0.1.1          scverse-misc==0.1.1
#> 36           seaborn       0.13.2              seaborn==0.13.2
#> 37     session-info2        0.4.2         session-info2==0.4.2
#> 38               six       1.17.0                  six==1.17.0
#> 39       statsmodels       0.14.6          statsmodels==0.14.6
#> 40     threadpoolctl        3.6.0         threadpoolctl==3.6.0
#> 41              tqdm       4.70.0                 tqdm==4.70.0
#> 42 typing_extensions       4.16.0    typing_extensions==4.16.0
#> 43            tzdata       2026.3               tzdata==2026.3
#> 44        umap-learn       0.5.12           umap-learn==0.5.12
#> 45              zarr        3.3.0                  zarr==3.3.0
```
