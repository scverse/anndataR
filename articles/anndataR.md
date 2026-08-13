# Using anndataR

## Introduction

*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* allows
users to work with `.h5ad` files and `.zarr` stores, interact with
`AnnData` objects and convert to/from `SingleCellExperiment` or `Seurat`
objects. This enables users to move data easily between the different
programming languages and analysis ecosystems needed to perform
single-cell data analysis.

This package builds on our experience developing and using other
interoperability packages and aims to provide a first-class R `AnnData`
experience.

### Relationship to other packages

Existing packages provide similar functionality to
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* but there
are some important differences:

- *[zellkonverter](https://bioconductor.org/packages/3.24/zellkonverter)*
  provides conversion of `SingleCellExperiment` objects to/from
  `AnnData` and reading/writing of `.h5ad` files. This is facilitated
  via *[reticulate](https://CRAN.R-project.org/package=reticulate)*
  using *[basilisk](https://bioconductor.org/packages/3.24/basilisk)* to
  manage Python environments (native reading of `.h5ad` is also
  possible). In contrast,
  *[anndataR](https://bioconductor.org/packages/3.24/anndataR)* provides
  a native R H5AD and Zarr interface, removing the need for Python
  dependencies. Conversion to/from `Seurat` objects is also supported.
- *[anndata](https://CRAN.R-project.org/package=anndata)* (on CRAN) is a
  wrapper around the Python *anndata* package. It provides a nicer
  interface from within R but still requires a Python environment.
  *[anndataR](https://bioconductor.org/packages/3.24/anndataR)* provides
  a pure R implementation of the `AnnData` data structure with different
  back ends and conversion to more common R objects.
- Other interoperability packages typically also require Python
  dependencies, have limited features or flexibility and/or are not
  available from a central package repository

There is significant overlap in functionality between these packages. We
anticipate that over time they will become more specialised and work
better together or be deprecated in favour of
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)*.

## Installation

Install *[anndataR](https://bioconductor.org/packages/3.24/anndataR)*
using *[BiocManager](https://CRAN.R-project.org/package=BiocManager)*:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("anndataR")
```

## Usage

These sections briefly show how to use
*[anndataR](https://bioconductor.org/packages/3.24/anndataR)*.

### Read from disk

First, we fetch an example `.h5ad` file included in the package:

``` r

library(anndataR)

h5ad_path <- system.file("extdata", "example.h5ad", package = "anndataR")
```

By default, a H5AD is read to an in-memory `AnnData` object:

``` r

adata <- read_h5ad(h5ad_path)
```

It can also be read as a `SingleCellExperiment` object (see
[`vignette("usage_singlecellexperiment")`](https://anndataR.scverse.org/articles/usage_singlecellexperiment.md)):

``` r

sce <- read_h5ad(h5ad_path, as = "SingleCellExperiment")
```

Or as a `Seurat` object (see
[`vignette("usage_seurat")`](https://anndataR.scverse.org/articles/usage_seurat.md)):

``` r

obj <- read_h5ad(h5ad_path, as = "Seurat")
```

There is also a HDF5-backed `AnnData` object:

``` r

adata <- read_h5ad(h5ad_path, as = "HDF5AnnData")
```

Similarly, these functionalities are provided for `.zarr` stores too.

``` r

zarr_path <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(zarr_path, exdir = td)
zarr_path <- file.path(td, "example_v2.zarr")
```

``` r

# in-memory
adata <- read_zarr(zarr_path)

# as SingleCellExperiment
sce <- read_zarr(zarr_path, as = "SingleCellExperiment")

# as Seurat
obj <- read_zarr(zarr_path, as = "Seurat")

# as Zarr-backed
adata <- read_zarr(zarr_path, as = "ZarrAnnData")
```

See for interacting with a Python `AnnData` via
*[reticulate](https://CRAN.R-project.org/package=reticulate)*.

### Using `AnnData` objects

View structure:

``` r

adata
#> ZarrAnnData object with n_obs × n_vars = 50 × 100
#>     obs: 'Float', 'FloatNA', 'Int', 'IntNA', 'Bool', 'BoolNA', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'leiden'
#>     var: 'String', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#>     uns: 'Bool', 'BoolNA', 'Category', 'DataFrameEmpty', 'hvg', 'Int', 'IntNA', 'IntScalar', 'leiden', 'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'Sparse1D', 'String', 'String2D', 'StringScalar', 'umap'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'PCs'
#>     layers: 'counts', 'csc_counts', 'dense_counts', 'dense_X'
#>     obsp: 'connectivities', 'distances'
#>     varp: 'test_varp'
```

Access `AnnData` slots:

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
#>          String n_cells_by_counts mean_counts log1p_mean_counts
#> Gene000 String0                44        1.94          1.078410
#> Gene001 String1                42        2.04          1.111858
#> Gene002 String2                43        2.12          1.137833
#> Gene003 String3                41        1.72          1.000632
#> Gene004 String4                42        2.06          1.118415
#>         pct_dropout_by_counts total_counts
#> Gene000                    12           97
#> Gene001                    16          102
#> Gene002                    14          106
#> Gene003                    18           86
#> Gene004                    16          103
```

[Subsetting](#subsetting) `AnnData` objects is covered below. See
[`` ?`AnnData-usage` ``](https://anndataR.scverse.org/reference/AnnData-usage.md)for
more details on how to work with `AnnData` objects.

### Interoperability

Convert the `AnnData` object to a `SingleCellExperiment` object (see
[`vignette("usage_singlecellexperiment")`](https://anndataR.scverse.org/articles/usage_singlecellexperiment.md)):

``` r

sce <- adata$as_SingleCellExperiment()
sce
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

Convert the `AnnData` object to a `Seurat` object (see
[`vignette("usage_seurat")`](https://anndataR.scverse.org/articles/usage_seurat.md)):

``` r

obj <- adata$as_Seurat()
obj
#> An object of class Seurat 
#> 100 features across 50 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  5 layers present: counts, csc_counts, dense_counts, dense_X, X
#>  2 dimensional reductions calculated: X_pca, X_umap
```

Convert a `SingleCellExperiment` object to an `AnnData` object (see
[`vignette("usage_singlecellexperiment")`](https://anndataR.scverse.org/articles/usage_singlecellexperiment.md)):

``` r

adata <- as_AnnData(sce)
adata
#> InMemoryAnnData object with n_obs × n_vars = 50 × 100
#>     obs: 'Float', 'FloatNA', 'Int', 'IntNA', 'Bool', 'BoolNA', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'leiden'
#>     var: 'String', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#>     uns: 'Bool', 'BoolNA', 'Category', 'DataFrameEmpty', 'hvg', 'Int', 'IntNA', 'IntScalar', 'leiden', 'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'Sparse1D', 'String', 'String2D', 'StringScalar', 'umap'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'X_pca'
#>     layers: 'counts', 'csc_counts', 'dense_counts', 'dense_X', 'X'
#>     obsp: 'connectivities', 'distances'
#>     varp: 'test_varp'
```

Convert a `Seurat` object to an `AnnData` object (see
[`vignette("usage_seurat")`](https://anndataR.scverse.org/articles/usage_seurat.md)):

``` r

adata <- as_AnnData(obj)
adata
#> InMemoryAnnData object with n_obs × n_vars = 50 × 100
#>     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Float', 'FloatNA', 'Int', 'IntNA', 'Bool', 'BoolNA', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'leiden'
#>     var: 'String', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#>     uns: 'Bool', 'BoolNA', 'Category', 'DataFrameEmpty', 'hvg', 'Int', 'IntNA', 'IntScalar', 'leiden', 'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'Sparse1D', 'String', 'String2D', 'StringScalar', 'umap'
#>     obsm: 'X_pca', 'X_umap'
#>     layers: 'counts', 'csc_counts', 'dense_counts', 'dense_X', 'X'
#>     obsp: 'connectivities', 'distances'
```

### Manually create an `AnnData` object

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
#> InMemoryAnnData object with n_obs × n_vars = 10 × 10
#>     obs: 'cell_type'
#>     var: 'gene_name'
```

### Write to disk

Write an `AnnData` object to disk:

``` r

tmpfile <- tempfile(fileext = ".h5ad")
adata$write_h5ad(tmpfile) # Alternatively, write_h5ad(adata, tmpfile)
```

Write a `SingleCellExperiment` object to disk:

``` r

tmpfile <- tempfile(fileext = ".h5ad")
write_h5ad(sce, tmpfile)
#> Warning: Could not write element 'obsp/connectivities' of type 'dgTMatrix': Unsupported
#> matrix format in obsp/connectivities ℹ Supported matrices inherit from
#> <RsparseMatrix> or <CsparseMatrix>
#> Warning: Could not write element 'obsp/distances' of type 'dgTMatrix': Unsupported
#> matrix format in obsp/distances ℹ Supported matrices inherit from
#> <RsparseMatrix> or <CsparseMatrix>
#> Warning: Could not write element 'varp/test_varp' of type 'dgTMatrix': Unsupported
#> matrix format in varp/test_varp ℹ Supported matrices inherit from
#> <RsparseMatrix> or <CsparseMatrix>
```

Write a `Seurat` object to disk:

``` r

tmpfile <- tempfile(fileext = ".h5ad")
write_h5ad(obj, tmpfile)
#> Warning: Matrix column names cannot be written to a <HDF5AnnData> object, they will be
#> lost
#> ℹ To write column names for obsm[['X_pca']], store it as <data.frame> instead
#>   of a double matrix
#> ℹ NOTE: obs_names and var_names are stored separately
#> Warning: Matrix column names cannot be written to a <HDF5AnnData> object, they will be
#> lost
#> ℹ To write column names for obsm[['X_umap']], store it as <data.frame> instead
#>   of a double matrix
#> ℹ NOTE: obs_names and var_names are stored separately
```

Similarly, we can write `AnnData` and other objects to `.zarr` stores
too.

``` r

tmpfile <- tempfile(fileext = ".zarr")
adata$write_zarr(tmpfile) # Alternatively, write_zarr(adata, tmpfile)

tmpfile <- tempfile(fileext = ".zarr")
write_zarr(sce, tmpfile)
#> Warning: Could not write element 'obsp/connectivities' of type 'dgTMatrix': Unsupported
#> matrix format in obsp/connectivities ℹ Supported matrices inherit from
#> <RsparseMatrix> or <CsparseMatrix>
#> Warning: Could not write element 'obsp/distances' of type 'dgTMatrix': Unsupported
#> matrix format in obsp/distances ℹ Supported matrices inherit from
#> <RsparseMatrix> or <CsparseMatrix>
#> Warning: Could not write element 'varp/test_varp' of type 'dgTMatrix': Unsupported
#> matrix format in varp/test_varp ℹ Supported matrices inherit from
#> <RsparseMatrix> or <CsparseMatrix>

tmpfile <- tempfile(fileext = ".zarr")
write_zarr(obj, tmpfile)
#> Warning: Matrix column names cannot be written to a <ZarrAnnData> object, they will be
#> lost
#> ℹ To write column names for obsm[['X_pca']], store it as <data.frame> instead
#>   of a double matrix
#> ℹ NOTE: obs_names and var_names are stored separately
#> Warning: Matrix column names cannot be written to a <ZarrAnnData> object, they will be
#> lost
#> ℹ To write column names for obsm[['X_umap']], store it as <data.frame> instead
#>   of a double matrix
#> ℹ NOTE: obs_names and var_names are stored separately
```

### Subsetting `AnnData` objects

*[anndataR](https://bioconductor.org/packages/3.24/anndataR)* provides
standard R subsetting methods that work with familiar bracket notation.
These methods return `AnnDataView` objects that provide lazy evaluation
for efficient memory usage.

#### Basic subsetting

Subset observations (rows) using logical conditions:

``` r

# Create some sample data
adata <- AnnData(
  X = matrix(rnorm(50), nrow = 10, ncol = 5),
  obs = data.frame(
    cell_type = factor(rep(c("A", "B", "C"), c(3, 4, 3))),
    score = runif(10)
  ),
  var = data.frame(
    gene_name = paste0("gene_", 1:5),
    highly_variable = c(TRUE, FALSE, TRUE, TRUE, FALSE)
  )
)

# Subset to cell type "A"
view1 <- adata[adata$obs$cell_type == "A", ]
view1
#> View of InMemoryAnnData object with n_obs × n_vars = 3 × 5
#>     obs: 'cell_type', 'score'
#>     var: 'gene_name', 'highly_variable'
```

Subset variables (columns) using logical conditions:

``` r

# Subset to highly variable genes
view2 <- adata[, adata$var$highly_variable]
view2
#> View of InMemoryAnnData object with n_obs × n_vars = 10 × 3
#>     obs: 'cell_type', 'score'
#>     var: 'gene_name', 'highly_variable'
```

#### Combined subsetting

Subset both observations and variables simultaneously:

``` r

# Subset to cell type "A" and highly variable genes
view3 <- adata[adata$obs$cell_type == "A", adata$var$highly_variable]
view3
#> View of InMemoryAnnData object with n_obs × n_vars = 3 × 3
#>     obs: 'cell_type', 'score'
#>     var: 'gene_name', 'highly_variable'
```

#### Using different index types

``` r

# Numeric indices
view4 <- adata[1:5, 1:3]
view4
#> View of InMemoryAnnData object with n_obs × n_vars = 5 × 3
#>     obs: 'cell_type', 'score'
#>     var: 'gene_name', 'highly_variable'

# Character names (if available)
rownames(adata) <- paste0("cell_", 1:10)
colnames(adata) <- paste0("gene_", 1:5)
view5 <- adata[c("cell_1", "cell_2"), c("gene_1", "gene_3")]
view5
#> View of InMemoryAnnData object with n_obs × n_vars = 2 × 2
#>     obs: 'cell_type', 'score'
#>     var: 'gene_name', 'highly_variable'
```

#### Working with views

Views maintain access to all original data slots:

``` r

# Access dimensions
dim(view3)
#> [1] 3 3
nrow(view3)
#> [1] 3
ncol(view3)
#> [1] 3

# Access names
rownames(view3)
#> [1] "cell_1" "cell_2" "cell_3"
colnames(view3)
#> [1] "gene_1" "gene_3" "gene_4"

# Access data matrices and metadata
view3$X
#>            gene_1      gene_3     gene_4
#> cell_1 -0.6905379  0.02229473  0.4543416
#> cell_2 -0.5585420  0.60361101 -0.8552025
#> cell_3 -0.5366633 -0.26265057 -0.2868952
view3$obs
#>        cell_type     score
#> cell_1         A 0.3492990
#> cell_2         A 0.9473183
#> cell_3         A 0.2161000
view3$var
#>        gene_name highly_variable
#> gene_1    gene_1            TRUE
#> gene_3    gene_3            TRUE
#> gene_4    gene_4            TRUE

# Views can be converted to concrete implementations
concrete <- view3$as_InMemoryAnnData()
concrete
#> InMemoryAnnData object with n_obs × n_vars = 3 × 3
#>     obs: 'cell_type', 'score'
#>     var: 'gene_name', 'highly_variable'
```

## Citing *anndataR*

If you use *[anndataR](https://bioconductor.org/packages/3.24/anndataR)*
in your work, please cite [*“anndataR improves interoperability between
R and Python in single-cell
transcriptomics”*](https://doi.org/10.1101/2025.08.18.669052):

``` r

citation("anndataR")
#> To cite package 'anndataR' in publications use:
#> 
#>   Deconinck L, Zappia L, Cannoodt R, Morgan M, Virshup I, Sang-aram C,
#>   Bredikhin D, Seurinck R, Saeys Y, scverse core (2026). "anndataR
#>   improves interoperability between R and Python in single-cell
#>   transcriptomics." _Bioinformatics_, *42*, btag288.
#>   doi:10.1093/bioinformatics/btag288
#>   <https://doi.org/10.1093/bioinformatics/btag288>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {anndataR improves interoperability between R and Python in single-cell transcriptomics},
#>     author = {Louise Deconinck and Luke Zappia and Robrecht Cannoodt and Martin Morgan and Isaac Virshup and Chananchida Sang-aram and Danila Bredikhin and Ruth Seurinck and Yvan Saeys and {scverse core}},
#>     journal = {Bioinformatics},
#>     volume = {42},
#>     issue = {6},
#>     month = {June},
#>     year = {2026},
#>     pages = {btag288},
#>     doi = {10.1093/bioinformatics/btag288},
#>   }
```

## Session info

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] anndataR_1.2.1              SingleCellExperiment_1.35.2
#>  [3] SummarizedExperiment_1.43.0 Biobase_2.73.2             
#>  [5] GenomicRanges_1.65.1        Seqinfo_1.3.0              
#>  [7] IRanges_2.47.2              S4Vectors_0.51.6           
#>  [9] BiocGenerics_0.59.12        generics_0.1.4             
#> [11] MatrixGenerics_1.25.0       matrixStats_1.5.0          
#> [13] SeuratObject_5.4.0          sp_2.2-3                   
#> [15] BiocStyle_2.41.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.5        
#>   [4] spatstat.utils_3.2-4   farver_2.1.2           rmarkdown_2.31        
#>   [7] fs_2.1.0               ragg_1.5.2             vctrs_0.7.3           
#>  [10] ROCR_1.0-12            spatstat.explore_3.8-2 htmltools_0.5.9       
#>  [13] S4Arrays_1.13.0        curl_7.1.0             Rhdf5lib_2.1.0        
#>  [16] SparseArray_1.13.2     rhdf5_2.57.9           sass_0.4.10           
#>  [19] sctransform_0.4.3      parallelly_1.48.0      KernSmooth_2.23-26    
#>  [22] bslib_0.12.0           htmlwidgets_1.6.4      desc_1.4.3            
#>  [25] ica_1.0-3              httr2_1.3.0            plyr_1.8.9            
#>  [28] plotly_4.12.1          zoo_1.9-0              cachem_1.1.0          
#>  [31] igraph_2.3.3           mime_0.13              lifecycle_1.0.5       
#>  [34] pkgconfig_2.0.3        Matrix_1.7-5           R6_2.6.1              
#>  [37] fastmap_1.2.0          fitdistrplus_1.2-6     future_1.75.0         
#>  [40] shiny_1.14.0           digest_0.6.39          paws.storage_0.10.0   
#>  [43] patchwork_1.3.2        tensor_1.5.1           Seurat_5.5.1          
#>  [46] RSpectra_0.16-2        irlba_2.3.7            textshaping_1.0.5     
#>  [49] Rarr_2.1.32            progressr_1.0.0        spatstat.sparse_3.2-0 
#>  [52] polyclip_1.10-7        httr_1.4.8             abind_1.4-8           
#>  [55] compiler_4.6.1         S7_0.2.2               fastDummies_1.7.6     
#>  [58] grumpy_0.1.1           R.utils_2.13.0         MASS_7.3-65           
#>  [61] DelayedArray_0.39.5    tools_4.6.1            lmtest_0.9-40         
#>  [64] otel_0.2.0             httpuv_1.6.17          future.apply_1.20.2   
#>  [67] goftest_1.2-3          R.oo_1.27.1            glue_1.8.1            
#>  [70] nlme_3.1-169           rhdf5filters_1.25.4    promises_1.5.0        
#>  [73] grid_4.6.1             Rtsne_0.17             cluster_2.1.8.2       
#>  [76] reshape2_1.4.5         gtable_0.3.6           spatstat.data_3.1-9   
#>  [79] R.methodsS3_1.8.2      tidyr_1.3.2            data.table_1.18.4     
#>  [82] XVector_0.53.0         spatstat.geom_3.8-2    RcppAnnoy_0.0.23      
#>  [85] ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1         
#>  [88] stringr_1.6.0          spam_2.11-4            RcppHNSW_0.7.0        
#>  [91] later_1.4.8            splines_4.6.1          dplyr_1.2.1           
#>  [94] lattice_0.22-9         deldir_2.0-4           survival_3.8-6        
#>  [97] paws.common_0.8.10     tidyselect_1.2.1       miniUI_0.1.2          
#> [100] pbapply_1.7-4          knitr_1.51             gridExtra_2.3.1       
#> [103] bookdown_0.47          scattermore_1.2        xfun_0.60             
#> [106] stringi_1.8.9          yaml_2.3.12            evaluate_1.0.5        
#> [109] codetools_0.2-20       tibble_3.3.1           BiocManager_1.30.27   
#> [112] cli_3.6.6              uwot_0.2.4             xtable_1.8-8          
#> [115] reticulate_1.46.0      systemfonts_1.3.2      jquerylib_0.1.4       
#> [118] Rcpp_1.1.2             spatstat.random_3.5-1  globals_0.19.1        
#> [121] png_0.1-9              spatstat.univar_3.2-0  parallel_4.6.1        
#> [124] pkgdown_2.2.1          ggplot2_4.0.3          dotCall64_1.2         
#> [127] listenv_1.0.0          viridisLite_0.4.3      scales_1.4.0          
#> [130] ggridges_0.5.7         crayon_1.5.3           purrr_1.2.2           
#> [133] rlang_1.3.0            cowplot_1.2.0
```
