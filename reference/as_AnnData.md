# Convert to an `AnnData` object

Convert other objects to an `AnnData` object. See the sections below for
details on how slots are mapped between objects. For more information on
the functionality of an `AnnData` object, see
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md).

## Usage

``` r
as_AnnData(
  x,
  x_mapping = NULL,
  layers_mapping = TRUE,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  obsm_mapping = TRUE,
  varm_mapping = TRUE,
  obsp_mapping = TRUE,
  varp_mapping = TRUE,
  uns_mapping = TRUE,
  assay_name = NULL,
  output_class = c("InMemory", "HDF5AnnData", "ZarrAnnData", "ReticulateAnnData"),
  ...
)

# S3 method for class 'SingleCellExperiment'
as_AnnData(
  x,
  x_mapping = NULL,
  layers_mapping = TRUE,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  obsm_mapping = TRUE,
  varm_mapping = TRUE,
  obsp_mapping = TRUE,
  varp_mapping = TRUE,
  uns_mapping = TRUE,
  assay_name = TRUE,
  output_class = c("InMemory", "HDF5AnnData", "ZarrAnnData", "ReticulateAnnData"),
  ...
)

# S3 method for class 'Seurat'
as_AnnData(
  x,
  x_mapping = NULL,
  layers_mapping = TRUE,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  obsm_mapping = TRUE,
  varm_mapping = TRUE,
  obsp_mapping = TRUE,
  varp_mapping = TRUE,
  uns_mapping = TRUE,
  assay_name = NULL,
  output_class = c("InMemory", "HDF5AnnData", "ZarrAnnData", "ReticulateAnnData"),
  ...
)
```

## Arguments

- x:

  The object to convert

- x_mapping:

  A string specifying the data to map to the `X` slot. If `NULL`, no
  data will be copied to the `X` slot. See below for details.

- layers_mapping:

  A named character vector where the names are keys of `layers` in the
  new `AnnData` object and values are the names of items in the
  corresponding slot of `x`. See below for details.

- obs_mapping:

  A named character vector where the names are names of `obs` columns in
  the new `AnnData` object and values are the names of columns in the
  corresponding slot of `x`. See below for details.

- var_mapping:

  A named character vector where the names are names of `var` columns in
  the new `AnnData` object and values are the names of columns in the
  corresponding slot of `x`. See below for details.

- obsm_mapping:

  A named character vector where the names are keys of `obsm` in the new
  `AnnData` object and values are the names of items in the
  corresponding slot of `x`. See below for details.

- varm_mapping:

  A named character vector where the names are keys of `varm` in the new
  `AnnData` object and values are the names of items in the
  corresponding slot of `x`. See below for details.

- obsp_mapping:

  A named character vector where the names are keys of `obsp` in the new
  `AnnData` object and values are the names of items in the
  corresponding slot of `x`. See below for details.

- varp_mapping:

  A named character vector where the names are keys of `varp` in the new
  `AnnData` object and values are the names of items in the
  corresponding slot of `x`. See below for details.

- uns_mapping:

  A named character vector where the names are keys of `uns` in the new
  `AnnData` object and values are the names of items in the
  corresponding slot of `x`. See below for details.

- assay_name:

  For
  [`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  objects, the name of the assay to be converted. If `NULL`, the default
  assay will be used
  ([`SeuratObject::DefaultAssay()`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html)).
  This is ignored for other objects.

- output_class:

  The `AnnData` class to convert to. Must be one of `"HDF5AnnData"` or
  `"InMemoryAnnData"`.

- ...:

  Additional arguments passed to the generator function for
  `output_class`

## Value

An `AnnData` object of the class requested by `output_class` containing
the data specified in the mapping arguments.

## Details of mapping arguments

All mapping arguments except for `x_mapping` expect a named character
vector where names are the keys of the slot in the `AnnData` object and
values are the names of items in the corresponding slot of `x`. If
`TRUE`, the conversion function will guess which items to copy as
described in the conversion tables for each object type. In most cases,
the default is to copy all items using the same names except where the
correspondence between objects is unclear. To avoid copying anything to
a slot, set the mapping argument to `FALSE`. Empty mapping arguments
(`NULL`, [`c()`](https://rdrr.io/r/base/c.html),
[`list()`](https://rdrr.io/r/base/list.html)) will be treated as `FALSE`
with a warning. If an unnamed vector is provided, the values will be
used as names.

- `TRUE` will guess which items to copy as described in the conversion
  tables for each object type

- `c(adata_item = "x_item")` will copy `x_item` from the slot in `x` to
  `adata_item` in the corresponding slot of new `AnnData` object

- `FALSE` will avoid copying anything to the slot

- `c("x_item")` is equivalent to `c(x_item = "x_item")`

## Converting from a `SingleCellExperiment` object

This table describes how slots in a
[`SingleCellExperiment::SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object to the new `AnnData` object.

|  |  |  |  |
|----|----|----|----|
| **From `SingleCellExperiment`** | **To `AnnData`** | **Example mapping argument** | **Default** |
| `assays(x)` | `adata$X` | `x_mapping = "counts"` | Nothing is copied to `X` |
| `assays(x)` | `adata$layers` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
| `colData(x)` | `adata$obs` | `obs_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
| `rowData(x)` | `adata$var` | `var_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
| `reducedDims(x)` | `adata$obsm` | `obsm_mapping = c(X_pca = "pca")` | All items are copied by name |
| `featureLoadings(reducedDims(x))` | `adata$varm` | `varm_mapping = c(PCs = "pca")` | Feature loadings from all [`SingleCellExperiment::LinearEmbeddingMatrix`](https://rdrr.io/pkg/SingleCellExperiment/man/LinearEmbeddingMatrix.html) objects in `reducedDims(x)` |
| `colPairs(x)` | `adata$obsp` | `obsp_mapping = c(connectivities = "RNA_nn")` | All items are copied by name |
| `rowPairs(x)` | `adata$varp` | `varp_mapping = c(similarities = "gene_overlaps")` | All items are copied by name |
| `metadata(x)` | `adata$uns` | `uns_mapping = c(metadata = "project_metadata")` | All items are copied by name |

### Unnamed assays

If `assayNames(x)` is `NULL` or any assay names are empty they will
automatically be named with a warning:

**Examples:**

- Old names: `NULL` -\> New names: `"assay1", "assay2", ...`

- Old names: `"counts"` -\> New names: `"counts", "assay2"`

## Converting from a `Seurat` object

Only one assay can be converted from a
[`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object to an `AnnData` object at a time. This can be controlled using
the `assay_name` argument. By default, the current default assay will be
used.

This table describes how slots in a
[`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object to the new `AnnData` object.

|  |  |  |  |
|----|----|----|----|
| **From `Seurat`** | **To `AnnData`** | **Example mapping argument** | **Default** |
| `Layers(x)` | `adata$X` | `x_mapping = "counts"` | Nothing is copied to `X` |
| `Layers(x)` | `adata$layers` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
| `x[[]]` | `adata$obs` | `obs_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
| `x[[assay_name]][[]]` | `adata$var` | `var_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
| `Reductions(x)` | `adata$obsm` | `obsm_mapping = c(X_pca = "pca")` | All embeddings matching `assay_name` are copied by name |
| `Loadings(x)` | `adata$varm` | `varm_mapping = c(PCs = "pca")` | All valid loadings are copied by name |
| `Graphs(x)` | `adata$obsp` | `obsp_mapping = c(connectivities = "RNA_nn")` | All graphs matching `assay_name` are copied by name |
| `Misc(x)` | `adata$varp` | `varp_mapping = c(similarities = "gene_overlaps")` | No data is copied to `varp` |
| `Misc(x)` | `adata$uns` | `uns_mapping = c(metadata = "project_metadata")` | All items are copied by name |

### Graph conversion

By default, all graphs in a
[`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object that match the assay being converted are copied to the `obsp`
slot of the new `AnnData` object. If a graph does not have an associated
assay:

- If `assay_name` is the default assay, they will be *converted* with a
  warning

- if `assay_name` is not the default assay, they will be *skipped* with
  a warning

To override this behaviour, provide a custom mapping using the
`obsp_mapping` argument.

### Unexpected dimensions

A
[`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
is more flexible in terms of the dimensions of items that can be stored
in various slots. For example, a `Layer` does not have to match the
dimensions of the whole object. If an item has unexpected dimensions, it
will be skipped with a warning.

## See also

Other AnnData creators:
[`AnnData()`](https://anndataR.scverse.org/reference/AnnData.md),
[`read_h5ad()`](https://anndataR.scverse.org/reference/read_h5ad.md),
[`read_zarr()`](https://anndataR.scverse.org/reference/read_zarr.md)

Other object converters:
[`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md),
[`as_ZarrAnnData()`](https://anndataR.scverse.org/reference/as_ZarrAnnData.md),
[`reticulate-helpers`](https://anndataR.scverse.org/reference/reticulate-helpers.md)

## Examples

``` r
# Convert a Seurat object to an AnnData object
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: ‘SeuratObject’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t

counts <- matrix(rbinom(20000, 1000, .001), nrow = 100)
obj <- CreateSeuratObject(counts = counts)
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
obj <- NormalizeData(obj)
#> Normalizing layer: counts
obj <- FindVariableFeatures(obj)
#> Finding variable features for layer counts
obj <- ScaleData(obj)
#> Centering and scaling data matrix
obj <- RunPCA(obj, npcs = 10L)
#> PC_ 1 
#> Positive:  Feature14, Feature88, Feature19, Feature51, Feature57, Feature29, Feature40, Feature31, Feature76, Feature56 
#>     Feature45, Feature17, Feature49, Feature66, Feature91, Feature83, Feature100, Feature77, Feature59, Feature54 
#>     Feature33, Feature96, Feature72, Feature22, Feature35, Feature21, Feature38, Feature81, Feature60, Feature50 
#> Negative:  Feature6, Feature47, Feature8, Feature82, Feature79, Feature18, Feature71, Feature95, Feature1, Feature86 
#>     Feature13, Feature2, Feature43, Feature9, Feature42, Feature39, Feature41, Feature48, Feature94, Feature36 
#>     Feature53, Feature69, Feature28, Feature93, Feature25, Feature30, Feature97, Feature84, Feature23, Feature74 
#> PC_ 2 
#> Positive:  Feature98, Feature77, Feature12, Feature33, Feature24, Feature70, Feature21, Feature11, Feature30, Feature65 
#>     Feature93, Feature87, Feature23, Feature90, Feature86, Feature7, Feature95, Feature42, Feature38, Feature15 
#>     Feature40, Feature52, Feature32, Feature35, Feature47, Feature79, Feature68, Feature1, Feature82, Feature20 
#> Negative:  Feature8, Feature89, Feature97, Feature28, Feature76, Feature73, Feature69, Feature83, Feature61, Feature63 
#>     Feature25, Feature27, Feature100, Feature2, Feature88, Feature59, Feature50, Feature99, Feature56, Feature4 
#>     Feature6, Feature19, Feature26, Feature5, Feature48, Feature57, Feature64, Feature36, Feature13, Feature29 
#> PC_ 3 
#> Positive:  Feature72, Feature87, Feature66, Feature27, Feature53, Feature58, Feature6, Feature18, Feature62, Feature31 
#>     Feature39, Feature73, Feature15, Feature28, Feature57, Feature79, Feature45, Feature51, Feature5, Feature100 
#>     Feature52, Feature30, Feature96, Feature95, Feature17, Feature19, Feature33, Feature99, Feature42, Feature29 
#> Negative:  Feature59, Feature40, Feature20, Feature83, Feature97, Feature88, Feature69, Feature9, Feature94, Feature13 
#>     Feature75, Feature34, Feature23, Feature38, Feature68, Feature35, Feature44, Feature32, Feature91, Feature14 
#>     Feature92, Feature65, Feature63, Feature50, Feature43, Feature1, Feature4, Feature47, Feature98, Feature71 
#> PC_ 4 
#> Positive:  Feature51, Feature56, Feature46, Feature23, Feature62, Feature42, Feature37, Feature81, Feature35, Feature55 
#>     Feature66, Feature69, Feature60, Feature53, Feature82, Feature1, Feature27, Feature4, Feature61, Feature22 
#>     Feature21, Feature86, Feature59, Feature17, Feature89, Feature25, Feature73, Feature16, Feature19, Feature31 
#> Negative:  Feature80, Feature58, Feature3, Feature15, Feature57, Feature45, Feature64, Feature78, Feature36, Feature8 
#>     Feature99, Feature77, Feature83, Feature63, Feature30, Feature84, Feature2, Feature96, Feature5, Feature7 
#>     Feature90, Feature9, Feature76, Feature24, Feature40, Feature70, Feature13, Feature92, Feature72, Feature34 
#> PC_ 5 
#> Positive:  Feature92, Feature30, Feature22, Feature89, Feature13, Feature62, Feature15, Feature83, Feature94, Feature53 
#>     Feature86, Feature18, Feature47, Feature11, Feature32, Feature38, Feature25, Feature40, Feature78, Feature31 
#>     Feature75, Feature33, Feature97, Feature29, Feature17, Feature100, Feature73, Feature20, Feature57, Feature68 
#> Negative:  Feature90, Feature9, Feature2, Feature58, Feature79, Feature35, Feature81, Feature14, Feature84, Feature95 
#>     Feature85, Feature60, Feature5, Feature16, Feature65, Feature96, Feature63, Feature28, Feature77, Feature26 
#>     Feature23, Feature34, Feature69, Feature42, Feature82, Feature49, Feature76, Feature7, Feature48, Feature51 
obj <- FindNeighbors(obj)
#> Computing nearest neighbor graph
#> Computing SNN
obj <- RunUMAP(obj, dims = 1:10)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
#> 17:46:30 UMAP embedding parameters a = 0.9922 b = 1.112
#> 17:46:30 Read 200 rows and found 10 numeric columns
#> 17:46:30 Using Annoy for neighbor search, n_neighbors = 30
#> 17:46:30 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> |
#> 17:46:30 Writing NN index file to temp file /tmp/RtmpxjP80R/file1f067cbbea0a
#> 17:46:30 Searching Annoy index using 1 thread, search_k = 3000
#> 17:46:30 Annoy recall = 100%
#> 17:46:31 Commencing smooth kNN distance calibration using 1 thread
#>  with target n_neighbors = 30
#> 17:46:32 Initializing from normalized Laplacian + noise (using RSpectra)
#> 17:46:32 Commencing optimization for 500 epochs, with 6146 positive edges
#> 17:46:32 Using rng type: pcg
#> 17:46:33 Optimization finished

as_AnnData(obj)
#> Warning: Row names of `Loadings(seurat_obj, "pca")` do not match the expected var names
#> ! The matrix will be expanded to include all var names.
#> InMemoryAnnData object with n_obs × n_vars = 200 × 100
#>     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA'
#>     var: 'vf_vst_counts_mean', 'vf_vst_counts_variance', 'vf_vst_counts_variance.expected', 'vf_vst_counts_variance.standardized', 'vf_vst_counts_variable', 'vf_vst_counts_rank', 'var.features', 'var.features.rank'
#>     obsm: 'pca', 'umap'
#>     varm: 'pca'
#>     layers: 'counts', 'data', 'scale.data'
#>     obsp: 'nn', 'snn'
# Convert a SingleCellExperiment object to an AnnData object
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
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
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> 
#> Attaching package: ‘IRanges’
#> The following object is masked from ‘package:sp’:
#> 
#>     %over%
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: ‘SummarizedExperiment’
#> The following object is masked from ‘package:Seurat’:
#> 
#>     Assays
#> The following object is masked from ‘package:SeuratObject’:
#> 
#>     Assays

sce <- SingleCellExperiment(
  assays = list(counts = matrix(1:5, 5L, 3L)),
  colData = DataFrame(cell = 1:3, row.names = paste0("Cell", 1:3)),
  rowData = DataFrame(gene = 1:5, row.names = paste0("Gene", 1:5))
)

as_AnnData(sce)
#> InMemoryAnnData object with n_obs × n_vars = 3 × 5
#>     obs: 'cell'
#>     var: 'gene'
#>     layers: 'counts'
```
