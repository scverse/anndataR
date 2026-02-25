# Convert to an `AnnData` object

Convert other objects to an `AnnData` object. See the sections below for
details on how slots are mapped between objects. For more information on
the functionality of an `AnnData` object, see
[AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md).

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
  output_class = c("InMemory", "HDF5AnnData", "ReticulateAnnData"),
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
  output_class = c("InMemory", "HDF5AnnData", "ReticulateAnnData"),
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
  output_class = c("InMemory", "HDF5AnnData", "ReticulateAnnData"),
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

|                                   |                  |                                                                  |                                                                                                                                                                                |
|-----------------------------------|------------------|------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **From `SingleCellExperiment`**   | **To `AnnData`** | **Example mapping argument**                                     | **Default**                                                                                                                                                                    |
| `assays(x)`                       | `adata$X`        | `x_mapping = "counts"`                                           | Nothing is copied to `X`                                                                                                                                                       |
| `assays(x)`                       | `adata$layers`   | `layers_mapping = c(counts = "counts")`                          | All items are copied by name                                                                                                                                                   |
| `colData(x)`                      | `adata$obs`      | `obs_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name                                                                                                                                                 |
| `rowData(x)`                      | `adata$var`      | `var_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")`     | All columns are copied by name                                                                                                                                                 |
| `reducedDims(x)`                  | `adata$obsm`     | `obsm_mapping = c(X_pca = "pca")`                                | All items are copied by name                                                                                                                                                   |
| `featureLoadings(reducedDims(x))` | `adata$varm`     | `varm_mapping = c(PCs = "pca")`                                  | Feature loadings from all [`SingleCellExperiment::LinearEmbeddingMatrix`](https://rdrr.io/pkg/SingleCellExperiment/man/LinearEmbeddingMatrix.html) objects in `reducedDims(x)` |
| `colPairs(x)`                     | `adata$obsp`     | `obsp_mapping = c(connectivities = "RNA_nn")`                    | All items are copied by name                                                                                                                                                   |
| `rowPairs(x)`                     | `adata$varp`     | `varp_mapping = c(similarities = "gene_overlaps")`               | All items are copied by name                                                                                                                                                   |
| `metadata(x)`                     | `adata$uns`      | `uns_mapping = c(metadata = "project_metadata")`                 | All items are copied by name                                                                                                                                                   |

## Converting from a `Seurat` object

Only one assay can be converted from a
[`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object to an `AnnData` object at a time. This can be controlled using
the `assay_name` argument. By default, the current default assay will be
used.

This table describes how slots in a
[`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object to the new `AnnData` object.

|                       |                  |                                                                  |                                                         |
|-----------------------|------------------|------------------------------------------------------------------|---------------------------------------------------------|
| **From `Seurat`**     | **To `AnnData`** | **Example mapping argument**                                     | **Default**                                             |
| `Layers(x)`           | `adata$X`        | `x_mapping = "counts"`                                           | Nothing is copied to `X`                                |
| `Layers(x)`           | `adata$layers`   | `layers_mapping = c(counts = "counts")`                          | All items are copied by name                            |
| `x[[]]`               | `adata$obs`      | `obs_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name                          |
| `x[[assay_name]][[]]` | `adata$var`      | `var_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")`     | All columns are copied by name                          |
| `Reductions(x)`       | `adata$obsm`     | `obsm_mapping = c(X_pca = "pca")`                                | All embeddings matching `assay_name` are copied by name |
| `Loadings(x)`         | `adata$varm`     | `varm_mapping = c(PCs = "pca")`                                  | All valid loadings are copied by name                   |
| `Graphs(x)`           | `adata$obsp`     | `obsp_mapping = c(connectivities = "RNA_nn")`                    | All graphs matching `assay_name` are copied by name     |
| `Misc(x)`             | `adata$varp`     | `varp_mapping = c(similarities = "gene_overlaps")`               | No data is copied to `varp`                             |
| `Misc(x)`             | `adata$uns`      | `uns_mapping = c(metadata = "project_metadata")`                 | All items are copied by name                            |

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
[`AnnData()`](https://anndataR.data-intuitive.com/reference/AnnData.md),
[`read_h5ad()`](https://anndataR.data-intuitive.com/reference/read_h5ad.md)

Other object converters:
[`as_HDF5AnnData()`](https://anndataR.data-intuitive.com/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.data-intuitive.com/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.data-intuitive.com/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.data-intuitive.com/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.data-intuitive.com/reference/as_SingleCellExperiment.md),
[`reticulate-helpers`](https://anndataR.data-intuitive.com/reference/reticulate-helpers.md)

## Examples

``` r
# Convert a Seurat object to an AnnData object
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> ‘SeuratObject’ was built under R 4.5.0 but the current version is
#> 4.5.2; it is recomended that you reinstall ‘SeuratObject’ as the ABI
#> for R may have changed
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
#> Positive:  Feature37, Feature76, Feature85, Feature72, Feature97, Feature14, Feature83, Feature69, Feature13, Feature92 
#>     Feature84, Feature95, Feature20, Feature5, Feature10, Feature33, Feature88, Feature32, Feature8, Feature60 
#>     Feature29, Feature93, Feature67, Feature1, Feature55, Feature61, Feature77, Feature74, Feature23, Feature28 
#> Negative:  Feature78, Feature9, Feature41, Feature4, Feature47, Feature66, Feature17, Feature51, Feature91, Feature46 
#>     Feature94, Feature21, Feature62, Feature90, Feature56, Feature49, Feature81, Feature19, Feature75, Feature50 
#>     Feature86, Feature63, Feature73, Feature48, Feature16, Feature53, Feature18, Feature15, Feature89, Feature39 
#> PC_ 2 
#> Positive:  Feature97, Feature67, Feature30, Feature5, Feature70, Feature23, Feature35, Feature22, Feature19, Feature2 
#>     Feature48, Feature77, Feature9, Feature41, Feature88, Feature34, Feature7, Feature21, Feature47, Feature58 
#>     Feature60, Feature78, Feature80, Feature4, Feature14, Feature62, Feature86, Feature96, Feature92, Feature12 
#> Negative:  Feature15, Feature59, Feature100, Feature79, Feature94, Feature98, Feature84, Feature37, Feature38, Feature18 
#>     Feature3, Feature65, Feature31, Feature64, Feature61, Feature87, Feature16, Feature24, Feature29, Feature99 
#>     Feature13, Feature72, Feature6, Feature93, Feature63, Feature36, Feature27, Feature33, Feature46, Feature75 
#> PC_ 3 
#> Positive:  Feature62, Feature77, Feature43, Feature56, Feature94, Feature17, Feature69, Feature52, Feature48, Feature18 
#>     Feature29, Feature85, Feature8, Feature32, Feature41, Feature11, Feature64, Feature63, Feature72, Feature42 
#>     Feature2, Feature99, Feature35, Feature86, Feature21, Feature31, Feature89, Feature16, Feature75, Feature20 
#> Negative:  Feature73, Feature30, Feature87, Feature49, Feature84, Feature78, Feature10, Feature40, Feature3, Feature34 
#>     Feature22, Feature70, Feature66, Feature53, Feature59, Feature65, Feature100, Feature82, Feature24, Feature97 
#>     Feature33, Feature28, Feature25, Feature93, Feature57, Feature81, Feature19, Feature91, Feature51, Feature4 
#> PC_ 4 
#> Positive:  Feature13, Feature55, Feature25, Feature4, Feature92, Feature39, Feature11, Feature98, Feature36, Feature94 
#>     Feature27, Feature2, Feature71, Feature24, Feature46, Feature67, Feature60, Feature50, Feature1, Feature28 
#>     Feature78, Feature41, Feature81, Feature80, Feature75, Feature10, Feature44, Feature58, Feature45, Feature14 
#> Negative:  Feature8, Feature3, Feature20, Feature5, Feature26, Feature43, Feature95, Feature63, Feature79, Feature87 
#>     Feature17, Feature70, Feature97, Feature90, Feature89, Feature18, Feature22, Feature51, Feature47, Feature61 
#>     Feature19, Feature52, Feature82, Feature76, Feature15, Feature38, Feature34, Feature12, Feature40, Feature84 
#> PC_ 5 
#> Positive:  Feature47, Feature68, Feature29, Feature35, Feature20, Feature95, Feature86, Feature66, Feature83, Feature11 
#>     Feature100, Feature70, Feature94, Feature1, Feature24, Feature10, Feature18, Feature62, Feature54, Feature39 
#>     Feature55, Feature4, Feature23, Feature16, Feature57, Feature48, Feature88, Feature15, Feature60, Feature82 
#> Negative:  Feature96, Feature32, Feature43, Feature22, Feature52, Feature59, Feature87, Feature71, Feature12, Feature56 
#>     Feature36, Feature2, Feature17, Feature72, Feature25, Feature41, Feature49, Feature92, Feature90, Feature34 
#>     Feature26, Feature37, Feature30, Feature69, Feature50, Feature99, Feature45, Feature79, Feature91, Feature13 
obj <- FindNeighbors(obj)
#> Computing nearest neighbor graph
#> Computing SNN
obj <- RunUMAP(obj, dims = 1:10)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
#> 05:21:58 UMAP embedding parameters a = 0.9922 b = 1.112
#> 05:21:58 Read 200 rows and found 10 numeric columns
#> 05:21:58 Using Annoy for neighbor search, n_neighbors = 30
#> 05:21:58 Building Annoy index with metric = cosine, n_trees = 50
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
#> 05:21:58 Writing NN index file to temp file /tmp/RtmpkWejZX/file1f1f598be27f
#> 05:21:58 Searching Annoy index using 1 thread, search_k = 3000
#> 05:21:58 Annoy recall = 100%
#> 05:21:58 Commencing smooth kNN distance calibration using 1 thread
#>  with target n_neighbors = 30
#> 05:21:59 Initializing from normalized Laplacian + noise (using RSpectra)
#> 05:21:59 Commencing optimization for 500 epochs, with 6160 positive edges
#> 05:21:59 Using rng type: pcg
#> 05:22:00 Optimization finished

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
