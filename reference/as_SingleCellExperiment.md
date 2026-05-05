# Convert an `AnnData` to a `SingleCellExperiment`

Convert an `AnnData` object to a `SingleCellExperiment` object

## Usage

``` r
as_SingleCellExperiment(
  adata,
  x_mapping = NULL,
  assays_mapping = TRUE,
  colData_mapping = TRUE,
  rowData_mapping = TRUE,
  reducedDims_mapping = TRUE,
  colPairs_mapping = TRUE,
  rowPairs_mapping = TRUE,
  metadata_mapping = TRUE
)
```

## Arguments

- adata:

  The `AnnData` object to convert

- x_mapping:

  A string specifying the name of the assay in the resulting
  `SingleCellExperiment` where the data in the `X` slot of `adata` will
  be mapped to

- assays_mapping:

  A named vector where names are names of `assays` in the resulting
  `SingleCellExperiment` object and values are keys of `layers` in
  `adata`. See below for details.

- colData_mapping:

  A named vector where names are columns of `colData` in the resulting
  `SingleCellExperiment` object and values are columns of `obs` in
  `adata`. See below for details.

- rowData_mapping:

  A named vector where names are columns of `rowData` in the resulting
  `SingleCellExperiment` object and values are columns of `var` in
  `adata`. See below for details.

- reducedDims_mapping:

  A named vector where names are names of `reducedDims` in the resulting
  `SingleCellExperiment` object and values are keys of `obsm` in
  `adata`. Alternatively, a named list where names are names of
  `reducedDims` in the resulting `SingleCellExperiment` object and
  values are vectors with the items `"sampleFactors"` and
  `"featureLoadings"` and/or `"metadata"`. See below for details.

- colPairs_mapping:

  A named vector where names are names of `colPairs` in the resulting
  `SingleCellExperiment` object and values are keys of `obsp` in
  `adata`. See below for details.

- rowPairs_mapping:

  A named vector where names are names of `rowPairs` in the resulting
  `SingleCellExperiment` object and values are keys of `varp` in
  `adata`. See below for details.

- metadata_mapping:

  A named vector where names are names of `metadata` in the resulting
  `SingleCellExperiment` object and values are keys of `uns` in `adata`.
  See below for details.

## Value

A `SingleCellExperiment` object containing the requested data from
`adata`

## Details

### Mapping arguments

All mapping arguments expect a named character vector where names are
the names of the slot in the `SingleCellExperiment` object and values
are the keys of the corresponding slot of `adata`. If `TRUE`, the
conversion function will guess which items to copy as described in the
conversion tables below. In most cases, the default is to copy all items
using the same names except where the correspondence between objects is
unclear. The `reducedDims_mapping` argument can also accept a more
complex list format, see below for details. To avoid copying anything to
a slot, set the mapping argument to `FALSE`. Empty mapping arguments
(`NULL`, [`c()`](https://rdrr.io/r/base/c.html),
[`list()`](https://rdrr.io/r/base/list.html)) will be treated as `FALSE`
with a warning. If an unnamed vector is provided, the values will be
used as names.

#### Examples:

- `TRUE` will guess which items to copy as described in the conversion
  table

- `c(sce_item = "adata_item")` will copy `adata_item` from the slot in
  `adata` to `sce_item` in the corresponding slot of the new
  `SingleCellExperiment` object

- `FALSE` will avoid copying anything to the slot

- `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`

### Conversion table

|  |  |  |  |
|----|----|----|----|
| **From `AnnData`** | **To `SingleCellExperiment`** | **Example mapping argument** | **Default** |
| `adata$X` | `assays(sce)` | `x_mapping = "counts"` | The data in `adata$X` is copied to the assay named `X` |
| `adata$layers` | `assays(sce)` | `assays_mapping = c(counts = "counts")` | All items are copied by name |
| `adata$obs` | `colData(sce)` | `colData_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
| `adata$var` | `rowData(sce)` | `rowData_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
| `adata$obsm` | `reducedDims(sce)` | `reducedDims_mapping = c(pca = "X_pca")` **OR** `reducedDims_mapping = list(pca = c(sampleFactors = "X_pca", featureLoadings = "PCs", metadata = "pca_metadata"))` | All items are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
| `adata$obsp` | `colPairs(sce)` | `colPairs_mapping = c(nn = "connectivities")` | All items are copied by name |
| `adata$varp` | `rowPairs(sce)` | `rowPairs_mapping = c(gene_overlaps = "similarities")` | All items are copied by name |
| `adata$uns` | `metadata(sce)` | `uns_mapping = c(project_metadata = "metadata")` | All items are copied by name |

### The `reducedDims_mapping` argument

For the simpler named vector format, the names should be the names of
`reducedDims` in the resulting `SingleCellExperiment` object, and the
values should be the keys of `obsm` in `adata`.

For more advanced mapping, use the list format where each item is a
vector with the following names used to create a
[`SingleCellExperiment::LinearEmbeddingMatrix`](https://rdrr.io/pkg/SingleCellExperiment/man/LinearEmbeddingMatrix.html)
(if `featureLoadings` or `metadata` is provided):

- `sampleFactors`: a key of the `obsm` slot in `adata`,
  `adata$obsm[[sampleFactors]]` is passed to the `sampleFactors`
  argument

- `featureLoadings`: a key of the `varm` slot in `adata` (optional),
  `adata$varm[[featureLoadings]]` is passed to the `featureLoadings`
  argument

- `metadata`: a key of the `uns` slot in `adata` (optional),
  `adata$uns[[metadata]]` is passed to the `metadata` argument

### The `x_mapping` and `assays_mapping` arguments

In order to specify where the data in `adata$X` will be stored in the
`assays(sce)` slot of the resulting object, you can use either the
`x_mapping` argument or the `assays_mapping` argument. If you use
`x_mapping`, it should be a string specifying the name of the layer in
`assays(sce)` where the data in `adata$X` will be stored. If you use
`assays_mapping`, it should be a named vector where names are names of
`assays(sce)` and values are keys of `layers` in `adata`. In order to
indicate the `adata$X` slot, you use `NA` as the value in the vector.
The name you provide for `x_mapping` may not be a name in
`assays_mapping`.

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md),
[`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md),
[`as_ZarrAnnData()`](https://anndataR.scverse.org/reference/as_ZarrAnnData.md),
[`reticulate-helpers`](https://anndataR.scverse.org/reference/reticulate-helpers.md)

## Examples

``` r
  ad <- AnnData(
    X = matrix(1:5, 3L, 5L),
    layers = list(A = matrix(5:1, 3L, 5L), B = matrix(letters[1:5], 3L, 5L)),
    obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
    var = data.frame(row.names = letters[1:5], gene = 1:5)
  )

  # Default usage
  sce <- ad$as_SingleCellExperiment(
    assays_mapping = TRUE,
    colData_mapping = TRUE,
    rowData_mapping = TRUE,
    reducedDims_mapping = TRUE,
    colPairs_mapping = TRUE,
    rowPairs_mapping = TRUE,
    metadata_mapping = TRUE
 )
```
