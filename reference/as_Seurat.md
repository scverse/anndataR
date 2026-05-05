# Convert an `AnnData` to a `Seurat`

Convert an `AnnData` object to a `Seurat` object

## Usage

``` r
as_Seurat(
  adata,
  assay_name = "RNA",
  x_mapping = NULL,
  layers_mapping = TRUE,
  object_metadata_mapping = TRUE,
  assay_metadata_mapping = TRUE,
  reduction_mapping = TRUE,
  graph_mapping = TRUE,
  misc_mapping = TRUE
)
```

## Arguments

- adata:

  The `AnnData` object to convert.

- assay_name:

  Name of the assay to be created in the new `Seurat` object

- x_mapping:

  A string specifying the name of the layer in the resulting `Seurat`
  object where the data in the `X` slot of `adata` will be mapped to

- layers_mapping:

  A named vector where names are names of `Layers` in the resulting
  `Seurat` object and values are keys of `layers` in `adata`. See below
  for details.

- object_metadata_mapping:

  A named vector where names are cell metadata columns in the resulting
  `Seurat` object and values are columns of `obs` in `adata`. See below
  for details.

- assay_metadata_mapping:

  A named vector where names are variable metadata columns in the assay
  of the resulting `Seurat` object and values are columns of `var` in
  `adata`. See below for details.

- reduction_mapping:

  A named vector where names are names of `Embeddings` in the resulting
  `Seurat` object and values are keys of `obsm` in `adata`.
  Alternatively, a named list where names are names of `Embeddings` in
  the resulting `Seurat` object and values are vectors with the items
  `"key"`, `"embeddings"` and (optionally) `"loadings"`. See below for
  details.

- graph_mapping:

  A named vector where names are names of `Graphs` in the resulting
  `Seurat` object and values are keys of `obsp` in `adata`. See below
  for details.

- misc_mapping:

  A named vector where names are names of `Misc` in the resulting
  `Seurat` object and values are keys of `uns` in `adata`. See below for
  details.

## Value

A `Seurat` object containing the requested data from `adata`

## Details

### Mapping arguments

All mapping arguments expect a named character vector where names are
the names of the slot in the `Seurat` object and values are the keys of
the corresponding slot of `adata`. If `TRUE`, the conversion function
will guess which items to copy as described in the conversion table
below. In most cases, the default is to copy all items using the same
names except where the correspondence between objects is unclear. The
`reduction_mapping` argument can also accept a more complex list format,
see below for details. To avoid copying anything to a slot, set the
mapping argument to `FALSE`. Empty mapping arguments (`NULL`,
[`c()`](https://rdrr.io/r/base/c.html),
[`list()`](https://rdrr.io/r/base/list.html)) will be treated as `FALSE`
with a warning. If an unnamed vector is provided, the values will be
used as names.

#### Examples:

- `TRUE` will guess which items to copy as described in the conversion
  table

- `c(seurat_item = "adata_item")` will copy `adata_item` from the slot
  in `adata` to `seurat_item` in the corresponding slot of the new
  `Seurat` object

- `FALSE` will avoid copying anything to the slot

- `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`

### Conversion table

|  |  |  |  |
|----|----|----|----|
| **From `AnnData`** | **To `Seurat`** | **Example mapping argument** | **Default** |
| `adata$X` | `Layers(seurat)` | `x_mapping = "counts"` *OR* `layers_mapping = c(counts = NA)` | The data in `adata$X` is copied to a layer named `X` |
| `adata$layers` | `Layers(seurat)` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
| `adata$obs` | `seurat[[]]` | `object_metadata_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
| `adata$var` | `seurat[[assay_name]][[]]` | `assay_metadata_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
| `adata$obsm` | `Reductions(seurat)` | `reduction_mapping = c(pca = "X_pca")` **OR** `reduction_mapping = list(pca = c(key = "PC_", embeddings = "X_pca", loadings = "PCs"))` | All items that can be coerced to a numeric matrix are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
| `adata$obsp` | `Graphs(seurat)` | `graph_mapping = c(nn = "connectivities")` | All items are copied by name |
| `adata$varp` | *NA* | *NA* | There is no corresponding slot for `varp` |
| `adata$uns` | `Misc(seurat)` | `misc_mapping = c(project_metadata = "metadata")` | All items are copied by name |

### The `reduction_mapping` argument

For the simpler named vector format, the names should be the names of
`Embeddings` in the resulting `Seurat` object, and the values should be
the keys of `obsm` in `adata`. A key will created from the `obsm` key.

For more advanced mapping, use the list format where each item is a
vector with the following names defining arguments to
[`SeuratObject::CreateDimReducObject()`](https://satijalab.github.io/seurat-object/reference/CreateDimReducObject.html):

- `key`: the key of the resulting
  [`SeuratObject::DimReduc`](https://satijalab.github.io/seurat-object/reference/DimReduc-class.html)
  object, passed to the `key` argument after processing with
  [`SeuratObject::Key()`](https://satijalab.github.io/seurat-object/reference/Key.html)

- `embeddings`: a key of the `obsm` slot in `adata`,
  `adata$obsm[[embeddings]]` is passed to the `embeddings` argument

- `loadings`: a key of the `varm` slot in `adata` (optional),
  `adata$varm[[loadings]]` is passed to the `loadings` argument

### The `x_mapping` and `layers_mapping` arguments

In order to specify where the data in `adata$X` will be stored in the
`Layers(seurat)` slot of the resulting object, you can use either the
`x_mapping` argument or the `layers_mapping` argument. If you use
`x_mapping`, it should be a string specifying the name of the layer in
`Layers(seurat)` where the data in `adata$X` will be stored. If you use
`layers_mapping`, it should be a named vector where names are names of
`Layers(seurat)` and values are keys of `layers` in `adata`. In order to
indicate the `adata$X` slot, you use `NA` as the value in the vector.
The name you provide for `x_mapping` may not be a name in
`layers_mapping`.

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md),
[`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md),
[`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md),
[`as_ZarrAnnData()`](https://anndataR.scverse.org/reference/as_ZarrAnnData.md),
[`reticulate-helpers`](https://anndataR.scverse.org/reference/reticulate-helpers.md)

## Examples

``` r
  ad <- AnnData(
    X = matrix(1:5, 3L, 5L),
    obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
    var = data.frame(row.names = letters[1:5], gene = 1:5)
  )

  # Default usage
  seurat <- ad$as_Seurat(
    assay_name = "RNA",
    x_mapping = "counts",
    layers_mapping = TRUE,
    object_metadata_mapping = TRUE,
    assay_metadata_mapping = TRUE,
    reduction_mapping = TRUE,
    graph_mapping = TRUE,
    misc_mapping = TRUE
  )
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
```
