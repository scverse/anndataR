# Convert an `AnnData` to a `ReticulateAnnData`

Convert another `AnnData` object to a
[`ReticulateAnnData`](https://anndataR.data-intuitive.com/reference/ReticulateAnnData.md)
object

## Usage

``` r
as_ReticulateAnnData(adata)
```

## Arguments

- adata:

  An `AnnData` object to be converted to
  [`ReticulateAnnData`](https://anndataR.data-intuitive.com/reference/ReticulateAnnData.md)

## Value

A
[`ReticulateAnnData`](https://anndataR.data-intuitive.com/reference/ReticulateAnnData.md)
object with the same data as the input `AnnData` object

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.data-intuitive.com/reference/as_AnnData.md),
[`as_HDF5AnnData()`](https://anndataR.data-intuitive.com/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.data-intuitive.com/reference/as_InMemoryAnnData.md),
[`as_Seurat()`](https://anndataR.data-intuitive.com/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.data-intuitive.com/reference/as_SingleCellExperiment.md),
[`reticulate-helpers`](https://anndataR.data-intuitive.com/reference/reticulate-helpers.md)

## Examples

``` r
# \donttest{
# Requires Python anndata to be installed
if (requireNamespace("reticulate", quietly = TRUE) &&
      reticulate::py_module_available("anndata")) {
  ad <- AnnData(
    X = matrix(1:5, 3L, 5L),
    layers = list(
      A = matrix(5:1, 3L, 5L),
      B = matrix(letters[1:5], 3L, 5L)
    ),
    obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
    var = data.frame(row.names = letters[1:5], gene = 1:5)
  )
  ad$as_ReticulateAnnData()
}
#> ReticulateAnnData object with n_obs × n_vars = 3 × 5
#>     layers: 'A', 'B'
# }
```
