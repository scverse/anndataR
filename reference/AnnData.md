# Create an in-memory AnnData object.

For more information on the functionality of an AnnData object, see
[AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

## Usage

``` r
AnnData(
  X = NULL,
  obs = NULL,
  var = NULL,
  layers = NULL,
  obsm = NULL,
  varm = NULL,
  obsp = NULL,
  varp = NULL,
  uns = NULL,
  shape = NULL
)
```

## Arguments

- X:

  See the `X` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- obs:

  See the `obs` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- var:

  See the `var` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- layers:

  See the `layers` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- obsm:

  See the `obsm` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- varm:

  See the `varm` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- obsp:

  See the `obsp` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- varp:

  See the `varp` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- uns:

  See the `uns` slot in
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- shape:

  Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if both `X` or
  `obs` and `var` are not provided.

## Value

An
[InMemoryAnnData](https://anndataR.data-intuitive.com/reference/InMemoryAnnData.md)
object

## See also

[AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)
for details of `AnnData` structure and usage

Other AnnData creators:
[`as_AnnData()`](https://anndataR.data-intuitive.com/reference/as_AnnData.md),
[`read_h5ad()`](https://anndataR.data-intuitive.com/reference/read_h5ad.md)

## Examples

``` r
adata <- AnnData(
  X = matrix(1:12, nrow = 3, ncol = 4),
  obs = data.frame(
    row.names = paste0("obs", 1:3),
    n_counts = c(1, 2, 3),
    n_cells = c(1, 2, 3)
  ),
  var = data.frame(
    row.names = paste0("var", 1:4),
    n_cells = c(1, 2, 3, 4)
  )
)

adata
#> InMemoryAnnData object with n_obs × n_vars = 3 × 4
#>     obs: 'n_counts', 'n_cells'
#>     var: 'n_cells'
```
