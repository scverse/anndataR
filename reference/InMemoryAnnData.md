# InMemoryAnnData

Implementation of an in-memory `AnnData` object where data is stored
within the R session. This is the simplest back end and will be most
familiar to users. It is want you will want to use in most cases where
you want to interact with an `AnnData` object.

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects.

## Value

An `InMemoryAnnData` object

## See also

[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects

Other AnnData classes:
[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md),
[`AnnDataView`](https://anndataR.scverse.org/reference/AnnDataView.md),
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md),
[`ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.md),
[`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)

## Super class

[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md)
-\> `InMemoryAnnData`

## Active bindings

- `X`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `layers`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obs`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `var`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obs_names`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `var_names`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obsm`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `varm`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obsp`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `varp`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `uns`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

## Methods

### Public methods

- [`InMemoryAnnData$new()`](#method-InMemoryAnnData-initialize)

- [`InMemoryAnnData$clone()`](#method-InMemoryAnnData-clone)

Inherited methods

- [`AbstractAnnData$as_HDF5AnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_HDF5AnnData)
- [`AbstractAnnData$as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_InMemoryAnnData)
- [`AbstractAnnData$as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_ReticulateAnnData)
- [`AbstractAnnData$as_Seurat()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_Seurat)
- [`AbstractAnnData$as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_SingleCellExperiment)
- [`AbstractAnnData$as_ZarrAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_ZarrAnnData)
- [`AbstractAnnData$layers_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-layers_keys)
- [`AbstractAnnData$n_obs()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-n_obs)
- [`AbstractAnnData$n_vars()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-n_vars)
- [`AbstractAnnData$obs_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-obs_keys)
- [`AbstractAnnData$obsm_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-obsm_keys)
- [`AbstractAnnData$obsp_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-obsp_keys)
- [`AbstractAnnData$print()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-print)
- [`AbstractAnnData$shape()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-shape)
- [`AbstractAnnData$uns_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-uns_keys)
- [`AbstractAnnData$var_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-var_keys)
- [`AbstractAnnData$varm_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-varm_keys)
- [`AbstractAnnData$varp_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-varp_keys)
- [`AbstractAnnData$write_h5ad()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-write_h5ad)
- [`AbstractAnnData$write_zarr()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-write_zarr)

------------------------------------------------------------------------

### `InMemoryAnnData$new()`

Creates a new instance of an in-memory `AnnData` object. Inherits from
[AbstractAnnData](https://anndataR.scverse.org/reference/AbstractAnnData.md).

#### Usage

    InMemoryAnnData$new(
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

#### Arguments

- `X`:

  See the `X` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obs`:

  See the `obs` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `var`:

  See the `var` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `layers`:

  See the `layers` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obsm`:

  See the `obsm` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `varm`:

  See the `varm` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `obsp`:

  See the `obsp` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `varp`:

  See the `varp` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `uns`:

  See the `uns` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `shape`:

  Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if both `X` or
  `obs` and `var` are not provided.

------------------------------------------------------------------------

### `InMemoryAnnData$clone()`

The objects of this class are cloneable with this method.

#### Usage

    InMemoryAnnData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## complete example
ad <- AnnData(
  X = matrix(1:15, 3L, 5L),
  layers = list(
    A = matrix(5:1, 3L, 5L),
    B = matrix(letters[1:5], 3L, 5L)
  ),
  obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
  var = data.frame(row.names = letters[1:5], gene = 1:5)
)
ad
#> InMemoryAnnData object with n_obs × n_vars = 3 × 5
#>     obs: 'cell'
#>     var: 'gene'
#>     layers: 'A', 'B'

## minimum example
AnnData(
  obs = data.frame(row.names = letters[1:10]),
  var = data.frame(row.names = LETTERS[1:5])
)
#> InMemoryAnnData object with n_obs × n_vars = 10 × 5
```
