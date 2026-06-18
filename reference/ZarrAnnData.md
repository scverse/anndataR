# ZarrAnnData

Implementation of a Zarr-backed `AnnData` object. This class provides an
interface to a Zarr file and minimal data is stored in memory until it
is requested by the user. It is primarily designed as an intermediate
object when reading/writing Zarr files but can be useful for accessing
parts of large files.

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects.

## Value

A `ZarrAnnData` object

## See also

[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects

Other AnnData classes:
[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md),
[`AnnDataView`](https://anndataR.scverse.org/reference/AnnDataView.md),
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md),
[`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md),
[`ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.md)

## Super class

[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md)
-\> `ZarrAnnData`

## Active bindings

- `X`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

- `layers`:

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

- `uns`:

  See
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

## Methods

### Public methods

- [`ZarrAnnData$new()`](#method-ZarrAnnData-initialize)

- [`ZarrAnnData$n_obs()`](#method-ZarrAnnData-n_obs)

- [`ZarrAnnData$n_vars()`](#method-ZarrAnnData-n_vars)

- [`ZarrAnnData$obs_keys()`](#method-ZarrAnnData-obs_keys)

- [`ZarrAnnData$var_keys()`](#method-ZarrAnnData-var_keys)

- [`ZarrAnnData$layers_keys()`](#method-ZarrAnnData-layers_keys)

- [`ZarrAnnData$obsm_keys()`](#method-ZarrAnnData-obsm_keys)

- [`ZarrAnnData$varm_keys()`](#method-ZarrAnnData-varm_keys)

- [`ZarrAnnData$obsp_keys()`](#method-ZarrAnnData-obsp_keys)

- [`ZarrAnnData$varp_keys()`](#method-ZarrAnnData-varp_keys)

- [`ZarrAnnData$uns_keys()`](#method-ZarrAnnData-uns_keys)

Inherited methods

- [`AbstractAnnData$as_HDF5AnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_HDF5AnnData)
- [`AbstractAnnData$as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_InMemoryAnnData)
- [`AbstractAnnData$as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_ReticulateAnnData)
- [`AbstractAnnData$as_Seurat()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_Seurat)
- [`AbstractAnnData$as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_SingleCellExperiment)
- [`AbstractAnnData$as_ZarrAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_ZarrAnnData)
- [`AbstractAnnData$print()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-print)
- [`AbstractAnnData$shape()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-shape)
- [`AbstractAnnData$write_h5ad()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-write_h5ad)
- [`AbstractAnnData$write_zarr()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-write_zarr)

------------------------------------------------------------------------

### `ZarrAnnData$new()`

`ZarrAnnData` constructor

#### Usage

    ZarrAnnData$new(
      file,
      X = NULL,
      obs = NULL,
      var = NULL,
      layers = NULL,
      obsm = NULL,
      varm = NULL,
      obsp = NULL,
      varp = NULL,
      uns = NULL,
      shape = NULL,
      mode = c("a", "r", "r+", "w", "w-", "x"),
      compression = c("none", "gzip", "blosc", "zstd", "lzma", "bz2", "zlib", "lz4")
    )

#### Arguments

- `file`:

  The file name (character) of the `.zarr` file. If this file already
  exits, other arguments must be `NULL`.

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

- `mode`:

  The mode to open the Zarr file. See
  [`as_ZarrAnnData()`](https://anndataR.scverse.org/reference/as_ZarrAnnData.md)
  for details

- `compression`:

  The compression algorithm to use. See
  [`as_ZarrAnnData()`](https://anndataR.scverse.org/reference/as_ZarrAnnData.md)
  for details

#### Details

The constructor creates a new Zarr `AnnData` interface object. This can
either be used to either connect to an existing `.zarr` file or to
create a new one. If any additional slot arguments are set an existing
file will be overwritten.

------------------------------------------------------------------------

### `ZarrAnnData$n_obs()`

See the `n_obs` field in
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$n_obs()

------------------------------------------------------------------------

### `ZarrAnnData$n_vars()`

See the `n_vars` field in
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$n_vars()

------------------------------------------------------------------------

### `ZarrAnnData$obs_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$obs_keys()

------------------------------------------------------------------------

### `ZarrAnnData$var_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$var_keys()

------------------------------------------------------------------------

### `ZarrAnnData$layers_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$layers_keys()

------------------------------------------------------------------------

### `ZarrAnnData$obsm_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$obsm_keys()

------------------------------------------------------------------------

### `ZarrAnnData$varm_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$varm_keys()

------------------------------------------------------------------------

### `ZarrAnnData$obsp_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$obsp_keys()

------------------------------------------------------------------------

### `ZarrAnnData$varp_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$varp_keys()

------------------------------------------------------------------------

### `ZarrAnnData$uns_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ZarrAnnData$uns_keys()
