# HDF5AnnData

Implementation of an HDF5-backed `AnnData` object. This class provides
an interface to a H5AD file and minimal data is stored in memory until
it is requested by the user. It is primarily designed as an intermediate
object when reading/writing H5AD files but can be useful for accessing
parts of large files.

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects.

## Value

An `HDF5AnnData` object

## See also

[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects

Other AnnData classes:
[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md),
[`AnnDataView`](https://anndataR.scverse.org/reference/AnnDataView.md),
[`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md),
[`ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.md),
[`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)

## Super class

[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md)
-\> `HDF5AnnData`

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

- [`HDF5AnnData$new()`](#method-HDF5AnnData-initialize)

- [`HDF5AnnData$obs_keys()`](#method-HDF5AnnData-obs_keys)

- [`HDF5AnnData$var_keys()`](#method-HDF5AnnData-var_keys)

- [`HDF5AnnData$layers_keys()`](#method-HDF5AnnData-layers_keys)

- [`HDF5AnnData$obsm_keys()`](#method-HDF5AnnData-obsm_keys)

- [`HDF5AnnData$varm_keys()`](#method-HDF5AnnData-varm_keys)

- [`HDF5AnnData$obsp_keys()`](#method-HDF5AnnData-obsp_keys)

- [`HDF5AnnData$varp_keys()`](#method-HDF5AnnData-varp_keys)

- [`HDF5AnnData$uns_keys()`](#method-HDF5AnnData-uns_keys)

- [`HDF5AnnData$close()`](#method-HDF5AnnData-close)

- [`HDF5AnnData$n_obs()`](#method-HDF5AnnData-n_obs)

- [`HDF5AnnData$n_vars()`](#method-HDF5AnnData-n_vars)

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

### `HDF5AnnData$new()`

Close the HDF5 file when the object is garbage collected

`HDF5AnnData` constructor

#### Usage

    HDF5AnnData$new(
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
      compression = c("none", "gzip", "lzf"),
      chunk_size = "auto"
    )

#### Arguments

- `file`:

  The file name (character) of the `.h5ad` file. If this file already
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

  The mode to open the HDF5 file. See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)
  for details

- `compression`:

  The compression algorithm to use. See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)
  for details

- `chunk_size`:

  The target chunk size in bytes. See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)
  for details

#### Details

The constructor creates a new HDF5 `AnnData` interface object. This can
either be used to either connect to an existing `.h5ad` file or to
create a new one. If any additional slot arguments are set an existing
file will be overwritten.

------------------------------------------------------------------------

### `HDF5AnnData$obs_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$obs_keys()

------------------------------------------------------------------------

### `HDF5AnnData$var_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$var_keys()

------------------------------------------------------------------------

### `HDF5AnnData$layers_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$layers_keys()

------------------------------------------------------------------------

### `HDF5AnnData$obsm_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$obsm_keys()

------------------------------------------------------------------------

### `HDF5AnnData$varm_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$varm_keys()

------------------------------------------------------------------------

### `HDF5AnnData$obsp_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$obsp_keys()

------------------------------------------------------------------------

### `HDF5AnnData$varp_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$varp_keys()

------------------------------------------------------------------------

### `HDF5AnnData$uns_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$uns_keys()

------------------------------------------------------------------------

### `HDF5AnnData$close()`

Close the HDF5 file

#### Usage

    HDF5AnnData$close()

------------------------------------------------------------------------

### `HDF5AnnData$n_obs()`

See the `n_obs` field in
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$n_obs()

------------------------------------------------------------------------

### `HDF5AnnData$n_vars()`

See the `n_vars` field in
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    HDF5AnnData$n_vars()
