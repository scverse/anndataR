# ReticulateAnnData

Implementation of an `AnnData` object that wraps a Python **anndata**
`AnnData` object using reticulate. This allows direct interaction with
Python `AnnData` objects while maintaining the R interface. It is useful
when you already have a Python `AnnData` or to access functionality that
has not yet been implemented in anndataR.

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects.

## Value

A `ReticulateAnnData` object

## See also

[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects

Other AnnData classes:
[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md),
[`AnnDataView`](https://anndataR.scverse.org/reference/AnnDataView.md),
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md),
[`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md)

## Super class

[`anndataR::AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md)
-\> `ReticulateAnnData`

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

- [`ReticulateAnnData$new()`](#method-ReticulateAnnData-new)

- [`ReticulateAnnData$n_obs()`](#method-ReticulateAnnData-n_obs)

- [`ReticulateAnnData$n_vars()`](#method-ReticulateAnnData-n_vars)

- [`ReticulateAnnData$py_anndata()`](#method-ReticulateAnnData-py_anndata)

Inherited methods

- [`anndataR::AbstractAnnData$as_HDF5AnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_HDF5AnnData)
- [`anndataR::AbstractAnnData$as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_InMemoryAnnData)
- [`anndataR::AbstractAnnData$as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_ReticulateAnnData)
- [`anndataR::AbstractAnnData$as_Seurat()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_Seurat)
- [`anndataR::AbstractAnnData$as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-as_SingleCellExperiment)
- [`anndataR::AbstractAnnData$layers_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-layers_keys)
- [`anndataR::AbstractAnnData$obs_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-obs_keys)
- [`anndataR::AbstractAnnData$obsm_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-obsm_keys)
- [`anndataR::AbstractAnnData$obsp_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-obsp_keys)
- [`anndataR::AbstractAnnData$print()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-print)
- [`anndataR::AbstractAnnData$shape()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-shape)
- [`anndataR::AbstractAnnData$uns_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-uns_keys)
- [`anndataR::AbstractAnnData$var_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-var_keys)
- [`anndataR::AbstractAnnData$varm_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-varm_keys)
- [`anndataR::AbstractAnnData$varp_keys()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-varp_keys)
- [`anndataR::AbstractAnnData$write_h5ad()`](https://anndataR.scverse.org/reference/AbstractAnnData.html#method-write_h5ad)

------------------------------------------------------------------------

### Method `new()`

`ReticulateAnnData` constructor

#### Usage

    ReticulateAnnData$new(
      py_anndata = NULL,
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

- `py_anndata`:

  A Python AnnData object created using reticulate, or NULL to create a
  new empty Python AnnData object

- `X`:

  See the `X` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `obs`:

  See the `obs` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `var`:

  See the `var` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `layers`:

  See the `layers` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `obsm`:

  See the `obsm` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `varm`:

  See the `varm` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `obsp`:

  See the `obsp` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `varp`:

  See the `varp` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `uns`:

  See the `uns` slot in
  [AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
  (only used if py_anndata is NULL)

- `shape`:

  Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if both `X` or
  `obs` and `var` are not provided. (only used if py_anndata is NULL)

#### Details

The constructor creates a new ReticulateAnnData interface object that
wraps a Python AnnData object. If `py_anndata` is provided, it must be a
valid Python AnnData object. If NULL, a new Python AnnData object will
be created using the other provided arguments.

------------------------------------------------------------------------

### Method `n_obs()`

See the `n_obs` field in
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ReticulateAnnData$n_obs()

------------------------------------------------------------------------

### Method `n_vars()`

See the `n_vars` field in
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    ReticulateAnnData$n_vars()

------------------------------------------------------------------------

### Method `py_anndata()`

Get the underlying Python AnnData object

#### Usage

    ReticulateAnnData$py_anndata()

#### Returns

The Python AnnData object wrapped by this ReticulateAnnData
