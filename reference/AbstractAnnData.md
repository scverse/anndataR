# Abstract AnnData class

This class is an abstract representation of an `AnnData` object. It is
intended to be used as a base class for concrete implementations of
`AnnData` objects, such as
[InMemoryAnnData](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
or [HDF5AnnData](https://anndataR.scverse.org/reference/HDF5AnnData.md).

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects.

## Value

An `AbstractAnnData` object

## See also

[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects

Other AnnData classes:
[`AnnDataView`](https://anndataR.scverse.org/reference/AnnDataView.md),
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md),
[`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md),
[`ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.md)

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

- [`AbstractAnnData$print()`](#method-AbstractAnnData-print)

- [`AbstractAnnData$shape()`](#method-AbstractAnnData-shape)

- [`AbstractAnnData$n_obs()`](#method-AbstractAnnData-n_obs)

- [`AbstractAnnData$n_vars()`](#method-AbstractAnnData-n_vars)

- [`AbstractAnnData$obs_keys()`](#method-AbstractAnnData-obs_keys)

- [`AbstractAnnData$var_keys()`](#method-AbstractAnnData-var_keys)

- [`AbstractAnnData$layers_keys()`](#method-AbstractAnnData-layers_keys)

- [`AbstractAnnData$obsm_keys()`](#method-AbstractAnnData-obsm_keys)

- [`AbstractAnnData$varm_keys()`](#method-AbstractAnnData-varm_keys)

- [`AbstractAnnData$obsp_keys()`](#method-AbstractAnnData-obsp_keys)

- [`AbstractAnnData$varp_keys()`](#method-AbstractAnnData-varp_keys)

- [`AbstractAnnData$uns_keys()`](#method-AbstractAnnData-uns_keys)

- [`AbstractAnnData$as_SingleCellExperiment()`](#method-AbstractAnnData-as_SingleCellExperiment)

- [`AbstractAnnData$as_Seurat()`](#method-AbstractAnnData-as_Seurat)

- [`AbstractAnnData$as_InMemoryAnnData()`](#method-AbstractAnnData-as_InMemoryAnnData)

- [`AbstractAnnData$as_ReticulateAnnData()`](#method-AbstractAnnData-as_ReticulateAnnData)

- [`AbstractAnnData$as_HDF5AnnData()`](#method-AbstractAnnData-as_HDF5AnnData)

- [`AbstractAnnData$write_h5ad()`](#method-AbstractAnnData-write_h5ad)

- [`AbstractAnnData$clone()`](#method-AbstractAnnData-clone)

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$print(...)

#### Arguments

- `...`:

  Optional arguments to print method

------------------------------------------------------------------------

### Method `shape()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$shape()

------------------------------------------------------------------------

### Method `n_obs()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$n_obs()

------------------------------------------------------------------------

### Method `n_vars()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$n_vars()

------------------------------------------------------------------------

### Method `obs_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$obs_keys()

------------------------------------------------------------------------

### Method `var_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$var_keys()

------------------------------------------------------------------------

### Method `layers_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$layers_keys()

------------------------------------------------------------------------

### Method `obsm_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$obsm_keys()

------------------------------------------------------------------------

### Method `varm_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$varm_keys()

------------------------------------------------------------------------

### Method `obsp_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$obsp_keys()

------------------------------------------------------------------------

### Method `varp_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$varp_keys()

------------------------------------------------------------------------

### Method `uns_keys()`

See
[AnnData-usage](https://anndataR.scverse.org/reference/AnnData-usage.md)

#### Usage

    AbstractAnnData$uns_keys()

------------------------------------------------------------------------

### Method [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

Convert to `SingleCellExperiment`

See
[`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)
for more details on the conversion

#### Usage

    AbstractAnnData$as_SingleCellExperiment(
      x_mapping = NULL,
      assays_mapping = TRUE,
      colData_mapping = TRUE,
      rowData_mapping = TRUE,
      reducedDims_mapping = TRUE,
      colPairs_mapping = TRUE,
      rowPairs_mapping = TRUE,
      metadata_mapping = TRUE
    )

#### Arguments

- `x_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `assays_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `colData_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `rowData_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `reducedDims_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `colPairs_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `rowPairs_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- `metadata_mapping`:

  See
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

#### Returns

A `SingleCellExperiment` object

------------------------------------------------------------------------

### Method [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

Convert to `Seurat`

See [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)
for more details on the conversion

#### Usage

    AbstractAnnData$as_Seurat(
      assay_name = "RNA",
      x_mapping = NULL,
      layers_mapping = TRUE,
      object_metadata_mapping = TRUE,
      assay_metadata_mapping = TRUE,
      reduction_mapping = TRUE,
      graph_mapping = TRUE,
      misc_mapping = TRUE
    )

#### Arguments

- `assay_name`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `x_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `layers_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `object_metadata_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `assay_metadata_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `reduction_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `graph_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- `misc_mapping`:

  See
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

#### Returns

A `Seurat` object

------------------------------------------------------------------------

### Method [`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md)

Convert to an
[`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md)

See
[`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md)
for more details on the conversion

#### Usage

    AbstractAnnData$as_InMemoryAnnData()

#### Returns

An
[InMemoryAnnData](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
object

------------------------------------------------------------------------

### Method [`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md)

Convert to a
[`ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.md)

See
[`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md)
for more details on the conversion

#### Usage

    AbstractAnnData$as_ReticulateAnnData()

#### Returns

A
[ReticulateAnnData](https://anndataR.scverse.org/reference/ReticulateAnnData.md)
object

------------------------------------------------------------------------

### Method [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)

Convert to an
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)

See
[`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)
for more details on the conversion

#### Usage

    AbstractAnnData$as_HDF5AnnData(
      file,
      compression = c("none", "gzip", "lzf"),
      chunk_size = "auto",
      mode = c("w-", "r", "r+", "a", "w", "x")
    )

#### Arguments

- `file`:

  See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)

- `compression`:

  See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)

- `chunk_size`:

  See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)

- `mode`:

  See
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)

#### Returns

An
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)
object

------------------------------------------------------------------------

### Method [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)

Write the `AnnData` object to an H5AD file

See
[`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)
for details

#### Usage

    AbstractAnnData$write_h5ad(
      path,
      compression = c("none", "gzip", "lzf"),
      chunk_size = "auto",
      mode = c("w-", "r", "r+", "a", "w", "x")
    )

#### Arguments

- `path`:

  See
  [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)

- `compression`:

  See
  [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)

- `chunk_size`:

  See
  [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)

- `mode`:

  See
  [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)

#### Returns

`path` invisibly

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AbstractAnnData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
