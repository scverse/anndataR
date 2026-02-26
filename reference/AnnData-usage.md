# AnnData structure and usage

The `AnnData` object stores a data matrix `X` together with annotations
of observations `obs` (`obsm`, `obsp`) and variables `var` (`varm`,
`varp`). Additional layers of data can be stored in `layers` and
unstructured annotations in `uns`.

### Back ends

There are different back ends for `AnnData` objects that inherit from
the abstract
[AbstractAnnData](https://anndataR.scverse.org/reference/AbstractAnnData.md)
class and store and access data in different ways:

- [InMemoryAnnData](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
  stores data in memory

- [HDF5AnnData](https://anndataR.scverse.org/reference/HDF5AnnData.md)
  provides an interface to a H5AD file

- [ReticulateAnnData](https://anndataR.scverse.org/reference/ReticulateAnnData.md)
  wraps a Python `AnnData` object via reticulate

See the class documentation for details.

### Usage

The items listed as **"Slots"** are elements of the `AnnData` object
that contain data and can be accessed or set. **"Fields"** return
information about the `AnnData` object but cannot be set directly. Both,
as well as methods, can be accessed using the `$` operator

For example:

- `adata$X` returns the `X` matrix

- `adata$X <- x` sets the `X` matrix

- `adata$method()` calls a method

## Value

An `AnnData` object inheriting from
[`AbstractAnnData`](https://anndataR.scverse.org/reference/AbstractAnnData.md)

## Fields

- `shape`:

  Dimensions (observations x variables) of the `AnnData` object

- `n_obs`:

  Number of observations

- `n_vars`:

  Number of variables

- `obs_keys`:

  Keys (column names) of `obs`

- `var_keys`:

  Keys (column names) of `var`

- `layers_keys`:

  Keys (element names) of `layers`

- `obsm_keys`:

  Keys (element names) of `obsm`

- `varm_keys`:

  Keys (element names) of `varm`

- `obsp_keys`:

  Keys (element names) of `obsp`

- `varp_keys`:

  Keys (element names) of `varp`

- `uns_keys`:

  Keys (element names) of `uns`

## Slots

- `X`:

  The main data matrix. Either `NULL` or an observation x variable
  matrix (without dimnames) with dimensions consistent with `n_obs` and
  `n_vars`.

- `layers`:

  Additional data layers. Must be `NULL` or a named list of matrices
  having dimensions consistent with `n_obs` and `n_vars`.

- `obs`:

  Observation annotations. A `data.frame` with columns containing
  information about observations. The number of rows of `obs` defines
  the observation dimension of the `AnnData` object (`n_obs`). If
  `NULL`, an `n_obs` × 0 `data.frame` will automatically be generated.

- `var`:

  Variable annotations. A `data.frame` with columns containing
  information about variables. The number of rows of `var` defines the
  variable dimension of the `AnnData` object (`n_vars`). If `NULL`, an
  `n_vars` × 0 `data.frame` will automatically be generated.

- `obs_names`:

  Observation names. Either `NULL` or a vector of unique identifiers
  used to identify each row of `obs` and to act as an index into the
  observation dimension of the `AnnData` object. For compatibility with
  *R* representations, `obs_names` should be a unique character vector.

- `var_names`:

  Variable names. Either `NULL` or a vector of unique identifiers used
  to identify each row of `var` and to act as an index into the variable
  dimension of the `AnnData` object. For compatibility with *R*
  representations, `var_names` should be a unique character vector.

- `obsm`:

  Multi-dimensional observation annotation. Must be `NULL` or a named
  list of array-like elements with number of rows equal to `n_obs`.

- `varm`:

  Multi-dimensional variable annotations. Must be `NULL` or a named list
  of array-like elements with number of rows equal to `n_vars`.

- `obsp`:

  Observation pairs. Must be `NULL` or a named list of array-like
  elements with number of rows and columns equal to `n_obs`.

- `varp`:

  Variable pairs. Must be `NULL` or a named list of array-like elements
  with number of rows and columns equal to `n_vars`.

- `uns`:

  Unstructured annotations. Must be `NULL` or a named list.

## Methods

### Conversion methods:

- [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)
  :

  Convert to
  [`SingleCellExperiment::SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html),
  see
  [`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md)

- [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md):

  Convert to
  [`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html),
  see
  [`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md)

- [`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md):

  Convert to
  [`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md),
  as
  [`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md)

- [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md):

  Convert to
  [`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md),
  see
  [`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md)

- [`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md):

  Convert to
  [`ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.md),
  see
  [`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md)

### Output methods:

- [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)
  :

  Write the `AnnData` object to an HDF5 file, see
  [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)

### General methods:

- [`print()`](https://rdrr.io/r/base/print.html):

  Print a summary of the `AnnData` object

## Functions that can be used to create AnnData objects

- [`AnnData()`](https://anndataR.scverse.org/reference/AnnData.md):

  Create an
  [InMemoryAnnData](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
  object

- [`read_h5ad()`](https://anndataR.scverse.org/reference/read_h5ad.md):

  Read an `AnnData` from a H5AD file

- [`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md):

  Convert other objects to an `AnnData` object

## See also

The documentation for the Python `anndata` package
<https://anndata.readthedocs.io/en/stable/>

[AbstractAnnData](https://anndataR.scverse.org/reference/AbstractAnnData.md)
for the abstract class that all `AnnData` objects inherit from

[InMemoryAnnData](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
for the in-memory implementation of `AnnData`

[HDF5AnnData](https://anndataR.scverse.org/reference/HDF5AnnData.md) for
the HDF5-backed implementation of `AnnData`

[ReticulateAnnData](https://anndataR.scverse.org/reference/ReticulateAnnData.md)
for the reticulate-based implementation that wraps Python AnnData
objects
