# Read H5AD

Read data from a H5AD file

## Usage

``` r
read_h5ad(
  path,
  as = c("InMemoryAnnData", "HDF5AnnData", "SingleCellExperiment", "Seurat"),
  mode = c("r", "r+", "a", "w", "w-", "x"),
  backed = FALSE,
  group = "/",
  ...
)
```

## Arguments

- path:

  Path to the H5AD file to read

- as:

  The type of object to return. One of:

  - `"InMemoryAnnData"`: Read the H5AD file into memory as an
    [`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
    object

  - `"HDF5AnnData"`: Read the H5AD file as an
    [`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)
    object

  - `"SingleCellExperiment"`: Read the H5AD file as a
    [`SingleCellExperiment::SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
    object

  - `"Seurat"`: Read the H5AD file as a
    [`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
    object

- mode:

  The mode to open the HDF5 file. When `group` is not `"/"`, these apply
  to the group `group` rather than the whole file:

  - `a` creates a new file or opens an existing one for read/write. The
    group at `group` is created if it doesn't exist yet.

  - `r` opens an existing file for reading, the group `group` must
    already exist.

  - `r+` opens an existing file for read/write, the group `group` must
    already exist.

  - `w` creates a file, truncating any existing ones. When
    `group != "/"`, the file is created if missing but otherwise left
    untouched, and only the group `group` is deleted and recreated. HDF5
    doesn't free the space used by the deleted group, see the "File
    size" section of
    [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md).

  - `w-`/`x` are synonyms, creating a file and failing if it already
    exists. When `group != "/"`, they instead fail if the group `group`
    already exists.

- backed:

  Whether to read the H5AD file in backed mode, returning an object
  containing
  [DelayedArray::DelayedMatrix](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  matrices. Which slots are backed depends on the value of `as`.

- group:

  The path to the group within the HDF5 file that the `AnnData` is
  stored at. Defaults to `"/"`, the file root. Can be used to read an
  `AnnData` from a group inside a file that also contains other content,
  e.g. a modality inside a `.h5mu` file.

- ...:

  Extra arguments provided to the `as_*` conversion function for the
  object specified by `as`

## Value

The object specified by `as`

## See also

Other AnnData creators:
[`AnnData()`](https://anndataR.scverse.org/reference/AnnData.md),
[`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md),
[`read_zarr()`](https://anndataR.scverse.org/reference/read_zarr.md)

## Examples

``` r
h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")

# Read the H5AD as a SingleCellExperiment object
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment")
}

# Read the H5AD as a Seurat object
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  seurat <- read_h5ad(h5ad_file, as = "Seurat")
}
```
