# Read H5AD

Read data from a H5AD file

## Usage

``` r
read_h5ad(
  path,
  as = c("InMemoryAnnData", "HDF5AnnData", "SingleCellExperiment", "Seurat"),
  mode = c("r", "r+", "a", "w", "w-", "x"),
  backed = FALSE,
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

  The mode to open the HDF5 file.

  - `a` creates a new file or opens an existing one for read/write.

  - `r` opens an existing file for reading.

  - `r+` opens an existing file for read/write.

  - `w` creates a file, truncating any existing ones.

  - `w-`/`x` are synonyms, creating a file and failing if it already
    exists.

- backed:

  Whether to read the H5AD file in backed mode, returning an object
  containing
  [DelayedArray::DelayedMatrix](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  matrices. Which slots are backed depends on the value of `as`.

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
