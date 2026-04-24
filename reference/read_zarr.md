# Read Zarr

Read data from a Zarr store

## Usage

``` r
read_zarr(
  path,
  as = c("InMemoryAnnData", "ZarrAnnData", "SingleCellExperiment", "Seurat"),
  mode = c("r", "r+", "a", "w", "w-", "x"),
  ...
)
```

## Arguments

- path:

  Path to the Zarr store to read

- as:

  The type of object to return. One of:

  - `"InMemoryAnnData"`: Read the Zarr store into memory as an
    [`InMemoryAnnData`](https://anndataR.scverse.org/reference/InMemoryAnnData.md)
    object

  - `"ZarrAnnData"`: Read the Zarr store as an
    [`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)
    object

  - `"SingleCellExperiment"`: Read the Zarr store as a
    [`SingleCellExperiment::SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
    object

  - `"Seurat"`: Read the Zarr store as a
    [`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
    object

- mode:

  The mode to open the Zarr file.

  - `a` creates a new file or opens an existing one for read/write.

  - `r` opens an existing file for reading.

  - `r+` opens an existing file for read/write.

  - `w` creates a file, truncating any existing ones.

  - `w-`/`x` are synonyms, creating a file and failing if it already
    exists.

- ...:

  Extra arguments provided to the `as_*` conversion function for the
  object specified by `as`

## Value

The object specified by `as`

## See also

Other AnnData creators:
[`AnnData()`](https://anndataR.scverse.org/reference/AnnData.md),
[`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md),
[`read_h5ad()`](https://anndataR.scverse.org/reference/read_h5ad.md)

## Examples

``` r
# Please use "example_v3.zarr.zip" for AnnData stored as Zarr version 3
zarr_dir <- system.file("extdata", "example_v2.zarr.zip", package = "anndataR")
td <- tempdir(check = TRUE)
unzip(zarr_dir, exdir = td)
zarr_store <- file.path(td, "example_v2.zarr")

# Read the Zarr as a SingleCellExperiment object
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  sce <- read_zarr(zarr_store, as = "SingleCellExperiment")
}

# Read the Zarr as a Seurat object
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  seurat <- read_zarr(zarr_store, as = "Seurat")
}
```
