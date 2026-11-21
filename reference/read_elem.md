# Read an element

Read a single element from an H5AD file or Zarr store

## Usage

``` r
read_elem(path, name, backed = FALSE, format = c("auto", "h5ad", "zarr"), ...)
```

## Arguments

- path:

  Path to the H5AD file or Zarr store to read from

- name:

  Path of the element inside the file or store, e.g. `"obs"` or
  `"obsm/X_pca"`

- backed:

  Whether to return array and matrix elements as
  [DelayedArray::DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  objects instead of reading them into memory.

- format:

  The format of `path`. One of:

  - `"auto"`: Detect the format from `path`. A file that is an HDF5 file
    is read as H5AD and a directory with Zarr metadata is read as Zarr.

  - `"h5ad"`: Read `path` as an H5AD file

  - `"zarr"`: Read `path` as a Zarr store

- ...:

  Not currently used. Unused arguments are an error.

## Value

The element, as the matching R type. An `anndata` group is returned as a
read-only
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)
or
[`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)
object.

## See also

Other element readers and writers:
[`write_elem()`](https://anndataR.scverse.org/reference/write_elem.md)

## Examples

``` r
h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")

if (requireNamespace("rhdf5", quietly = TRUE)) {
  # Read the obs data frame
  obs <- read_elem(h5ad_file, "obs")

  # Read a single obsm matrix
  pca <- read_elem(h5ad_file, "obsm/X_pca")
}
```
