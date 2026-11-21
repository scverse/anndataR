# Write an element

Write a single element to an H5AD file or Zarr store

## Usage

``` r
write_elem(
  path,
  name,
  value,
  overwrite = FALSE,
  format = c("auto", "h5ad", "zarr"),
  ...
)
```

## Arguments

- path:

  Path to the H5AD file or Zarr store to write to. Created if it doesn't
  exist.

- name:

  Path of the element inside the file or store. Use `"/"` to write an
  `AnnData` object to the whole file or store.

- value:

  The value to write. An `AnnData` object is written as an `anndata`
  group.

- overwrite:

  Whether to replace an element that already exists at `name`. If
  `FALSE` (default), an existing element is an error.

- format:

  The format of `path`. One of:

  - `"auto"`: Detect the format from `path`. An existing HDF5 file is
    written as H5AD and an existing directory with Zarr metadata as
    Zarr. A new path is detected from its extension.

  - `"h5ad"`: Write `path` as an H5AD file

  - `"zarr"`: Write `path` as a Zarr store

- ...:

  Options for writing the element:

  - `compression`: The compression to use. See
    [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md)
    and
    [`write_zarr()`](https://anndataR.scverse.org/reference/write_zarr.md).

  - `chunk_size`: The target chunk size. See
    [`write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.md).

## Value

`path` (invisibly)

## See also

Other element readers and writers:
[`read_elem()`](https://anndataR.scverse.org/reference/read_elem.md)

## Examples

``` r
h5ad_file <- tempfile(fileext = ".h5ad")

if (requireNamespace("rhdf5", quietly = TRUE)) {
  # Write a data frame to a new file
  write_elem(h5ad_file, "uns/df", data.frame(a = 1:3, b = letters[1:3]))

  # Replace it with a vector
  write_elem(h5ad_file, "uns/df", 1:10, overwrite = TRUE)

  read_elem(h5ad_file, "uns/df")
}
#>  [1]  1  2  3  4  5  6  7  8  9 10
```
