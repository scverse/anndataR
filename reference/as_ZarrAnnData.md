# Convert an `AnnData` to an `ZarrAnnData`

Convert another `AnnData` object to an
[`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)
object

## Usage

``` r
as_ZarrAnnData(
  adata,
  file,
  compression = c("none", "gzip", "blosc", "zstd", "lzma", "bz2", "zlib", "lz4"),
  mode = c("w-", "r", "r+", "a", "w", "x"),
  zarr_format = NULL
)
```

## Arguments

- adata:

  An `AnnData` object to be converted to
  [`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)

- file:

  The file name (character) of the `.zarr` file

- compression:

  The compression algorithm to use when writing the Zarr file. Can be
  one of `"none"`, `"gzip"`, `"blosc"`, `"zstd"`, `"lzma"`, `"bz2"`,
  `"zlib"` or `"lz4"`. Defaults to `"none"`.

- mode:

  The mode to open the Zarr file:

  - `a` creates a new file or opens an existing one for read/write

  - `r` opens an existing file for reading

  - `r+` opens an existing file for read/write

  - `w` creates a file, truncating any existing ones

  - `w-`/`x` are synonyms, creating a file and failing if it already
    exists

- zarr_format:

  The format to use when creating the Zarr store. Should be either 2 or
  3 for Zarr v2 or v3 formats, respectively. The default can also be set
  using the `anndataR.zarr_format` option, e.g.
  `options(anndataR.zarr_format = 2)`. If neither is set, Zarr v3 is
  used.

  When an existing store is opened, that store's format is used instead
  and setting `zarr_format` to a different format is an error. Use
  `mode = "w"` to rewrite an existing store in another format.

## Value

A [`ZarrAnnData`](https://anndataR.scverse.org/reference/ZarrAnnData.md)
object with the same data as the input `AnnData` object.

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md),
[`as_HDF5AnnData()`](https://anndataR.scverse.org/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md),
[`reticulate-helpers`](https://anndataR.scverse.org/reference/reticulate-helpers.md)
