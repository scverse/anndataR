# Convert an `AnnData` to an `HDF5AnnData`

Convert another `AnnData` object to an
[`HDF5AnnData`](https://anndataR.data-intuitive.com/reference/HDF5AnnData.md)
object

## Usage

``` r
as_HDF5AnnData(
  adata,
  file,
  compression = c("none", "gzip", "lzf"),
  mode = c("w-", "r", "r+", "a", "w", "x")
)
```

## Arguments

- adata:

  An `AnnData` object to be converted to
  [`HDF5AnnData`](https://anndataR.data-intuitive.com/reference/HDF5AnnData.md)

- file:

  The file name (character) of the `.h5ad` file

- compression:

  The compression algorithm to use when writing the HDF5 file. Can be
  one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.

- mode:

  The mode to open the HDF5 file:

  - `a` creates a new file or opens an existing one for read/write

  - `r` opens an existing file for reading

  - `r+` opens an existing file for read/write

  - `w` creates a file, truncating any existing ones

  - `w-`/`x` are synonyms, creating a file and failing if it already
    exists

## Value

An
[`HDF5AnnData`](https://anndataR.data-intuitive.com/reference/HDF5AnnData.md)
object with the same data as the input `AnnData` object.

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.data-intuitive.com/reference/as_AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.data-intuitive.com/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.data-intuitive.com/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.data-intuitive.com/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.data-intuitive.com/reference/as_SingleCellExperiment.md),
[`reticulate-helpers`](https://anndataR.data-intuitive.com/reference/reticulate-helpers.md)
