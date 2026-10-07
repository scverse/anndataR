# Convert an `AnnData` to an `HDF5AnnData`

Convert another `AnnData` object to an
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)
object

## Usage

``` r
as_HDF5AnnData(
  adata,
  file,
  compression = c("none", "gzip", "lzf"),
  chunk_size = "auto",
  mode = c("w-", "r", "r+", "a", "w", "x")
)
```

## Arguments

- adata:

  An `AnnData` object to be converted to
  [`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)

- file:

  The file name (character) of the `.h5ad` file

- compression:

  The compression algorithm to use when writing the HDF5 file. Can be
  one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.

- chunk_size:

  The target chunk size in bytes to use when writing HDF5 datasets. When
  `"auto"` (default), the chunk size is determined automatically using
  an algorithm that mimics h5py's auto-chunking behaviour. Set to `NULL`
  to disable chunking (contiguous storage, the rhdf5 default), or a
  number to use a specific target size in bytes.

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
[`HDF5AnnData`](https://anndataR.scverse.org/reference/HDF5AnnData.md)
object with the same data as the input `AnnData` object.

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.scverse.org/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.scverse.org/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.scverse.org/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.scverse.org/reference/as_SingleCellExperiment.md),
[`as_ZarrAnnData()`](https://anndataR.scverse.org/reference/as_ZarrAnnData.md),
[`reticulate-helpers`](https://anndataR.scverse.org/reference/reticulate-helpers.md)
