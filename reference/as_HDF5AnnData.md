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
  mode = c("w-", "r", "r+", "a", "w", "x"),
  backed = FALSE,
  group = "/"
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

  Whether the object is disk backed and returns
  [DelayedArray::DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  object for matrix data. Can only be `TRUE` when `mode == "r"`.

- group:

  The path to the group within the HDF5 file that the `AnnData` is
  stored at. Defaults to `"/"`, the file root. Can be used to write an
  `AnnData` to a group inside a file that also contains other content,
  e.g. a modality inside a `.h5mu` file.

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
