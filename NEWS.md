# anndataR 1.3.1

- Add support for reading backed objects using `DelayedArray` matrices, including a `backed` argument for `HDF5AnnData` and `InMemoryAnnData`, and conversion of backed `AnnData` objects to `SingleCellExperiment`/`Seurat` (PR #387).
- Add support for writing Zarr v3 stores, selected with the `zarr_format` argument of `write_zarr()`/`as_ZarrAnnData()` or the `anndataR.zarr_format` option. New stores are written as Zarr v3 by default (PR #455).
- Write strings as VLen-UTF8 rather than as fixed-width, so that string arrays read back as variable length strings in Python `anndata` (PR #455).
- Use the format of an existing Zarr store when writing to it and truncate the store when `mode = "w"`, so that a store can no longer end up with a mix of Zarr v2 and v3 nodes (PR #455).
- Enable additional linters and optimise several suboptimal code patterns (PR #453).
- Fix `create_zarr()` so that Zarr stores can be created at relative paths and at paths containing regex metacharacters (PR #478).
- Fix use of `&` instead of `&&` when validating column names in `AbstractAnnData` (PR #487).
- Update tests to take into account the `Rarr` structured datatype breaking change (PR #462).
- Temporarily require Python `anndata < 0.13.0` in tests, workflows and the Python usage vignette (PR #481).
- Add Zarr to the benchmarks (PR #446).
- Add Bioconductor R-universe checks (PR #470).
- Fix CI runners (PR #457).
- Upgrade to roxygen2 8.0.0.
- Simplify roxygen documentation for the H5AD and Zarr helpers by using `@inheritParams` (PR #479).
- Update `CITATION` with the published paper (PR #464).
- Mark Zarr support as implemented in the class diagram (PR #460).

# anndataR 1.3.0

- Bioconductor 3.24 devel

# anndataR 1.2.0

- Bioconductor 3.23 release

# anndataR 1.1.3

- Prepare for Bioconductor 3.23 release (PR #449).
- Add initial Zarr support for reading and writing `.zarr` stores (PR #190).

# anndataR 1.1.2

- Fix broken links (PR #429).
- Fix the software design vignette (PR #438).
- Fix GitHub Actions MacOS/Windows setup issue (PR #447).

# anndataR 1.1.1

- Fix CI (PR #418).
- Add continuous benchmarking using bencher (PR #423, PR #425).
- Handle unnamed `SingleCellExperiment` assays in `as_AnnData()` by automatically assigning names with a warning (PR #420).
- Optimise sparse matrix reading performance by avoiding `Matrix::sparseMatrix` and constructing objects manually (PR #417).
- Allow manually setting chunk size for HDF5 writes (PR #424).
- Add auto-chunking for HDF5 writes to improve performance (PR #424).

# anndataR 1.1.0

- Bioconductor 3.23 devel
- Separate reading and writing keys in `HDF5AnnData` (PR #392).
- Add a check for 0 dimensions when converting to `Seurat` (PR #393).
- Warn that matrix dim names will not be written to H5AD files (PR #396).
- Add a warning to `.from_SCE_convert()` (PR #402).
- Warn when a rotation is provided but not stored in a `LinearEmbeddingMatrix` (PR #407).
- Fix linting issues (PR #389).
- Add Bioconductor badges (PR #375).
- Update **{pkgdown}** rendering (PR #379).
- Force loading the environment in the Python vignette (PR #394).
- Add issue templates (PR #395).
- Update GitHub Actions based on those used by **{rhdf5}** (PR #378).
- Set different cache versions in GitHub Actions (PR #403).
- Update CI (PR #416).

# anndataR 1.0.2

- Fix broken links (PR #430).

# anndataR 1.0.1

- Update author email address (PR #410).
- Fix linting issues (PR #412).
- Fix roundtrip test; ensure that `X` is always 2-dimensional (PR #413).
- Fix HDF5 write test; make sure to only pass CsparseMatrix to Seurat (PR #413).

# anndataR 1.0.0

- Bioconductor 3.22 release (October 2025)

# anndataR 0.99.6

- Use the correct mapping figure in the `SingleCellExperiment` usage vignette (PR #372)
- Fix accessing layers keys for `ReticulateAnnData` (PR #372)
- Minor updates to the Python usage vignette (PR #372)
- Add the Python vignette to the **{pkgdown}** index (PR #372)

# anndataR 0.99.4

- Address minor issues from Bioconductor checks (PR #371)
  - Combine the `opts_chunks$set()` calls in the Python usage vignette
  - Replace `\donttrun` with `\donttest` in man pages

# anndataR 0.99.3

## New functionality

- Implemented an `AnnDataView` class, which provides a lazy view of an `AnnData` object without copying data (PR #324)
- Implemented S3 methods for `AbstractAnnData` objects: `dim`, `nrow`, `ncol`, `dimnames`, `rownames`, `colnames`, and `[` (PR #324)
- Add `ReticulateAnnData` class for seamless Python integration via **{reticulate}** (PR #322)
- Add `AGENTS.md` with instructions for developers (PR #367).

## Major changes

- Refactor obs/var_names handling for improved data consistency (PR #328)
  - `InMemoryAnnData` now stores `obs_names` and `var_names` as separate private
    fields instead of relying on rownames of `obs`/`var` data.frames
  - `HDF5AnnData` maintains separate obs/var names management to ensure
    consistency between obs/var data.frames and dimnames
  - All matrix data (`X`, `layers`, `obsm`, `varm`, `obsp`, `varp`) is now
    stored internally **without** dimnames for consistency
  - Dimnames are added **on-the-fly** when users access data, ensuring proper
    obs/var name display

## Minor changes

- Refactor setter methods in `HDF5AnnData` and `InMemoryAnnData` to use pipe 
  operators for cleaner code (PR #328)
- Add explanatory comments for matrix generation alignment with Python
  **dummy-anndata** (PR #328)
- Add `get_generator_types()` function return allowed/example types for
  `generate_dataset()` (PR #354)
- Add checks for type arguments to `generate_dataset()` (PR #354)
- Generalise the layers created by `generate_dataset()` when `format = "Seurat"`
  (PR #354)

## Bug fixes

- Add compression parameter to additional write operations in `HDF5AnnData` 
  for consistency (PR #328)
- Directly use `obs_names` and `var_names` properties instead of corresponding
  indirect S3 methods `rownames` and `colnames` (PR #328)
- Fix error message variable name in `.validate_aligned_array()` method 
  (expected_colnames → expected_rownames) (PR #328)
- Fix Seurat conversion for PCA loadings with variable feature subsets (PR #328)
  - Seurat PCA loadings only contain variable features, not all genes
  - Now properly expands loadings matrix to include all genes with zeros for
    non-variable features
  - Adds warning when rownames don't match var_names during conversion
- Avoid writing character datasets to H5AD files with LZF compression as it
  causes R to crash (PR #356)
- Handle slots that may have incomplete dimensions when converting from
  `Seurat`. These are now skipped with a warning instead of indirectly raising
  an error. (PR #369)

## Documentation

- `citation("anndataR")` now returns details of the **{anndataR}** preprint
  (PR #351)
- Update vignettes to clarify and expand text and improve formatting (PR #360)

# anndataR 0.99.2

* Add `@return` to R6 object man pages to address **{BiocCheck}** warning
  (PR #319)

# anndataR 0.99.1

* Resubmit because of failing Bioc checks (PR #318)

# anndataR 0.99.0

- Bump required R version to 4.5 in preparation for Bioconductor submission
  (PR #309)
- Remove deprecated functions and arguments (PR #311, PR #313)
- Minor clean-ups and improvements (PR #313).
- Add more tests to increase coverage (PR #315)

# anndataR 0.2.0

## Breaking changes

- Switch the HDF5 back end to use the **{rhdf5}** package instead of **{hdf5r}**
  (PR #283, Fixes #272, #175, #299)
  - This addresses various issues related to H5AD files and allows better
    integration with Bioconductor. Most of the previous known issues have now
    been resolved.
  - It also greatly improves compatibility with H5AD files written by Python
    **anndata**
  - **NOTE:** Make sure to install **{rhdf5}** instead of **{hdf5r}** to be able
    to read and write H5AD files!

## Major changes

- Updates for compatibility with Python **anndata** >= 0.12.0 (PR #305,
  Fixes #304)
  - Add helpers for reading/writing `NULL` values to/from H5AD files
  - Writing of `NULL` values can be disabled by setting
    `option(anndataR.write_null = FALSE)` to allow the files to be read by
    Python **anndata** < 0.12.0
- A `counts` or `data` layer is no longer required during `Seurat` conversion
  (PR #284)
  - There will still be a warning if neither of this is present as it may
    affect compatibility with **{Seurat}** functions
    
## Minor changes

- Use accessor functions/methods instead of direct slot access where possible
  (PR #291)
- Refactor superfluous for loops (PR #298)
- Change uses of `sapply()` to `vapply()` (PR #294)
- Ignore `development_status.Rmd` vignette when building package (PR #296)
- Remove `anndataR.Rproj` file from repository (PR #292)

## Bug fixes

- Fix a bug where string arrays were not transposed correctly when writing to
  H5AD files (PR #305)
- Fix a bug where the dimensions of dense arrays were not properly conserved
  when reading from H5AD (PR #305)

## Documentation

- Simplify and update vignettes (PR #282)
- Add Bioconductor installation instructions in preparation for submission (PR #297)

## Testing

- Improvements to round trip testing (PR #283, PR #293, PR #305)
  - Most round trip tests are now enabled and pass successfully
  - Conversion helpers have been added to assist with **{reticulate}** tests

# anndataR 0.1.0 (initial release candidate)

Initial release candidate of **{anndataR}** including:

- Native reading and writing of H5AD files
- R implementations of `InMemoryAnnData` and `HDF5AnnData` objects
- Conversion between `AnnData` and `SingleCellExperiment` or `Seurat` objects
- Extensive function documentation and vignettes demonstrating usage
- Comprehensive unit testing and identification of known issues
