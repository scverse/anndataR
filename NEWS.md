# anndataR unreleased

## New functionality

* Implemented an `AnnDataView` class, which provides a lazy view of an `AnnData` object without copying data (PR #1096)
* Implemented S3 methods for `AbstractAnnData` objects: `dim`, `nrow`, `ncol`, `dimnames`, `rownames`, `colnames`, and `[` (PR #1096)
* Add `ReticulateAnnData` class for seamless Python integration via **{reticulate}** (PR #322)

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

## Documentation

- `citation("anndataR")` now returns details of the **{anndataR}** preprint
  (PR #351)

# anndataR 0.99.2

* Add `@return` to R6 object man pages to address **{BiocCheck}** warning
  (PR #319)

# anndataR 0.99.1

* Resubmit because of failing Bioc checks (PR #318)

# anndataR 0.99.0

- Bump required R version to 4.5 in preparation for Bioconductor submission
  (PR #309)
- Remove deprecated functions and arguments (PR #311, PR #313)
- Minor cleanups and improvements (PR #313).
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

# anndataR 0.1.0 (inital release candidate)

Initial release candidate of **{anndataR}** including:

- Native reading and writing of H5AD files
- R implementations of `InMemoryAnnData` and `HDF5AnnData` objects
- Conversion between `AnnData` and `SingleCellExperiment` or `Seurat` objects
- Extensive function documentation and vignettes demonstrating usage
- Comprehensive unit testing and identification of known issues
