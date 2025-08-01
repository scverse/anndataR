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
