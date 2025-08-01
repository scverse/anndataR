# anndataR devel

## anndataR 0.1.0.9011

- Updates for compatibility with Python **anndata** >= 0.12.0 (PR #305, Fixes #304)
  - Add helpers for reading/writing `NULL` values to/from H5AD files
  - Writing of `NULL` values can be disabled by setting
    `option(anndataR.write_null = FALSE)` to allow the files to be read by
    Python **anndata** < 0.12.0
- Fix a bug where string arrays were not transposed correctly when writing to
  H5AD files (PR #305)
- Fix a bug where the dimenions of dense arrays were not properly conserved
  when reading from H5AD (PR #305)
- Remove workarounds and skipping of `none` values in roundtrip tests (PR #305)

## anndataR 0.1.0.9010

- Switch HDF5 back end from **{hdf5r}** to **{rhdf5}** (PR #283, Fixes #272, #175, #299)
  - Includes improved compatibility with H5AD files written by Python **anndata**
- Improvements to roundtrip testing (PR #283)

## anndataR 0.1.0.9009

- Fix execution of roundtrip tests (PR #293)

## anndataR 0.1.0.9008

- Add Bioconductor installation instructions in preparation of submission (PR #297)

## anndataR 0.1.0.9007

- Refactor superfluous for loops (PR #298)
- 
## anndataR 0.1.0.9006

- Ignore `development_status.Rmd` vignette when building package (PR #296)

## anndataR 0.1.0.9005

- Bypass requiring a `counts` or `data` layer during `Seurat` conversion (PR #284)

## anndataR 0.1.0.9004

- Use accessors instead of direct slot access where possible (PR #291)

## anndataR 0.1.0.9003

- Simplify and update vignettes (PR #282)

## anndataR 0.1.0.9002

- Remove `anndataR.Rproj` file from repository (PR #292)

## anndataR 0.1.0.9001

- Change uses of `sapply()` to `vapply()` (PR #294)

# anndataR 0.1.0 (inital release candidate)

Initial release candidate of **{anndataR}** including:

- Native reading and writing of H5AD files
- R implementations of `InMemoryAnnData` and `HDF5AnnData` objects
- Conversion between `AnnData` and `SingleCellExperiment` or `Seurat` objects
- Extensive function documentation and vignettes demonstrating usage
- Comprehensive unit testing and identification of known issues
