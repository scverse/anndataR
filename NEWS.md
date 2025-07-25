# anndataR 0.1.0.9009

- Fix execution of roundtrip tests
- 
# anndataR 0.1.0.9008

- Add Bioconductor installation instructions in preparation of submission

# anndataR 0.1.0.9007

- Refactor superfluous for loops
- 
# anndataR 0.1.0.9006

- ignore `development_status.Rmd` vignette when building package

# anndataR 0.1.0.9005

- Bypass requiring a `counts` or `data` layer during `Seurat` conversion

# anndataR 0.1.0.9004

- Use accessors instead of direct slot access where possible

# anndataR 0.1.0.9003

- Simplify & update vignetttes

# anndataR 0.1.0.9002

- remove `anndataR.Rproj` file from repository

# anndataR 0.1.0.9001

- change uses of `sapply` to `vapply`

# anndataR 0.1.0 (inital release candidate)

Initial release candidate of **{anndataR}** including:

- Native reading and writing of H5AD files
- R implementations of `InMemoryAnnData` and `HDF5AnnData` objects
- Conversion between `AnnData` and `SingleCellExperiment` or `Seurat` objects
- Extensive function documentation and vignettes demonstrating usage
- Comprehensive unit testing and identification of known issues
