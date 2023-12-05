# anndataR 0.99.0

## New features

* PR #158: Change package version from 0.0.0.9000 --> 0.99.0 to align with Bioc devel 
    versioning standards.
  - Modify NEWS file to reflect changes in the package versioning.

* Various PRs: Initial release of anndataR, providing support for working with
  AnnData objects in R. Feature list:
  - Slots
    - X
    - layers
    - obs
    - obs_names
    - var
    - var_names
  - Backends:
    - HDF5AnnData
    - InMemoryAnnData
  - Converters:
    - SingleCellExperiment
    - Seurat
    
  ## Bug fixes
  