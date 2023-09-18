# anndataR 0.99.0

## New features

* Add `example_data` function. 

## Bug fixes

* *DESCRIPTION*:
  - Add "email" tag to `Authors`.
  - Add `BiocViews`

# anndataR 0.1.0

## New features

Initial release of anndataR, providing support for working with AnnData objects in R.

Feature list:

* Slots:
  - X
  - layers
  - obs
  - obs_names
  - var
  - var_names

* Backends:
  - HDF5AnnData
  - InMemoryAnnData

* Converters:
  - SingleCellExperiment
  - Seurat