# anndataR 0.99.1

## New features

* Add `example_data` function. 
* Harmonise all references to `h5ad_file` (`file`,`h5ad_path`,`ad`)
* Use `@returns` in Roxygen notes to allow multi-line notes.

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