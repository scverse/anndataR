# Contributing to anndataR

Thank you for your interest in contributing to anndataR!
Contributions of all kinds are welcome, including bug reports, documentation improvements and code.

## Reporting issues and asking questions

Please use the [issue tracker](https://github.com/scverse/anndataR/issues) and pick the matching template.

## Contributing code

1. Open an issue first for anything larger than a small fix
2. Fork the repository and create a branch from `devel`
3. Make your changes, following the conventions below
4. Open a pull request against `devel` and fill in the PR template

### Conventions

- Format code with [Air](https://posit-dev.github.io/air/) (`air format .`) and vignettes with `styler::style_dir("vignettes")`, then check with `lintr::lint_package()`
- Document functions with roxygen2 and regenerate the man pages with `devtools::document()`, using the roxygen2 version pinned in `DESCRIPTION` (`Config/roxygen2/version`)
- Add or update tests in `tests/testthat/`
- Add or update examples in the vignettes where relevant
- Check that new behaviour matches Python **anndata** where applicable (see the roundtrip tests)
- Keep lines to 120 characters or fewer (use `paste()` to wrap long strings in tests)
- Use `check_requires()` before using an optional package and add a matching skip in tests
- Use `cli` for messages and errors

### Running checks

```r
devtools::test()
devtools::test(filter = "roundtrip-X")  # Run a single test file
devtools::check()
BiocCheck::BiocCheck()
```

Some tests require Python **anndata** and are skipped if it is not available.
Set `ANNDATAR_TEST_REQUIRE_PYTHON=true` to make them fail instead.

### Testing

- Use the `skip_if_no_*()` helpers in `tests/testthat/` for tests needing optional dependencies (Python anndata, dummy_anndata, zarr, `h5diff`)
- Use `generate_dataset()` to create mock data with different matrix types and edge cases
- Add a roundtrip test (`test-roundtrip-*.R`) when changing anything that affects compatibility with Python
- Verify HDF5 changes against Python output with `h5diff`

## Package design

See the [software design article](https://anndataR.scverse.org/articles/software_design.html) for details.
Some key points to be aware of:

- **AnnData classes**: `AbstractAnnData` is the base class for `InMemoryAnnData`, `HDF5AnnData`, `ZarrAnnData`, `ReticulateAnnData` and `AnnDataView`.
  Changes to the interface usually need to be made in all of them.
- **Names are stored separately**: `obs_names` and `var_names` are stored apart from `obs`/`var`, and dimnames are added when data is accessed.
  Internal code should use `adata$obs_names` and `adata$var_names` rather than `rownames(adata$obs)` or `colnames(adata$X)`.
- **Views are read-only**: `AnnDataView` (returned by `adata[i, j]`) subsets lazily and can't be modified.
  Convert it with one of the `$as_*AnnData()` methods first when needed.
- **Validation**: setters use the validation helpers in `AbstractAnnData`.
  Update these if you change the rules for a slot.
- **Conversion**: `from_*()` functions convert other formats to AnnData and `as_*()` functions convert AnnData to other formats.
  Update both when adding a conversion feature.
- **Python differences**: Python AnnData allows duplicate names and distinguishes `None` from an empty dict, which R doesn't.
  Known issues are tracked in `inst/known_issues.yaml`.

### Where things live

anndataR implements the AnnData specification and the conversions, and leaves the file formats to other Bioconductor packages:

- **rhdf5** reads and writes HDF5, **Rarr** reads and writes Zarr.
  A bug in how data lands on disk (chunking, codecs, string encoding, metadata) belongs in their issue trackers; link it from the anndataR issue and bump the version pin in `DESCRIPTION` once it is fixed.
- **DelayedArray** provides backed (lazy) matrices, with **HDF5Array** seeds for H5AD and **ZarrArray** seeds for Zarr.
- **Matrix** provides the in-memory sparse matrices.
- **reticulate** and Python **anndata** are only needed for `ReticulateAnnData` and for the tests that compare against Python.

All of these are optional dependencies: guard new calls with `check_requires()` and skip the matching tests when the package isn't installed.

## Package authors

New contributors are welcome to add themselves to the package authors in `DESCRIPTION`.
Add a `person()` entry with the `ctb` role in your pull request, including your ORCID and GitHub username in the `comment` field.

## Code of conduct

Please be respectful and constructive in all interactions.
This project follows the [scverse code of conduct](https://github.com/scverse/.github/blob/main/CODE_OF_CONDUCT.md).
