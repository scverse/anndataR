# Agent instructions for anndataR

Read [CONTRIBUTING.md](CONTRIBUTING.md) first; the conventions, checks and key design points live there and are not repeated here.
This file adds what an agent needs on top: which package owns what, where the direction is recorded, and the things in this codebase that aren't obvious from reading it.

## What anndataR is

The native R implementation of AnnData: it reads and writes `.h5ad` files and Zarr stores without Python, and converts to and from SingleCellExperiment and Seurat.
It is released on Bioconductor and part of scverse.
Python **anndata** is the reference implementation: when R and Python disagree about a file, Python's behaviour is the spec unless an issue says otherwise.

## Before reporting a finding

Most of what looks wrong when reading the code is already tracked, and some of it is being replaced by an open PR.
Before raising something, check:

- open issues: `gh issue list --search "<keyword>"`; useful labels are `h5ad`, `zarr`, `backed`, `Seurat`, `SingleCellExperiment`, `performance`, `design`, `triage`
- open PRs: `gh pr list`; several change architecture, and a finding in code a PR rewrites should be made against that PR, not against `devel`
- the [project board](https://github.com/orgs/scverse/projects/67): this is the live roadmap, and this file deliberately doesn't carry one
- `inst/known_issues.yaml`: expected test failures, skipped with `message_if_known()`

If it is tracked, link the issue or PR instead of re-reporting it.
Upstream trackers matter too: a bug in a dependency is often already filed on `Huber-group-EMBL/Rarr` or `Huber-group-EMBL/rhdf5`.

## Who owns what

anndataR sits on a stack of Bioconductor packages, and more and more of what it used to do itself is their job.
The split:

| Concern | Owner | Where it touches anndataR |
| --- | --- | --- |
| The AnnData spec: element encodings, `encoding-type`/`encoding-version`, slot validation | anndataR | `R/read_*_helpers.R`, `R/write_*_helpers.R`, the validators in `AbstractAnnData` |
| The AnnData classes and lazy views | anndataR | `R/*AnnData.R`, `R/AnnDataView.R` |
| Conversion to and from SingleCellExperiment and Seurat | anndataR | `R/as_*.R`, `R/from_*.R` |
| Raw HDF5 reading and writing | **rhdf5** | `R/HDF5File.R` wraps the handle, `R/write_hdf5_helpers.R` |
| Raw Zarr reading and writing: chunks, codecs, string encoding, consolidated metadata | **Rarr** | `R/Rarr_utils.R`, `zarr_write_compressed()` |
| Backed (lazy, on-disk) matrices | **DelayedArray**, with **HDF5Array** seeds for H5AD and **ZarrArray** seeds for Zarr (#507) | `backed = TRUE` in the readers; `as_SingleCellExperiment()` and `as_Seurat()` pass the seeds through |
| In-memory sparse matrices | **Matrix** (`dgCMatrix`, `dgRMatrix`); SparseArray is being evaluated in #382 | everywhere |
| Python interop | **reticulate** and Python **anndata** | `ReticulateAnnData`, `R/reticulate_helpers.R` |

What follows from this:

- **Format-level bugs go upstream.** If the problem is how bytes land on disk (a chunk layout, a codec name, string encoding, metadata), file it on rhdf5 or Rarr, link it from the anndataR issue, work around it here with a `TODO` naming the upstream issue, and bump the version pin in `DESCRIPTION` once the fix is released. Recent examples: Rarr#262 (the C string encoder, merged), Rarr#111 (`patch_zarr_vlen_utf8()` waits on it), rhdf5#241 (a spurious `NA` warning, tracked in `inst/known_issues.yaml`).
- **Prefer the owner's implementation.** anndataR decides *what* to write (the spec) and passes *how* through. When the owner doesn't provide something yet, as with auto-chunking in Rarr, implement it here pinned to Python's behaviour and propose it upstream.
- **Everything below anndataR is optional.** These packages are in `Suggests` and loaded through `check_requires(what, requires, where, version)`. A new call into one of them needs a `check_requires()` guard and a matching `skip_if_*()` in the tests.

## Direction

These are decided; assume them when a change touches the area.

- **Backed mode goes through DelayedArray seeds for both formats.** HDF5Array today, ZarrArray for Zarr once #507 lands. Backed mode covers every slot, not only `X` as in Python: a `DelayedArray` in `obsm` or `uns` is intentional.
- **H5AD and Zarr share one encoding layer.** `R/read_h5ad_helpers.R` and `R/read_zarr_helpers.R` (and their `write_` counterparts) mirror each other, and #445 is about merging them. A new encoding is added to both; don't let them drift further apart.
- **Parity with Python anndata 0.13** (#514): sharded Zarr v3 stores, new encodings, and the missing API (`concat` #325, `copy` #352, `read_elem`/`write_elem` #466).
- **Performance is measured, not guessed.** Bencher runs on every PR (`.github/workflows/bench-*.yaml`); #329 and #382 are the long-running threads on memory use.

Open questions, not direction; don't build on them yet:

- `S4Vectors::DataFrame` instead of `data.frame` for `obs`/`var` (#554).
- Moving the conversions to zellkonverter (#383): if zellkonverter takes ownership of the AnnData <> SingleCellExperiment/Seurat conversions and builds on `AbstractAnnData`, rather than on its own Python stack.

## Things that aren't obvious from the code

- `obs_names` and `var_names` are stored apart from `obs`/`var` and from every matrix; dimnames are added on access by `.add_matrix_dimnames()` and friends. Internal code uses `adata$obs_names`, never `rownames(adata$obs)` or `colnames(adata$X)`, and setters strip row names before storing. Python keeps the names separate too, which is what makes round trips exact.
- `adata[i, j]` returns an `AnnDataView`: lazy, read-only, converted with `$as_InMemoryAnnData()` and the other `$as_*AnnData()` methods.
- The S3 methods (`dim()`, `dimnames()`, `[`) for every class live in `R/AbstractAnnData-s3methods.R`; `n_obs()`/`n_vars()` are the R6 equivalents.
- `HDF5AnnData` opens the file per operation (`open_and_defer_close()`) and closes it when the object is garbage collected; `test-h5ad-fileclosure.R` guards against leaked handles.
- R data frames can't hold duplicate row names, Python can. `generate_dataset(obs_has_row_names = FALSE)` covers that edge.
- `NULL` in `uns` is written as an empty dataset, as anndata >= 0.12 does; `options(anndataR.write_null = FALSE)` turns that off.
- Every getter and setter carries a `# trackstatus: class=..., feature=..., status=...` comment; the development-status article is generated from them, so a new feature gets one.
- `man/` is generated. Never hand-edit it; run `devtools::document()` with the roxygen2 version pinned in `DESCRIPTION`.

## Testing against Python

- The `skip_if_no_*()` helpers in `tests/testthat/helper-skip_if_*.R` declare the Python requirement with `reticulate::py_require()` and skip when it can't be imported. `ANNDATAR_TEST_REQUIRE_PYTHON=true` turns the skips into errors; CI and scverse/integration-testing set it, so a broken Python environment can't silently skip every Python test.
- `test-roundtrip-*.R` write with R and read with Python and vice versa, over every element type of `dummy_anndata`, for H5AD, Zarr v2 and Zarr v3. A change to any encoding gets a roundtrip test.
- `expect_equal_py()` compares R and Python objects; `generate_dataset()` builds the R fixtures.
- `h5diff` (skipped when absent) checks that H5AD files written by R and Python are the same on disk.
