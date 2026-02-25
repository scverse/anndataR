# anndataR Benchmarks

Performance benchmarks for anndataR using [`bench::mark()`](https://bench.r-lib.org/) with continuous tracking on [bencher.dev](https://bencher.dev/perf/anndatar).

## Prerequisites

**R packages:**

```r
install.packages(c("bench", "jsonlite", "reticulate"))
```

**Python packages** (used to generate canonical test H5AD files via [dummy-anndata](https://pypi.org/project/dummy-anndata/)):

```bash
pip install anndata dummy-anndata
```

## Running locally

From the repository root:

```bash
# Run all suites with default settings (2000 obs × 1000 vars, 3 iterations)
Rscript benchmarks/run_benchmarks.R

# Quick smoke test
Rscript benchmarks/run_benchmarks.R --n-obs 100 --n-vars 50 --iterations 1

# Run a single suite
Rscript benchmarks/run_benchmarks.R --suite read

# Custom output path
Rscript benchmarks/run_benchmarks.R --output /tmp/results.json

# Show all options
Rscript benchmarks/run_benchmarks.R --help
```

### CLI options

| Option           | Default                   | Description                                                             |
| ---------------- | ------------------------- | ----------------------------------------------------------------------- |
| `--n-obs N`      | 2000                      | Number of observations (cells)                                          |
| `--n-vars N`     | 1000                      | Number of variables (genes)                                             |
| `--iterations N` | 3                         | Iterations per benchmark                                                |
| `--suite NAME`   | all                       | Suite to run: `all`, `read`, `write`, `get`, `set`, `convert`, `subset` |
| `--output FILE`  | `benchmarks/results.json` | Output JSON path                                                        |

## Output format

Results are written as [Bencher Metric Format (BMF)](https://bencher.dev/docs/reference/bencher-metric-format/) JSON:

```json
{
  "read_InMemory_float_csparse": {
    "latency": {
      "value": 266384594.08,
      "lower_value": 254901470.03,
      "upper_value": 266384594.08
    },
    "memory": {
      "value": 1770408
    }
  }
}
```

- **latency**: median / min / max in nanoseconds
- **memory**: peak allocation in bytes

## Benchmark suites

### `read` — Reading H5AD files

Benchmarks `read_h5ad()` across X matrix types and backends.

| Benchmark                | Description                      |
| ------------------------ | -------------------------------- |
| `read_InMemory_{x_type}` | Read H5AD into `InMemoryAnnData` |
| `read_HDF5_{x_type}`     | Read H5AD into `HDF5AnnData`     |

X types: `float_matrix`, `float_csparse`, `float_rsparse`, `integer_csparse`

### `write` — Writing H5AD files

Benchmarks `write_h5ad()` from both backends with compression variants.

| Benchmark                                | Description   |
| ---------------------------------------- | ------------- |
| `write_{backend}_{x_type}_{compression}` | Write to H5AD |

Backends: `InMemory`, `HDF5`. Compressions: `none`, `gzip`.

### `get` — Slot getters

Benchmarks reading every AnnData slot, key methods, dimension methods, and S3 methods.

| Benchmark                   | Description                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `get_{backend}_{slot}`      | Access a slot (`X`, `obs`, `var`, `layers`, `obsm`, `varm`, `obsp`, `varp`, `uns`, `obs_names`, `var_names`)                                      |
| `method_{backend}_{method}` | Call a method (`obs_keys`, `var_keys`, `obsm_keys`, `varm_keys`, `obsp_keys`, `varp_keys`, `layers_keys`, `uns_keys`, `n_obs`, `n_vars`, `shape`) |
| `s3_{backend}_{method}`     | S3 generic (`dim`, `nrow`, `ncol`, `dimnames`, `rownames`, `colnames`)                                                                            |

### `set` — Slot setters

Benchmarks writing every AnnData slot on both backends.

| Benchmark              | Description      |
| ---------------------- | ---------------- |
| `set_{backend}_{slot}` | Set a slot value |

### `convert` — Format conversions

| Benchmark                           | Description                 |
| ----------------------------------- | --------------------------- |
| `convert_HDF5_to_InMemory_{x_type}` | `as_InMemoryAnnData()`      |
| `convert_InMemory_to_HDF5_{x_type}` | `as_HDF5AnnData()`          |
| `convert_InMemory_to_SCE`           | `as_SingleCellExperiment()` |
| `convert_SCE_to_InMemory`           | `as_AnnData()` from SCE     |
| `convert_InMemory_to_Seurat`        | `as_Seurat()`               |
| `convert_Seurat_to_InMemory`        | `as_AnnData()` from Seurat  |

### `subset` — View creation and materialisation

| Benchmark                           | Description                        |
| ----------------------------------- | ---------------------------------- |
| `subset_obs_int_{backend}`          | `ad[1:n, ]` integer index          |
| `subset_both_int_{backend}`         | `ad[1:n, 1:p]` integer index       |
| `subset_obs_char_{backend}`         | `ad[names, ]` character index      |
| `subset_both_char_{backend}`        | `ad[names, names]` character index |
| `subset_obs_logical_{backend}`      | `ad[logical, ]` logical index      |
| `materialize_to_InMemory_{backend}` | `view$as_InMemoryAnnData()`        |
| `materialize_to_HDF5_{backend}`     | `view$as_HDF5AnnData()`            |

## Test data generation

H5AD files are generated using Python's `dummy_anndata` package and written by Python `anndata`. This ensures the benchmark input files have canonical HDF5 encoding, independent of anndataR's own writer — so any encoding bugs in anndataR cannot silently affect benchmark data.

Each file contains:

- **X**: matrix of the specified type
- **layers**: `float_csparse`, `integer_csparse`
- **obs**: `categorical`, `dense_array`, `string_array`
- **var**: `string_array`, `boolean_array`, `dense_array`
- **obsm**: `float_matrix`, `float_csparse`
- **varm**: `float_matrix`
- **obsp / varp**: `float_csparse`
- **uns**: `string`, `integer`, `float` (with nested dict)

Generated files are cached in `$TMPDIR/anndatar_bench/` for the session.

## CI integration

Three GitHub Actions workflows handle continuous benchmarking:

### `bench-base.yaml` — Baseline tracking

- **Trigger**: push to `devel`
- **Purpose**: Tracks baseline performance on [bencher.dev](https://bencher.dev/perf/anndatar)
- Uses a t-test with upper boundary 0.99 to detect regressions

### `bench-pr-run.yaml` — PR benchmark runner

- **Trigger**: pull request (opened, reopened, synchronize)
- **Purpose**: Runs benchmarks on PR code. Fork-safe — no secrets are needed.
- Uploads `benchmark_results.json` and `event.json` as artifacts

### `bench-pr-track.yaml` — PR result tracker

- **Trigger**: `workflow_run` (after PR run completes)
- **Purpose**: Downloads artifacts and submits results to bencher.dev. Runs in the default branch context so repository secrets are available.
- Posts regression alerts as PR comments

### Required secrets

| Secret              | Description                                       |
| ------------------- | ------------------------------------------------- |
| `BENCHER_API_TOKEN` | API token from [bencher.dev](https://bencher.dev) |

`GITHUB_TOKEN` is provided automatically.

## File structure

    benchmarks/
    ├── README.md              # This file
    ├── run_benchmarks.R       # Main entry point / CLI runner
    ├── lib/
    │   └── helpers.R          # Data generation, bench→BMF conversion, CLI parsing
    └── suites/
        ├── bench_read.R       # read_h5ad benchmarks
        ├── bench_write.R      # write_h5ad benchmarks
        ├── bench_get.R        # Slot getter benchmarks
        ├── bench_set.R        # Slot setter benchmarks
        ├── bench_convert.R    # Conversion benchmarks
        └── bench_subset.R     # Subsetting / view benchmarks
