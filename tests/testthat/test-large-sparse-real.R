# The other large matrix tests lower the limit on the number of non-zero values,
# see test-large-sparse.R. These use a matrix that really has more than 2^31 - 1
# non-zero values. This needs about 60 GB of memory and 15 GB of disk space, and
# takes about half an hour, so they only run if asked to:
#
#   ANNDATAR_TEST_LARGE=true ANNDATAR_TEST_LARGE_DIR=/path/to/cache \
#     Rscript -e 'devtools::test(filter = "large-sparse-real")'
#
# R doesn't give memory back to the system between tests, so run them in separate
# processes by setting ANNDATAR_TEST_LARGE_PART to one of "read_csr", "read_csc",
# "read_zarr_v2", "read_zarr_v3", "read_backed", "write_h5ad", "reread_h5ad",
# "write_zarr_v2", "reread_zarr_v2", "write_zarr_v3", "reread_zarr_v3",
# "convert_sce", "convert_seurat", "python_to_r" or "r_to_python". Without it, all
# of them run, which needs more memory than the above.
skip_if(
  !identical(Sys.getenv("ANNDATAR_TEST_LARGE"), "true"),
  "ANNDATAR_TEST_LARGE is not 'true'"
)
skip_if_not_installed("SparseArray")
skip_if_not_installed("BiocGenerics")
skip_if_not_installed("MatrixGenerics")
skip_if_not_installed("rhdf5")
skip_if_not_installed("Rarr", minimum_version = "2.1.45")
skip_if_no_anndata_py()
skip_if_no_zarr()

large_part <- Sys.getenv("ANNDATAR_TEST_LARGE_PART", "all")
skip_unless_part <- function(part) {
  skip_if(!large_part %in% c("all", part), paste("Not running part", part))
}

cache_dir <- Sys.getenv("ANNDATAR_TEST_LARGE_DIR", tempfile("anndataR-large"))
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Row i of the CSR matrix has `per_major` values, at columns 3k + i %% 3 and with
# values (i %% 7) + 1 + (k %% 5), so that any value can be checked afterwards.
# The CSC files hold the same arrays with the axes swapped.
n_major <- 250000L
n_minor <- 30000L
per_major <- 8800L
nnz <- n_major * as.double(per_major)
stopifnot(nnz > .Machine$integer.max)

expected_major <- function(i) {
  # i is the 1-based major element, the returned value holds the first 9 minor
  # elements
  i0 <- i - 1L
  out <- numeric(9L)
  k <- 0:8
  columns <- 3L * k + i0 %% 3L
  in_range <- columns < 9L
  out[columns[in_range] + 1L] <- (i0 %% 7L) + 1 + (k[in_range] %% 5L)
  out
}
# The first and last elements, and those around offset 2^31 in the arrays
test_majors <- c(1L, 2L, 244032L, 244033L, 244034L, 249999L, 250000L)

make_large_fixtures <- function(dir) {
  reticulate::py_run_string(
    "
import numpy as np, h5py, zarr, pandas as pd, anndata as ad, warnings

N_MAJOR, N_MINOR, PER = 250_000, 30_000, 8_800
NNZ = N_MAJOR * PER
BLOCK = 2100
CHUNK = 2**22

k = np.arange(PER, dtype=np.int32)
blk_i = np.concatenate([3 * k + i % 3 for i in range(BLOCK)]).astype(np.int32)
blk_d = np.concatenate([(i % 7) + 1.0 + (k % 5) for i in range(BLOCK)])
indptr = np.arange(N_MAJOR + 1, dtype=np.int64) * PER

def blocks():
    done = 0
    while done < N_MAJOR:
        n = min(BLOCK, N_MAJOR - done)
        yield n * PER, blk_i[: n * PER], blk_d[: n * PER]
        done += n

def empty(n_obs, n_var):
    return ad.AnnData(
        obs=pd.DataFrame(index=[f'c{i}' for i in range(n_obs)]),
        var=pd.DataFrame(index=[f'g{i}' for i in range(n_var)]),
    )

def make_h5ad(path, kind):
    n_obs, n_var = (N_MAJOR, N_MINOR) if kind == 'csr' else (N_MINOR, N_MAJOR)
    empty(n_obs, n_var).write_h5ad(path)
    with h5py.File(path, 'a') as f:
        if 'X' in f:
            del f['X']
        g = f.create_group('X')
        g.attrs.update({
            'encoding-type': f'{kind}_matrix',
            'encoding-version': '0.1.0',
            'shape': np.array([n_obs, n_var], dtype=np.int64),
        })
        opts = dict(compression='gzip', compression_opts=1)
        g.create_dataset('indptr', data=indptr, **opts)
        di = g.create_dataset('indices', shape=(NNZ,), dtype='i4', chunks=(CHUNK,), **opts)
        dd = g.create_dataset('data', shape=(NNZ,), dtype='f8', chunks=(CHUNK,), **opts)
        pos = 0
        for n, i_, d_ in blocks():
            di[pos:pos + n] = i_
            dd[pos:pos + n] = d_
            pos += n

def make_zarr(path, zarr_format):
    ad.settings.zarr_write_format = zarr_format
    empty(N_MAJOR, N_MINOR).write_zarr(path)
    root = zarr.open_group(path, mode='a')
    if 'X' in root:
        del root['X']
    g = root.create_group('X')
    g.attrs.update({
        'encoding-type': 'csr_matrix',
        'encoding-version': '0.1.0',
        'shape': [N_MAJOR, N_MINOR],
    })
    g.create_array('indptr', data=indptr, chunks=(CHUNK,))
    ai = g.create_array('indices', shape=(NNZ,), dtype='i4', chunks=(CHUNK,))
    ax = g.create_array('data', shape=(NNZ,), dtype='f8', chunks=(CHUNK,))
    pos = 0
    for n, i_, d_ in blocks():
        ai[pos:pos + n] = i_
        ax[pos:pos + n] = d_
        pos += n
"
  )
  main <- reticulate::import_main(convert = FALSE)
  paths <- fixture_paths(dir)
  if (!file.exists(paths$csr)) {
    main$make_h5ad(paths$csr, "csr")
  }
  if (!file.exists(paths$csc)) {
    main$make_h5ad(paths$csc, "csc")
  }
  if (!dir.exists(paths$zarr_v2)) {
    main$make_zarr(paths$zarr_v2, 2L)
  }
  if (!dir.exists(paths$zarr_v3)) {
    main$make_zarr(paths$zarr_v3, 3L)
  }
  invisible()
}

fixture_paths <- function(dir) {
  list(
    csr = file.path(dir, "csr.h5ad"),
    csc = file.path(dir, "csc.h5ad"),
    zarr_v2 = file.path(dir, "csr_v2.zarr"),
    zarr_v3 = file.path(dir, "csr_v3.zarr")
  )
}

make_large_fixtures(cache_dir)
paths <- fixture_paths(cache_dir)

expect_large_csr <- function(X) {
  expect_identical(dim(X), c(n_major, n_minor))
  expect_identical(SparseArray::nzcount(X), nnz)
  for (i in test_majors) {
    expect_equal(
      as.vector(as.matrix(X[i, 1:9, drop = FALSE])),
      expected_major(i),
      info = i
    )
  }
}

expect_large_csc <- function(X) {
  expect_identical(dim(X), c(n_minor, n_major))
  expect_identical(SparseArray::nzcount(X), nnz)
  for (i in test_majors) {
    expect_equal(
      as.vector(as.matrix(X[1:9, i, drop = FALSE])),
      expected_major(i),
      info = i
    )
  }
}

readers <- list(
  read_csr = function() read_h5ad(paths$csr),
  read_zarr_v2 = function() read_zarr(paths$zarr_v2),
  read_zarr_v3 = function() read_zarr(paths$zarr_v3)
)
for (part in names(readers)) {
  test_that(paste(part, "reads more than 2^31 - 1 non-zero values"), {
    skip_unless_part(part)
    expect_no_warning(adata <- readers[[part]]())
    expect_s4_class(adata$X, "SVT_SparseMatrix")
    expect_large_csr(adata$X)
  })
}

test_that("read_csc reads more than 2^31 - 1 non-zero values", {
  skip_unless_part("read_csc")
  expect_no_warning(adata <- read_h5ad(paths$csc))
  expect_s4_class(adata$X, "SVT_SparseMatrix")
  expect_large_csc(adata$X)
})

test_that("backed reads of more than 2^31 - 1 non-zero values work", {
  skip_unless_part("read_backed")
  skip_if_not_installed("HDF5Array")

  X <- read_h5ad(paths$csr, as = "HDF5AnnData", backed = TRUE)$X
  expect_s4_class(X, "DelayedMatrix")
  expect_identical(dim(X), c(n_major, n_minor))
  for (i in test_majors) {
    expect_equal(
      as.vector(as.matrix(X[i, 1:9, drop = FALSE])),
      expected_major(i),
      info = i
    )
  }
})

test_that("a large Zarr can be read backed", {
  skip_unless_part("read_backed")
  skip_if_not_installed("ZarrArray")

  for (zarr_format in c("zarr_v2", "zarr_v3")) {
    X <- read_zarr(paths[[zarr_format]], as = "ZarrAnnData", backed = TRUE)$X
    expect_s4_class(X, "DelayedMatrix")
    expect_identical(dim(X), c(n_major, n_minor))
    for (i in test_majors) {
      expect_equal(
        as.vector(as.matrix(X[i, 1:9, drop = FALSE])),
        expected_major(i),
        info = paste(zarr_format, i)
      )
    }
  }
})

# Written files are kept for the "reread_*" parts, which remove them
writers <- list(
  h5ad = list(
    path = file.path(cache_dir, "written.h5ad"),
    write = function(adata, path) {
      write_h5ad(adata, path, mode = "w", compression = "gzip")
    },
    read = function(path) read_h5ad(path),
    check = "check_h5ad"
  ),
  zarr_v2 = list(
    path = file.path(cache_dir, "written_v2.zarr"),
    write = function(adata, path) write_zarr(adata, path, zarr_format = 2L),
    read = function(path) read_zarr(path),
    check = "check_zarr"
  ),
  zarr_v3 = list(
    path = file.path(cache_dir, "written_v3.zarr"),
    write = function(adata, path) write_zarr(adata, path, zarr_format = 3L),
    read = function(path) read_zarr(path),
    check = "check_zarr"
  )
)

# Python only reads the end of the arrays, as it would need too much memory
# on top of what R is using otherwise
check_written <- function() {
  reticulate::py_run_string(
    "
import h5py, zarr

def check_store(f):
    return dict(
        encoding=f['X'].attrs['encoding-type'],
        shape=[int(s) for s in f['X'].attrs['shape']],
        indptr_dtype=str(f['X/indptr'].dtype),
        indices_dtype=str(f['X/indices'].dtype),
        n=float(f['X/data'].shape[0]),
        last=float(f['X/indptr'][-1]),
        tail_data=[float(x) for x in f['X/data'][-9:]],
        tail_indices=[int(x) for x in f['X/indices'][-9:]],
    )

def check_h5ad(path):
    with h5py.File(path, 'r') as f:
        return check_store(f)

def check_zarr(path):
    return check_store(zarr.open_group(path, mode='r'))
"
  )
}

for (format in names(writers)) {
  writer <- writers[[format]]

  test_that(paste("write", format, "more than 2^31 - 1 non-zero values"), {
    skip_unless_part(paste0("write_", format))
    adata <- read_h5ad(paths$csc)
    unlink(writer$path, recursive = TRUE)
    expect_no_warning(writer$write(adata, writer$path))

    tail_k <- (per_major - 9L):(per_major - 1L)
    result <- check_written()[[writer$check]](writer$path)
    expect_identical(result$encoding, "csc_matrix")
    expect_equal(unlist(result$shape), c(n_minor, n_major))
    expect_identical(result$indptr_dtype, "int64")
    expect_identical(result$indices_dtype, "int32")
    expect_equal(result$n, nnz)
    expect_equal(result$last, nnz)
    expect_equal(
      unlist(result$tail_indices),
      3L * tail_k + (n_major - 1L) %% 3L
    )
    expect_equal(
      unlist(result$tail_data),
      ((n_major - 1L) %% 7L) + 1 + (tail_k %% 5L)
    )
  })

  test_that(paste("reread", format, "more than 2^31 - 1 non-zero values"), {
    skip_unless_part(paste0("reread_", format))
    skip_if_not(file.exists(writer$path), "Run the matching write part first")
    withr::defer(unlink(writer$path, recursive = TRUE))
    adata <- writer$read(writer$path)
    expect_s4_class(adata$X, "SVT_SparseMatrix")
    expect_large_csc(adata$X)
  })
}

test_that("convert_sce: more than 2^31 - 1 non-zero values to SCE and back", {
  skip_unless_part("convert_sce")
  skip_if_not_installed("SingleCellExperiment")

  adata <- read_h5ad(paths$csc)
  sce <- as_SingleCellExperiment(adata)
  rm(adata)
  gc(full = TRUE)

  assay <- SummarizedExperiment::assay(sce, "X")
  expect_s4_class(assay, "SVT_SparseMatrix")
  expect_identical(dim(assay), c(n_major, n_minor))
  expect_identical(SparseArray::nzcount(assay), nnz)
  rm(assay)

  back <- as_AnnData(sce, x_mapping = "X")
  expect_s4_class(back$X, "SVT_SparseMatrix")
  expect_large_csc(back$X)
})

test_that("convert_seurat: more than 2^31 - 1 non-zero values to Seurat and back", {
  skip_unless_part("convert_seurat")
  skip_if_not_installed("Seurat")

  adata <- read_h5ad(paths$csc)
  expect_message(
    seurat <- as_Seurat(adata, x_mapping = "counts"),
    "DelayedMatrix"
  )
  rm(adata)
  gc(full = TRUE)

  layer <- SeuratObject::LayerData(seurat, "counts")
  expect_s4_class(layer, "DelayedMatrix")
  expect_identical(dim(layer), c(n_major, n_minor))
  expect_equal(
    as.vector(as.matrix(layer[n_major, 1:9, drop = FALSE])),
    expected_major(n_major)
  )

  # DelayedArray only transposes when it is accessed, this adds no memory
  back <- as_AnnData(seurat, x_mapping = "counts")
  expect_s4_class(back$X, "DelayedMatrix")
  expect_identical(dim(back$X), c(n_minor, n_major))
})

# Python and R both hold the matrix in these, so Python keeps its arrays in a
# memory map to stay within the memory that R alone needs
py_large_csc <- function(dir) {
  reticulate::py_run_string(
    "
import numpy as np, h5py, os, scipy.sparse as sp

def large_csc(path, dir):
    with h5py.File(path, 'r') as f:
        out = {}
        for name, dtype in (('data', 'f8'), ('indices', 'i4')):
            file = os.path.join(dir, name + '.bin')
            n = f['X'][name].shape[0]
            if not (os.path.exists(file) and os.path.getsize(file) == n * np.dtype(dtype).itemsize):
                mm = np.memmap(file, dtype=dtype, mode='w+', shape=(n,))
                step = 2**26
                for start in range(0, n, step):
                    mm[start:start + step] = f['X'][name][start:start + step]
                mm.flush()
                del mm
            out[name] = np.memmap(file, dtype=dtype, mode='r', shape=(n,))
        indptr = f['X/indptr'][:]
        shape = tuple(int(s) for s in f['X'].attrs['shape'])
    return sp.csc_matrix((out['data'], out['indices'], indptr), shape=shape, copy=False)
",
    convert = FALSE
  )$large_csc(paths$csc, dir)
}

test_that("python_to_r: more than 2^31 - 1 non-zero values from scipy", {
  skip_unless_part("python_to_r")
  X <- py_to_r_matrix(py_large_csc(cache_dir))
  expect_s4_class(X, "SVT_SparseMatrix")
  expect_large_csc(X)
})

test_that("r_to_python: more than 2^31 - 1 non-zero values to scipy", {
  skip_unless_part("r_to_python")
  bi <- reticulate::import_builtins(convert = FALSE)

  py_x <- reticulate::r_to_py(read_h5ad(paths$csc)$X)
  expect_equal(reticulate::py_to_r(bi$float(py_x$nnz)), nnz)
  expect_equal(unlist(reticulate::py_to_r(py_x$shape)), c(n_minor, n_major))
  last <- reticulate::py_to_r(py_x$getcol(n_major - 1L)$toarray())
  expect_equal(as.vector(last)[1:9], expected_major(n_major))
})
