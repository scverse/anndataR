# /// script
# requires-python = "==3.14.4"
# dependencies = [
#   "anndata==0.12.10",
#   "igraph==1.0.0",
#   "leidenalg==0.11.0",
#   "scanpy==1.12.1",
#   "scipy==1.17.1",
#   "zarr==3.1.6",
# ]
# ///
import os
import shutil
import zipfile

import anndata
import numpy
import pandas
import scanpy
import scipy.sparse

# This script uses Python to create example H5AD and Zarr files for testing
# interoperability between languages. It is designed to be a small but
# relatively complex file that tests reading of different types and data
# structures.
#
# In order to run the script, install uv (https://docs.astral.sh/uv/) and run:
#
# uv run inst/scripts/example_files.py
#
# The standard scanpy workflow has also been applied to populate
# some of the most common information from real analyses. It should be updated
# to test new issues as they are discovered.
#
# NOTE: When updating this script for the {anndataR} example H5AD file please
# update the package versions used above, update the script version, date and
# changelog below and format the file using Ruff (https://docs.astral.sh/ruff/):
#
# ruff format inst/scripts/example_files.py && ruff check --select I --fix inst/scripts/example_files.py
#
# Version: 0.4.1
# Date: 2026-04-15
#
# CHANGELOG
#
# v0.4.1 (2026-04-15)
# - Replace requirements.yml with uv dependency comments
# - Update package versions to latest stable versions
# - Add progress messages
# - Use Ruff for formatting
# v0.4.0 (2025-11-24)
# - Add zarr example
# - Add requirements.yml
# v0.3.0 (2025-08-04)
# - Add adata.varp["test_varp"] to test reading of varp
# - Update package versions to latest stable versions
# v0.3.0 (2025-08-04)
# - Add adata.varp["test_varp"] to test reading of varp
# - Update package versions to latest stable versions
# v0.2.0 (2023-05-11)
# - Add 1D sparse matrix to `adata.uns["Sparse1D"]
# - Reduce the size of `adata.uns["String2D"]` and add columns to values
# v0.1.1 (2023-05-09)
# - Reduce the size of `adata.uns["String2D"]` to save space
# - Reduce dimension to 50 x 100 to save space
# v0.1.0 (2023-05-08)
# - Initial version

numpy.random.seed(0)

print(">>> Creating AnnData...")

# Randomly generate a counts matrix
counts = numpy.random.poisson(2, size=(50, 100))

# Create an AnnData
adata = anndata.AnnData(scipy.sparse.csr_matrix(counts.copy(), dtype=numpy.float32))
adata.obs_names = [f"Cell{i:03d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene{i:03d}" for i in range(adata.n_vars)]

# Populate layers with different matrix types
adata.layers["counts"] = adata.X.copy()
adata.layers["dense_counts"] = counts.copy()
adata.layers["csc_counts"] = scipy.sparse.csc_matrix(counts.copy(), dtype=numpy.float32)

# Populate adata.var with different types
adata.var["String"] = [f"String{i}" for i in range(adata.n_vars)]

# Populate adata.obs with different types
adata.obs["Float"] = 42.42
adata.obs["FloatNA"] = adata.obs["Float"]
adata.obs["FloatNA"][0] = float("nan")
adata.obs["Int"] = numpy.arange(adata.n_obs)
adata.obs["IntNA"] = pandas.array([None] + [42] * (adata.n_obs - 1))
adata.obs["Bool"] = pandas.array([False] + [True] * (adata.n_obs - 1))
adata.obs["BoolNA"] = pandas.array([False, None] + [True] * (adata.n_obs - 2))

# Populate adata.uns with different types
adata.uns["Category"] = pandas.array(["a", "b", None], dtype="category")
adata.uns["Bool"] = [True, True, False]
adata.uns["BoolNA"] = pandas.array([True, False, None])
adata.uns["Int"] = [1, 2, 3]
adata.uns["IntNA"] = pandas.array([1, 2, None])
adata.uns["IntScalar"] = 1
adata.uns["Sparse1D"] = scipy.sparse.csc_matrix([1, 2, 0, 0, 0, 3])
adata.uns["StringScalar"] = "A string"
adata.uns["String"] = [f"String {i}" for i in range(10)]
adata.uns["String2D"] = [[f"row{i}col{j}" for i in range(10)] for j in range(5)]
adata.uns["DataFrameEmpty"] = pandas.DataFrame(index=adata.obs.index)

print("\n>>> Running scanpy workflow...")

# Run the standard scanpy workflow
print("Calculating QC metrics...")
scanpy.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
print("Normalizing..")
scanpy.pp.normalize_total(adata, inplace=True)
adata.layers["dense_X"] = adata.X.copy().toarray()
scanpy.pp.log1p(adata)
print("Finding highly variable genes...")
scanpy.pp.highly_variable_genes(adata)
print("Calculating PCA...")
scanpy.tl.pca(adata)
print("Finding neighbors...")
scanpy.pp.neighbors(adata)
print("Calculating UMAP...")
scanpy.tl.umap(adata)
print("Calculating Leiden clusters...")
scanpy.tl.leiden(adata)
print("Calculating marker genes...")
scanpy.tl.rank_genes_groups(adata, "leiden")

# add varp to test reading of varp
adata.varp["test_varp"] = numpy.random.rand(adata.n_vars, adata.n_vars)

print("\n>>> Writing H5AD file...")
adata.write_h5ad("inst/extdata/example.h5ad", compression="gzip")

# Write Zarr files in both v2 and v3 formats and zip them
os.chdir("inst/extdata/")
for fmt in (2, 3):
    anndata.settings.zarr_write_format = fmt

    zarr_dir = f"example_v{fmt}.zarr"
    zip_path = f"{zarr_dir}.zip"

    print(f"\n>>> Writing Zarr v{fmt} file...")
    adata.write_zarr(zarr_dir)

    print(f"Zipping Zarr v{fmt} file...")

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for root, dirs, files in os.walk(zarr_dir):
            for file in files:
                z.write(os.path.join(root, file))

    shutil.rmtree(zarr_dir)

print("\n>>> Done!")
