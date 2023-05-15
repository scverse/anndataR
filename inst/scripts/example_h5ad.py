# python v3.10.10
import anndata  # anndata v0.8.0
import scanpy  # scanpy v1.9.3
import numpy  # numpy v1.23.5
import pandas  # pandas v2.0.0
import scipy.sparse  # scipy v1.10.1

# This script uses Python to create an example H5AD file for testing
# interoperability between languages. It is designed to be a small but
# relatively complex file that tests reading of different types and data
# structures. The standard scanpy workflow has also been applied to populate
# some of the most common information from real analyses. It should be updated
# to test new issues as they are discovered.
#
# NOTE: When updating this script for the {anndataR} example H5AD file please
# update the package versions used above, update the script version, date and
# changelog below and format the file using Python Black
# (https://black.readthedocs.io/en/stable/).
#
# Version: 0.2.0
# Date: 2023-05-11
#
# CHANGELOG
#
# v0.2.0 (2023-05-11)
# - Add 1D sparse matrix to `adata.uns["Sparse1D"]
# - Reduce the size of `adata.uns["String2D"]` and add columns to values
# v0.1.1 (2023-05-09)
# - Reduce the size of `adata.uns["String2D"]` to save space
# - Reduce dimension to 50 x 100 to save space
# v0.1.0 (2023-05-08)
# - Initial version

numpy.random.seed(0)

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

# Run the standard scanpy workflow
scanpy.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
scanpy.pp.normalize_total(adata, inplace=True)
adata.layers["dense_X"] = adata.X.copy().toarray()
scanpy.pp.log1p(adata)
scanpy.pp.highly_variable_genes(adata)
scanpy.tl.pca(adata)
scanpy.pp.neighbors(adata)
scanpy.tl.umap(adata)
scanpy.tl.leiden(adata)
scanpy.tl.rank_genes_groups(adata, "leiden")

# Write the H5AD file
adata.write("example.h5ad")
