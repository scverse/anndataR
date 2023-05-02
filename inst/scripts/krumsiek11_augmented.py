# This script was used to create the `krumsiek11_augmented_v0-8.h5ad`
# file. It adds some extra data to the previous `krumsiek11.h5ad`
# dataset to cover some additional cases for testing (NAs, booleans,
# etc). The data was saved in AnnData=0.8.0 format.
#
# Key package versions:
#  - anndata=0.8.0
#  - h5py=3.8.0
#  - hdf5=1.14.0
#  - numpy=1.23.5
#  - pandas=1.5.3
#  - python=3.9.16
#  - scanpy=1.9.2

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

adata = ad.read_h5ad("inst/extdata/krumsiek11.h5ad")

adata.obsm["dense"] = np.add(*np.indices((adata.shape[0], 10)))
adata.obsm["sparse_csr"] = sparse.csr_matrix(np.add(*np.indices((adata.shape[0], 10))))
adata.obsm["sparse_csc"] = sparse.csc_matrix(np.add(*np.indices((adata.shape[0], 10))))
adata.obsm["sparse_csr_1"] = sparse.csr_matrix(np.add(*np.indices((adata.shape[0], 1))))
adata.obsm["sparse_csc_1"] = sparse.csc_matrix(np.add(*np.indices((adata.shape[0], 1))))

# add string column to rowData/var. Make the entries unique so it's
# saved as str instead of factor
adata.var["dummy_str"] = [f"row{i}" for i in range(adata.shape[1])]

# add float column to colData/obs
adata.obs["dummy_num"] = 42.42

# float column with NA
adata.obs["dummy_num2"] = adata.obs["dummy_num"]
adata.obs["dummy_num2"][0] = float("nan")

# int column
adata.obs["dummy_int"] = np.arange(adata.shape[0])

# int column with NA
adata.obs["dummy_int2"] = pd.array([None] + [42] * (adata.shape[0] - 1))

# bool column
adata.obs["dummy_bool"] = True
adata.obs["dummy_bool"][0] = False

# bool column with NA
adata.obs["dummy_bool2"] = pd.array([False, None] + [True] * (adata.shape[0] - 2))

# also add some entries to the metadata/uns
adata.uns["dummy_category"] = pd.array(["a", "b", None], dtype="category")

adata.uns["dummy_bool"] = [True, True, False]
adata.uns["dummy_bool2"] = pd.array([True, False, None])

adata.uns["dummy_int"] = [1,2,3]
adata.uns["dummy_int2"] = pd.array([1,2,None])

adata.uns["dummy_int_scalar"] = 1
adata.uns["dummy_string_scalar"] = "foo"

adata.uns["dummy_string_array"] = [[f"row{i}{j}" for i in range(adata.shape[1])] for j in range(adata.shape[0])]

# for testing obs with no content
adata.obs = pd.DataFrame(index = adata.obs.index) # empty pd.DataFrame

adata.write("inst/extdata/krumsiek11_augmented_sparse_empty_v0-8.h5ad")
