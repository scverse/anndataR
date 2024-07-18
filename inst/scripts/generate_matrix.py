import numpy as np
import scipy as sp


def float_mtx(n_obs, n_vars, NAs=False):
    mtx = np.random.random((n_obs, n_vars))
    if NAs:  # numpy matrices do no support pd.NA
        mtx[0, 0] = np.nan
        mtx[0, 1] = None  # gets converted to np.nan
    return mtx


def int_mtx(n_obs, n_vars):
    mtx = np.random.randint(0, 101, size=(n_obs, n_vars))
    return mtx


# Possible matrix generators
# integer matrices do not support NAs in Python
matrix_generators = {
    "generate_float_matrix": lambda n_obs, n_vars: float_mtx(n_obs, n_vars),
    "generate_float_matrix_nas": lambda n_obs, n_vars: float_mtx(
        n_obs, n_vars, NAs=True
    ),
    "generate_float_csparse": lambda n_obs, n_vars: sp.sparse.csc_matrix(
        float_mtx(n_obs, n_vars)
    ),
    "generate_float_csparse_nas": lambda n_obs, n_vars: sp.sparse.csc_matrix(
        float_mtx(n_obs, n_vars, NAs=True)
    ),
    "generate_float_rsparse": lambda n_obs, n_vars: sp.sparse.csr_matrix(
        float_mtx(n_obs, n_vars)
    ),
    "generate_float_rsparse_nas": lambda n_obs, n_vars: sp.sparse.csr_matrix(
        float_mtx(n_obs, n_vars, NAs=True)
    ),
    "generate_integer_matrix": lambda n_obs, n_vars: int_mtx(n_obs, n_vars),
    "generate_integer_csparse": lambda n_obs, n_vars: sp.sparse.csc_matrix(
        int_mtx(n_obs, n_vars)
    ),
    "generate_integer_rsparse": lambda n_obs, n_vars: sp.sparse.csr_matrix(
        int_mtx(n_obs, n_vars)
    ),
}


def generate_matrix(n_obs, n_vars, matrix_type):
    """
    Generate a matrix of given dimensions and type.

    Parameters:
        n_obs (int): The number of observations (rows) in the matrix.
        n_vars (int): The number of variables (columns) in the matrix.
        matrix_type (str): The type of matrix to generate.

    Returns:
        The generated matrix, either numpy.ndarray or scipy.sparse.csc_matrix or scipy.sparse.csr_matrix.

    Raises:
        AssertionError: If the matrix_type is unknown.

    """
    assert matrix_type in matrix_generators, f"Unknown matrix type: {matrix_type}"

    return matrix_generators[matrix_type](n_obs, n_vars)
