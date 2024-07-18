from generate_vector import vector_generators
from generate_matrix import matrix_generators

import pandas as pd
import numpy as np

scalar_generators = {
    "string": "version",
    "char": "a",
    "integer": 1,
    "float": 1.0,
    "boolean": True,
    "none": None,
    "NA": pd.NA,
    "nan": np.nan,
}


def generate_scalar(scalar_type):
    if scalar_type[:7] == "scalar_":
        return vector_generators[scalar_type[7:]](1)
    return scalar_generators[scalar_type]


def generate_type(type, n_rows, n_cols):
    if type in scalar_generators or type[:7] == "scalar_":
        return generate_scalar(type)
    if type in vector_generators:
        return vector_generators[type](n_rows)
    if type in matrix_generators:
        return matrix_generators[type](n_rows, n_cols)
    return None


def generate_dict(n_rows, n_cols, types=None, nested=True):
    if types is None:  # types are all vectors and all matrices
        scalar_types = list(scalar_generators.keys()) + [
            f"scalar_{t}" for t in vector_generators.keys()
        ]
        types = (
            scalar_types
            + list(vector_generators.keys())
            + list(matrix_generators.keys())
        )

    data = {t: generate_type(t, n_rows, n_cols) for t in types}
    if nested:
        data["nested"] = generate_dict(n_rows, n_cols, types, False)

    return data
