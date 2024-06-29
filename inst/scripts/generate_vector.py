import random as rd
import pandas as pd
import numpy as np

# TODO: I don't like the randomness here, find a different way
# TODO: check what numpy & pandas does with np.nan and pd.NA and None? Should these be included as well?
# TODO: check as what these values are written to the h5ad files
vector_generators = {
    "categorical": lambda n: pd.array(rd.choices(["Value1", "Value2", "Value3"], k=n), dtype="category"),
    "categorical_ordered": lambda n: pd.array(rd.choices(["Value1", "Value2", "Value3"], k=n), 
                                              dtype=pd.CategoricalDtype(["Value1", "Value2", "Value3"], ordered=True)),
    "categorical_missing_values": lambda n: pd.array(rd.choices(["Value1", "Value2", "Value3", np.nan], k=n), dtype="category"), 
    "categorical_ordered_missing_values": lambda n: pd.array(rd.choices(["Value1", "Value2", "Value3", np.nan], k=n), 
                                                             dtype=pd.CategoricalDtype(["Value1", "Value2", "Value3"], ordered=True)),
    
    "string_array": lambda n: [f"value" for _ in range(n)],

    # should we also check a 1d sparse array? We should probably leave it for the matrix generation
    "dense_array": lambda n: [rd.random() for _ in range(n)],

    "integer_array": lambda n: [1 for _ in range(n)],
    "nullable_integer_array": lambda n: pd.array(rd.choices(["Value1", "Value2", np.nan], k=n), dtype="Int64"),

    "boolean_array": lambda n: [rd.choices([True, False], k=n)],
    # should we also test np.nan in a boolean array (which should always evaluate to False)?
    "nullable_boolean_array": lambda n: pd.array(rd.choices(["Value1", "Value2", pd.NA], k=n), dtype="boolean"),
    "boolean_array_nan": lambda n: pd.array(rd.choices([True, False, np.nan], k=n), dtype="boolean"),

}

def generate_vector(vector_type, n):
    # check if vector_type is valid
    if vector_type not in vector_generators:
        raise ValueError(f"Invalid vector type: {vector_type}")
    return vector_generators[vector_type](n)
