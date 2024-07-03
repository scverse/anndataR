import random as rd
import pandas as pd
import numpy as np


def nullable_integer_array(n):
    nullable_array = [1 for _ in range(n)]
    # np.nan, pd.NA and None should all end up as null values, masked in the h5ad file
    nullable_array[0] = pd.NA
    nullable_array[1] = np.nan
    nullable_array[2] = None
    return pd.array(nullable_array, dtype="Int64")

def nullable_boolean_array(n):
    nullable_array = [rd.choices([True, False], k=n)]
    # np.nan, pd.NA and None should all end up as null values, masked in the h5ad file
    nullable_array[0] = pd.NA
    nullable_array[1] = np.nan
    nullable_array[2] = None
    return pd.array(nullable_array, dtype="boolean")

def missing_values_categorical(n, ordered = True):
    missing_values = pd.Categorical([["Value1", "Value2", "Value3"][i % 3] for i in range(n)], 
                              categories = ["Value1", "Value2", "Value3"],
                              ordered=ordered)
    # They should all end up as code -1 in the h5ad file
    missing_values[0] = np.nan
    missing_values[1] = pd.NA
    missing_values[2] = None
    return missing_values

vector_generators = {
    "categorical": lambda n: pd.Categorical([["Value1", "Value2", "Value3"][i % 2] for i in range(n)]),
    "categorical_ordered": lambda n: pd.Categorical([["Value1", "Value2", "Value3"][i % 2] for i in range(n)], ordered=True),
    "categorical_missing_values": lambda n: missing_values_categorical(n, ordered=False), 
    "categorical_ordered_missing_values": lambda n: missing_values_categorical(n, ordered=True),

    "string_array": lambda n: [f"value_{i}" for i in range(n)],

    # should we also check a 1d sparse array? We should probably leave it for the matrix generation?
    "dense_array": lambda n: [rd.random() for _ in range(n)],

    "integer_array": lambda n: [1 for _ in range(n)],
    "nullable_integer_array": nullable_integer_array,

    "boolean_array": lambda n: [rd.choices([True, False], k=n)],
    "nullable_boolean_array": nullable_boolean_array,
}

def generate_vector(n, vector_type):
    """
    Generate a vector of a specified type.

    Parameters:
    vector_type (str): The type of vector to generate.
    n (int): The length of the vector.

    Returns:
    list: The generated vector.

    Raises:
    AssertionError: If the vector_type is unknown.
    """
    # check if vector_type is valid
    assert vector_type in vector_generators, f"Unknown vector type: {vector_type}"

    return vector_generators[vector_type](n)
