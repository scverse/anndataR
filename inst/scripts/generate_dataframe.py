import pandas as pd
from generate_vector import vector_generators


def generate_dataframe(n_rows, types=None):
    """
    Generate a pandas DataFrame with specified number of rows and column types.

    Parameters:
        n_rows (int): The number of rows in the DataFrame.
        types (list, optional): A list of column types to include in the DataFrame.
                                Choose from the list of vector_generators keys.
                                If not provided, all available column types will be included.

    Returns:
        pandas.DataFrame: The generated DataFrame.

    """
    if types is None:
        types = list(vector_generators.keys())

    data = {t: vector_generators[t](n_rows) for t in types}
    return pd.DataFrame(data)
