import pandas as pd
from generate_vector import vector_generators

def generate_dataframe(n_rows, types = None):
    if types is None:
        types = list(vector_generators.keys())

    data = {t: vector_generators[t](n_rows) for t in types}
    return pd.DataFrame(data)
