import anndata as ad

from generate_matrix import matrix_generators
from generate_vector import vector_generators
from generate_dataframe import generate_dataframe
from generate_dict import scalar_generators, generate_dict


def generate_dataset(
    n_obs=10,
    n_vars=20,
    x_type="generate_integer_matrix",
    layer_types=None,
    obs_types=None,
    var_types=None,
    obsm_types=None,
    varm_types=None,
    obsp_types=None,
    varp_types=None,
    uns_types=None,
):

    assert x_type in matrix_generators, f"Unknown matrix type: {x_type}"
    assert layer_types is None or all(
        t in matrix_generators.keys() for t in layer_types
    ), "Unknown layer type"
    assert obs_types is None or all(
        t in vector_generators.keys() for t in obs_types
    ), "Unknown obs type"
    assert var_types is None or all(
        t in vector_generators.keys() for t in var_types
    ), "Unknown var type"
    assert obsm_types is None or all(
        t in matrix_generators.keys() or t in vector_generators.keys()
        for t in obsm_types
    ), "Unknown obsm type"
    assert varm_types is None or all(
        t in matrix_generators.keys() or t in vector_generators.keys()
        for t in varm_types
    ), "Unknown varm type"
    assert obsp_types is None or all(
        t in matrix_generators.keys() for t in obsp_types
    ), "Unknown obsp type"
    assert varp_types is None or all(
        t in matrix_generators.keys() for t in varp_types
    ), "Unknown varp type"
    # TODO uns types

    if layer_types is None:  # layer_types are all matrices
        layer_types = list(matrix_generators.keys())
    if obs_types is None:  # obs_types are all vectors
        obs_types = list(vector_generators.keys())
    if var_types is None:  # var_types are all vectors
        var_types = list(vector_generators.keys())
    if obsm_types is None:  # obsm_types are all matrices or vectors
        obsm_types = list(matrix_generators.keys()) + list(vector_generators.keys())
    if varm_types is None:  # varm_types are all matrices or vectors
        varm_types = list(matrix_generators.keys()) + list(vector_generators.keys())
    if obsp_types is None:  # obsp_types are all matrices
        obsp_types = list(matrix_generators.keys())
    if varp_types is None:  # varp_types are all matrices
        varp_types = list(matrix_generators.keys())
    if uns_types is None:
        uns_types = (
            list(vector_generators.keys())
            + list(matrix_generators.keys())
            + list(scalar_generators.keys())
        )

    X = matrix_generators[x_type](n_obs, n_vars)
    layers = {t: matrix_generators[t](n_obs, n_vars) for t in layer_types}

    obs_names = [f"Cell{i:03d}" for i in range(n_obs)]
    var_names = [f"Gene{i:03d}" for i in range(n_vars)]

    obs = generate_dataframe(n_obs, obs_types)
    var = generate_dataframe(n_vars, var_types)
    obs.index = obs_names
    var.index = var_names

    obsm = {}
    for t in obsm_types:
        if t in matrix_generators.keys():
            obsm[t] = matrix_generators[t](n_obs, n_obs)
        elif t in vector_generators.keys():
            obsm[t] = vector_generators[t](n_obs)

    varm = {}
    for t in varm_types:
        if t in matrix_generators.keys():
            varm[t] = matrix_generators[t](n_vars, n_vars)
        elif t in vector_generators.keys():
            varm[t] = vector_generators[t](n_vars)

    obsp = {t: matrix_generators[t](n_obs, n_obs) for t in obsp_types}
    varp = {t: matrix_generators[t](n_vars, n_vars) for t in varp_types}

    uns = generate_dict(n_obs, n_vars, uns_types)

    #TODO figure out why it cannot write NAs in a GROUP

    return ad.AnnData(
        X,
        layers=layers,
        obs=obs,
        var=var,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        uns=uns,
    )
