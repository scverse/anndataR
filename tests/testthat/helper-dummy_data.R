dummy_data <- function(n_obs = 10L, n_vars = 20L) {
  # generate X
  X <- Matrix::rsparsematrix(nrow = n_obs, ncol = n_vars, density = .1)

  # generate layers
  layers <- list(
    X2 = X * 2,
    X3 = X * 3
  )

  # generate obs
  obs <- data.frame(
    cell_type = sample(c("tcell", "bcell"), n_obs, replace = TRUE),
    cluster = sample.int(3, n_obs, replace = TRUE)
  )

  # generate var
  var <- data.frame(
    geneinfo = sample(c("a", "b", "c"), n_vars, replace = TRUE)
  )

  # generate obs_names
  obs_names <- paste0("cell", seq_len(n_obs))

  # generate var_names
  var_names <- paste0("gene", seq_len(n_vars))

  list(
    X = X,
    obs = obs,
    obs_names = obs_names,
    var = var,
    var_names = var_names,
    layers = layers
  )
}