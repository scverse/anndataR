to_InMemory <- function(ad) {
  InMemoryAnnData$new(
    X = ad$X,
    obs = ad$obs,
    var = ad$var,
    obs_names = ad$obs_names,
    var_names = ad$var_names,
    layers = ad$layers
  )
}