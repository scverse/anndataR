skip_if_no_anndata()

# construct dummy objects
dummy <- dummy_data(10L, 20L)

# create anndata in python
obs_ <- dummy$obs
rownames(obs_) <- dummy$obs_names
var_ <- dummy$var
rownames(var_) <- dummy$var_names
ad <- anndata::AnnData(
  X = dummy$X,
  obs = obs_,
  var = var_
)

for (slot in c("X", "obs", "var", "obs_names", "var_names", "layers")) {
  test_label <- paste0("test Python->R read ", slot)

  test_that(test_label, {
    # write to file
    filename <- withr::local_file(paste0("read_", slot, ".h5ad"))
    ad$write_h5ad(filename)

    # read from file
    ad_new <- HDF5AnnData$new(filename)

    expect_equal(ad_new[[slot]], dummy[[slot]], tolerance = 1e-10)
  })
}



# It's not possible to do this roundtrip yet
# nolint start
# test_that("test R+HDF5AnnData -> h5ad -> reticulate+Python", {
#   # create anndata in R
#   ad <- InMemoryAnnData$new(
#     X = X,
#     obs = obs,
#     var = var,
#     obs_names = obs_names,
#     var_names = var_names
#   )

#   # convert and write to file
#   ad$to_HDF5AnnData("file.h5ad")

#   # read from file, hopefully Python anndata is able to read this
#   ad_new <- anndata::read_h5ad("file.h5ad")

#   # check whether ad_new$obs == obs and so on
#   expect_equal(ad_new$n_obs(), nrow(obs))
#   expect_equal(ad_new$n_vars(), nrow(var))
# })
# nolint end
