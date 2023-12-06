# TODO: re-enable
# nolint start
# skip_if_no_anndata()
# skip_if_not_installed("rhdf5")

# data <- generate_dataset(10L, 20L)

# for (name in names(data$obsp)) {
#   test_that(paste0("roundtrip with obsp and varp '", name, "'"), {
#     # create anndata
#     ad <- AnnData(
#       obsp = data$obsp[name],
#       varp = data$varp[name],
#       obs_names = data$obs_names,
#       var_names = data$var_names
#     )

#     # write to file
#     filename <- withr::local_file(paste0("roundtrip_obspvarp_", name, ".h5ad"))
#     write_h5ad(ad, filename)

#     # read from file
#     ad_new <- read_h5ad(filename, to = "HDF5AnnData")

#     # expect slots are unchanged
#     expect_equal(
#       ad_new$obsp[[name]],
#       data$obsp[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#     expect_equal(
#       ad_new$varp[[name]],
#       data$varp[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#   })
# }

# for (name in names(data$obsp)) {
#   test_that(paste0("reticulate->hdf5 with obsp and varp '", name, "'"), {
#     # add rownames
#     obsp <- data$obsp[name]
#     varp <- data$varp[name]
#     rownames(obsp[[name]]) <- colnames(obsp[[name]]) <- data$obs_names
#     rownames(varp[[name]]) <- colnames(varp[[name]]) <- data$var_names

#     #  create anndata
#     ad <- anndata::AnnData(
#       obsp = obsp,
#       varp = varp,
#       shape = dim(data$X),
#       obs = data.frame(row.names = data$obs_names),
#       var = data.frame(row.names = data$var_names)
#     )

#     # write to file
#     filename <- withr::local_file(paste0("reticulate_to_hdf5_obspvarp_", name, ".h5ad"))
#     ad$write_h5ad(filename)

#     # read from file
#     ad_new <- HDF5AnnData$new(filename)

#     # expect slots are unchanged
#     expect_equal(
#       ad_new$obsp[[name]],
#       data$obsp[[name]],
#       tolerance = 1e-6
#     )
#     expect_equal(
#       ad_new$varp[[name]],
#       data$varp[[name]],
#       tolerance = 1e-6
#     )
#   })
# }

# for (name in names(data$obsp)) {
#   test_that(paste0("hdf5->reticulate with obsp and varp '", name, "'"), {
#     # write to file
#     filename <- withr::local_file(paste0("hdf5_to_reticulate_obspvarp_", name, ".h5ad"))

#     # strip rownames
#     obsp <- data$obsp[name]
#     varp <- data$varp[name]
#     rownames(obsp[[name]]) <- NULL
#     rownames(varp[[name]]) <- NULL

#     # make anndata
#     ad <- HDF5AnnData$new(
#       file = filename,
#       obsp = obsp,
#       varp = varp,
#       obs_names = data$obs_names,
#       var_names = data$var_names
#     )

#     # read from file
#     ad_new <- anndata::read_h5ad(filename)

#     # expect slots are unchanged
#     obsp_ <- ad_new$obsp[[name]]
#     if (!is.null(obsp_)) {
#       rownames(obsp_) <- colnames(obsp_) <- NULL
#     }
#     expect_equal(
#       obsp_,
#       data$obsp[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#     varp_ <- ad_new$varp[[name]]
#     if (!is.null(varp_)) {
#       rownames(varp_) <- colnames(varp_) <- NULL
#     }
#     expect_equal(
#       varp_,
#       data$varp[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#   })
# }
# nolint end
