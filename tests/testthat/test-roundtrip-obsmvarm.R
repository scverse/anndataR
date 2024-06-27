# TODO: re-enable
# nolint start
# skip_if_no_anndata()
# skip_if_not_installed("rhdf5")

# data <- generate_dataset(10L, 20L)

# obsm_names <- names(data$obsm)
# # TODO: Add denseMatrix support to anndata and anndataR
# obsm_names <- obsm_names[!grepl("_dense", obsm_names)]

# for (name in obsm_names) {
#   test_that(paste0("roundtrip with obsm and varm '", name, "'"), {
#     # create anndata
#     ad <- AnnData(
#       obsm = data$obsm[name],
#       varm = data$varm[name],
#       obs_names = data$obs_names,
#       var_names = data$var_names
#     )

#     # write to file
#     filename <- withr::local_file(paste0("roundtrip_obsmvarm_", name, ".h5ad"))
#     write_h5ad(ad, filename)

#     # read from file
#     ad_new <- read_h5ad(filename, to = "HDF5AnnData")

#     # expect slots are unchanged
#     expect_equal(
#       ad_new$obsm[[name]],
#       data$obsm[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#     expect_equal(
#       ad_new$varm[[name]],
#       data$varm[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#   })
# }

# r2py_names <- names(data$obsm)

# # TODO: remove this when https://github.com/scverse/anndata/issues/1146 is fixed
# r2py_names <- r2py_names[!grepl("_with_nas", r2py_names)]

# for (name in r2py_names) {
#   test_that(paste0("reticulate->hdf5 with obsm and varm '", name, "'"), {
#     # add rownames
#     obsm <- data$obsm[name]
#     varm <- data$varm[name]
#     rownames(obsm[[name]]) <- data$obs_names
#     rownames(varm[[name]]) <- data$var_names

#     # create anndata
#     ad <- anndata::AnnData(
#       obsm = obsm,
#       varm = varm,
#       shape = dim(data$X),
#       obs = data.frame(row.names = data$obs_names),
#       var = data.frame(row.names = data$var_names)
#     )

#     # write to file
#     filename <- withr::local_file(paste0("reticulate_to_hdf5_obsmvarm_", name, ".h5ad"))
#     ad$write_h5ad(filename)

#     # read from file
#     ad_new <- HDF5AnnData$new(filename)

#     # expect slots are unchanged
#     expect_equal(
#       ad_new$obsm[[name]],
#       data$obsm[[name]],
#       tolerance = 1e-6
#     )
#     expect_equal(
#       ad_new$varm[[name]],
#       data$varm[[name]],
#       tolerance = 1e-6
#     )
#   })
# }

# for (name in r2py_names) {
#   test_that(paste0("hdf5->reticulate with obsm and varm '", name, "'"), {
#     # write to file
#     filename <- withr::local_file(paste0("hdf5_to_reticulate_obsmvarm_", name, ".h5ad"))

#     # strip rownames
#     obsm <- data$obsm[name]
#     varm <- data$varm[name]
#     rownames(obsm[[name]]) <- NULL
#     rownames(varm[[name]]) <- NULL

#     # make anndata
#     ad <- HDF5AnnData$new(
#       file = filename,
#       obsm = obsm,
#       varm = varm,
#       obs_names = data$obs_names,
#       var_names = data$var_names
#     )

#     # read from file
#     ad_new <- anndata::read_h5ad(filename)

#     # expect slots are unchanged
#     obsm_ <- ad_new$obsm[[name]]
#     if (!is.null(obsm_)) {
#       rownames(obsm_) <- NULL
#       colnames(obsm_) <- NULL
#     }
#     expect_equal(
#       obsm_,
#       data$obsm[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#     varm_ <- ad_new$varm[[name]]
#     if (!is.null(varm_)) {
#       rownames(varm_) <- NULL
#       colnames(varm_) <- NULL
#     }
#     expect_equal(
#       varm_,
#       data$varm[[name]],
#       ignore_attr = TRUE,
#       tolerance = 1e-6
#     )
#   })
# }
# nolint end
