h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
adata <- read_h5ad(h5ad_file)
adata$obs <- adata$obs["Float"]
adata$var <- adata$var[c("String", "mean_counts")]
adata$uns <- NULL

adata <- generate_dataset(
  n_obs = 10L, 
  n_vars = 20L,
  layer_types = c("numeric_matrix", "integer_dense"),
  # obs_types = c("character", "logical"),
  var_types = c("character", "numeric"),
  obsm_types = c("numeric_matrix", "integer_dense"),
  varm_types = c("numeric_matrix", "integer_dense"),
  obsp_types = c("numeric_matrix", "integer_dense"),
  varp_types = c("numeric_matrix", "integer_dense"),
  uns_types = c("df_numeric", "scalar_integer"), 
  format = "AnnData"
)

# write_h5ad(adata, file.path(tempdir(), "dingdin.h5ad"), mode = "w")
write_h5ad(adata, tempfile(fileext = ".h5ad"), mode = "w")

# hashes <- purrr::map(6:10, \(i) {
hashes <- lapply(6:10, \(i) {
  # tmp <- tempdir()
  fp <- file.path(paste0("test5", i, ".h5ad"))
  write_h5ad(adata, fp, mode = "w")
  hash <- rlang::hash_file(fp)
  tmp2 <- file.path(paste0("test2", i, "-ls.txt"))
  res <- system(paste0("h5ls -r -S --data -v ", fp, " > ", tmp2))
  resawk <- system(paste0(r"(awk '/^\/[^\s]*/{obj=$0} /Location:/{if(obj!=""){split($2,a,":"); print a[2] "\t" obj "\t" $0; obj=""}}' )", tmp2), intern = TRUE)
  
  print(resawk)
  
  # awkr <- awk '/^\/[^\s]*/{obj=$0} /Location:/{if(obj!=""){split($2,a,":"); print a[2] "\t" obj "\t" $0; obj=""}}' empty-r1-ls.txt | sort -n | cut -f2,3
  Sys.sleep(1)
  c(hash, resawk)
})
hashes

h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
adata <- read_h5ad(h5ad_file)
adata2 <- adata
adata2$obs <- adata$obs["Float"]
adata2$var <- adata$var[c("String", "mean_counts")]
adata2$uns <- NULL
adata2

tmp <- tempdir()
write_h5ad(adata2, file.path(tmp, "test1.h5ad"), mode = "w")
Sys.sleep(5)
write_h5ad(adata2, file.path(tmp, "test2.h5ad"), mode = "w")

hash1 <- rlang::hash_file(paste0(tmp, "/test1.h5ad"))
hash2 <- rlang::hash_file(paste0(tmp, "/test2.h5ad"))

print(c(hash1, hash2))

system2("sha256sum",  paste0(tmp, "/*.h5ad"))
#> c95255dff8948e1539ee5415438604042f525bfff80099e6cff37f6fc9dd2672  /tmp/RtmpXzMfme/test1.h5ad
#> 531d8f95bcd0550fd1def7fab835ff729b8168c11d3e040246e3a411d58975f0  /tmp/RtmpXzMfme/test2.h5ad
system2("h5diff",  c(paste0(tmp, "/test1.h5ad"), paste0(tmp, "/test2.h5ad")))











library(hdf5r)

h5_path <- tempfile(fileext = ".h5")

# Writing just a dataset gives a unique hash
hashes <- purrr::map_chr(1:5, \(i) {
  dataset_properties <- hdf5r::H5P_DATASET_CREATE$new()
  dataset_properties$set_obj_track_times(track_times = FALSE)
  file <- H5File$new(h5_path, mode = "a", file_create_pl = dataset_properties)
  file[["testdataset"]] <- 1:10
  hash <- rlang::hash_file(h5_path)
  file$close_all()
  fs::file_delete(h5_path)
  
  Sys.sleep(1)
  
  hash
})

unique(hashes)
#> [1] "8a9c898ee3fdc72b868aa4e9b4d8cf56"

# Writing just an attribute gives a unique hash
hashes <- purrr::map_chr(1:5, \(i) {
  file <- H5File$new(h5_path, mode = "a")
  h5attr(file, "testattrib") <- LETTERS[1:10]
  hash <- rlang::hash_file(h5_path)
  file$close_all()
  fs::file_delete(h5_path)
  
  Sys.sleep(1)
  
  hash
})

unique(hashes)
#> [1] "b249756aeccfe9ae8982fdf696e491fc"

  file <- H5File$new(h5_path, mode = "a")
  dataset_properties <- hdf5r::H5P_DATASET_CREATE$new()
  dataset_properties$set_obj_track_times(track_times = FALSE)
  robj1 <- 1:10
  robj2 <- LETTERS[1:10]
  file$create_dataset(
    name = "testdataset",
    dims = 10L,
    dtype = NULL,
    robj = robj1,
    dataset_create_pl = dataset_properties
  )
  file$create_attr_by_name(
    attr_name = "testattrib",
    obj_name = "testdataset",
    robj = robj2,
  )


  hash <- rlang::hash_file(h5_path)
  file$close_all()
  fs::file_delete(h5_path)
  
  Sys.sleep(1)

# Writing both a dataset and an attribute gives different hashes
hashes <- purrr::map_chr(1:5, \(i) {
  file <- H5File$new(h5_path, mode = "a")
  dataset_properties <- hdf5r::H5P_DATASET_CREATE$new()
  dataset_properties$set_obj_track_times(track_times = FALSE)
  file$create_dataset(
    name = "testdataset",
    dims = 10L,
    dtype = NULL,
    robj = 1:10,
    dataset_create_pl = dataset_properties
  )
  file$create_attr_by_name(
    attr_name = "testattrib",
    obj_name = "testdataset",
    robj = LETTERS[1:10],
  )
  hash <- rlang::hash_file(h5_path)
  file$close_all()
  fs::file_delete(h5_path)
  
  Sys.sleep(1)
  
  hash
})

unique(hashes)
#> [1] "b50e1f459822fa24fffeba7a0bd3c782" "78a732d903399382dfcadffa88252912"
#> [3] "513382f4abf7a63f15d395581195114f" "16481c5fd067e3223ebc7b6a7e290013"
#> [5] "112d7700f6134854f1485215eac751dd"