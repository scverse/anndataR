#' create_zarr_group
#'
#' create zarr groups
#'
#' @param store the location of (zarr) store
#' @param name name of the group
#' @param version zarr version
#' @export
create_zarr_group <- function(store, name, version = "v2"){
  split.name <- strsplit(name, split = "\\/")[[1]]
  if(length(split.name) > 1){
    split.name <- vapply(seq_len(length(split.name)),
                         function(x) paste(split.name[seq_len(x)], collapse = "/"),
                         FUN.VALUE = character(1))
    split.name <- rev(tail(split.name,2))
    if(!dir.exists(file.path(store,split.name[2])))
      create_zarr_group(store = store, name = split.name[2])
  }
  dir.create(file.path(store, split.name[1]), showWarnings = FALSE)
  switch(version,
         v2 = {
           write("{\"zarr_format\":2}", file = file.path(store, split.name[1], ".zgroup"))},
         v3 = {
           stop("Currently only zarr v2 is supported!")
         },
         stop("only zarr v2 is supported. Use version = 'v2'")
  )
}

#' create_zarr
#'
#' create zarr store
#'
#' @param dir the location of zarr store
#' @param prefix prefix of the zarr store
#' @param version zarr version
#'
#' @examples
#' dir.create(td <- tempfile())
#' zarr_name <- "test"
#' create_zarr(dir = td, prefix = "test")
#' dir.exists(file.path(td, "test.zarr"))
#'
#' @export
create_zarr <- function(dir, prefix, version = "v2"){
  create_zarr_group(store = dir, name = paste0(prefix, ".zarr"), version = version)
}

#' Read the .zattrs file associated with a Zarr array or group
#'
#' @param path A character vector of length 1. This provides the
#'   path to a Zarr array or group. This can either be on a local file
#'   system or on S3 storage.
#' @param s3_client A list representing an S3 client.  This should be produced
#' by [paws.storage::s3()].
#'
#' @returns A list containing the .zattrs elements
#'
#' @importFrom jsonlite read_json fromJSON
#' @importFrom stringr str_extract str_remove
#'
#' @export
read_zattrs <- function(path, s3_client = NULL) {
  path <- .normalize_array_path(path)
  zattrs_path <- paste0(path, ".zattrs")

  if(!file.exists(zattrs_path))
    stop("The group or array does not contain attributes (.zattrs)")

  if (!is.null(s3_client)) {

    parsed_url <- parse_s3_path(zattrs_path)

    s3_object <- s3_client$get_object(Bucket = parsed_url$bucket,
                                      Key = parsed_url$object)

    zattrs <- fromJSON(rawToChar(s3_object$Body))
  } else {
    zattrs <- read_json(zattrs_path)
  }
  return(zattrs)
}

#' Read the .zattrs file associated with a Zarr array or group
#'
#' @param path A character vector of length 1. This provides the
#'   path to a Zarr array or group.
#' @param new.zattrs a list inserted to .zattrs at the \code{path}.
#' @param overwrite if TRUE, existing .zattrs elements will be overwritten by \code{new.zattrs}.
#'
#' @importFrom jsonlite toJSON
#'
#' @export
write_zattrs <- function(path, new.zattrs = list(), overwrite = TRUE){
  path <- .normalize_array_path(path)
  zattrs_path <- paste0(path, ".zattrs")

  if(is.null(names(new.zattrs)))
    stop("list elements should be named")

  if("" %in% names(new.zattrs)){
    message("Ignoring unnamed list elements")
    new.zattrs <- new.zattrs[which(names(new.zattrs == ""))]
  }

  if(file.exists(zattrs_path)){
    old.zattrs <- read_json(zattrs_path)
    if(overwrite){
      old.zattrs <- old.zattrs[setdiff(names(old.zattrs), names(new.zattrs))]
    } else {
      new.zattrs <- new.zattrs[setdiff(names(new.zattrs), names(old.zattrs))]
    }
    new.zattrs <- c(old.zattrs, new.zattrs)
  }

  json <- .format_json(toJSON(new.zattrs, auto_unbox = TRUE, pretty = TRUE, null = "null"))
  write(x = json, file = zattrs_path)
}

#' Normalize a Zarr array path
#'
#' Taken from https://zarr.readthedocs.io/en/stable/spec/v2.html#logical-storage-paths
#'
#' @param path Character vector of length 1 giving the path to be normalised.
#'
#' @returns A character vector of length 1 containing the normalised path.
#'
#' @keywords Internal
.normalize_array_path <- function(path) {
  ## we strip the protocol because it gets messed up by the slash removal later
  if (grepl(x = path, pattern = "^((https?://)|(s3://)).*$")) {
    root <- gsub(x = path, pattern = "^((https?://)|(s3://)).*$",
                 replacement = "\\1")
    path <- gsub(x = path, pattern = "^((https?://)|(s3://))(.*$)",
                 replacement = "\\4")
  } else {
    ## Replace all backward slash ("\\") with forward slash ("/")
    path <- gsub(x = path, pattern = "\\", replacement = "/", fixed = TRUE)
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    root <- gsub(x = path, "(^[[:alnum:]:.]*/)(.*)", replacement = "\\1")
    path <- gsub(x = path, "(^[[:alnum:]:.]*/)(.*)", replacement = "\\2")
  }

  ## Strip any leading "/" characters
  path <- gsub(x = path, pattern = "^/", replacement = "", fixed = FALSE)
  ## Strip any trailing "/" characters
  path <- gsub(x = path, pattern = "/$", replacement = "", fixed = FALSE)
  ## Collapse any sequence of more than one "/" character into a single "/"
  path <- gsub(x = path, pattern = "//*", replacement = "/", fixed = FALSE)
  ## The key prefix is then obtained by appending a single "/" character to
  ## the normalized logical path.
  path <- paste0(root, path, "/")

  return(path)
}

.format_json <- function(json) {
  json <- gsub(x = json, pattern = "[", replacement = "[\n    ", fixed = TRUE)
  json <- gsub(x = json, pattern = "],", replacement = "\n  ],", fixed = TRUE)
  json <- gsub(x = json, pattern = ", ", replacement = ",\n    ", fixed = TRUE)
  return(json)
}
