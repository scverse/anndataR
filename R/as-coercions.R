#' Coercion helpers for `as()`
#'
#' These helper registrations wire up S4-style `as()` conversions so that
#' AnnData implementations (including [`InMemoryAnnData`], [`HDF5AnnData`], and
#' [`ReticulateAnnData`]) as well as [`SingleCellExperiment`] and
#' [`SeuratObject::Seurat`] objects can be coerced
#' between one another without the caller needing to know the underlying helper
#' functions. Because `as()` cannot accept additional arguments, conversions
#' that require them (such as writing HDF5-backed AnnData objects) raise an
#' informative error pointing to the richer interface.
#'
#' @keywords internal
#' @noRd
NULL

# Class compatibility registrations -----------------------------------------

.register_oldclass <- function(class, super = character()) {
  if (!methods::isClass(class)) {
    methods::setOldClass(c(class, super))
  }
}

.register_oldclass("AbstractAnnData", "R6")
.register_oldclass("InMemoryAnnData", c("AbstractAnnData", "R6"))
.register_oldclass("HDF5AnnData", c("AbstractAnnData", "R6"))
.register_oldclass("ReticulateAnnData", c("AbstractAnnData", "R6"))
.register_oldclass("AnnDataView", c("AbstractAnnData", "R6"))

.as_abort_extra_args <- function(from, to, helper) {
  cli::cli_abort(
    c(
      "Can't coerce {.cls {from}} to {.cls {to}} with {.fun as()}.",
      "i" = helper
    ),
    call = rlang::caller_env()
  )
}

.warn_as_limited <- function(recommendation) {
  cli::cli_warn(
    c(
      "Using {.fun as} limits control over data mapping.",
      "i" = recommendation
    ),
    call = rlang::caller_env()
  )
}

# Handler constructors -------------------------------------------------------

.make_convert_handler <- function(converter, warn = NULL, pre = NULL) {
  force(warn)
  force(pre)

  function(from) {
    if (!is.null(pre)) {
      pre(from)
    }
    if (!is.null(warn)) {
      .warn_as_limited(warn)
    }
    converter(from)
  }
}

.make_abort_handler <- function(from_class, to_class, helper) {
  force(from_class)
  force(to_class)
  force(helper)

  function(from) {
    .as_abort_extra_args(from_class, to_class, helper)
  }
}

.register_set_as_rules <- function(rules) {
  for (rule in rules) {
    methods::setAs(rule$from, rule$to, rule$handler)
  }
}

.format_control_recommendation <- function(call_expr) {
  sprintf("Prefer `%s` for fine-grained control over data mapping.", call_expr)
}

# AnnData <-> AnnData coercion rules ----------------------------------------
warn_ann_inmemory <- .format_control_recommendation(
  "adata$as_InMemoryAnnData(...)"
)
warn_ann_reticulate <- .format_control_recommendation(
  "adata$as_ReticulateAnnData(...)"
)

ann_data_rules <- list(
  list(
    from = "AbstractAnnData",
    to = "InMemoryAnnData",
    handler = .make_convert_handler(
      converter = as_InMemoryAnnData,
      warn = warn_ann_inmemory
    )
  ),
  list(
    from = "AbstractAnnData",
    to = "ReticulateAnnData",
    handler = .make_convert_handler(
      converter = as_ReticulateAnnData,
      warn = warn_ann_reticulate
    )
  ),
  list(
    from = "AbstractAnnData",
    to = "HDF5AnnData",
    handler = .make_abort_handler(
      from_class = "AbstractAnnData",
      to_class = "HDF5AnnData",
      helper = "Use `adata$as_HDF5AnnData(file = <path>)` to provide the output file."
    )
  )
)

.register_set_as_rules(ann_data_rules)

# SingleCellExperiment coercion rules ---------------------------------------

if (rlang::is_installed("SingleCellExperiment")) {
  warn_customise <- .format_control_recommendation("as_AnnData(...)")
  warn_sce <- .format_control_recommendation(
    "adata$as_SingleCellExperiment(...)"
  )

  single_cell_rules <- list(
    list(
      from = "SingleCellExperiment",
      to = "InMemoryAnnData",
      handler = .make_convert_handler(
        converter = function(from) {
          as_AnnData(from, output_class = "InMemoryAnnData")
        },
        warn = warn_customise
      )
    ),
    list(
      from = "SingleCellExperiment",
      to = "ReticulateAnnData",
      handler = .make_convert_handler(
        converter = function(from) {
          as_AnnData(from, output_class = "ReticulateAnnData")
        },
        warn = warn_customise
      )
    ),
    list(
      from = "SingleCellExperiment",
      to = "HDF5AnnData",
      handler = .make_abort_handler(
        from_class = "SingleCellExperiment",
        to_class = "HDF5AnnData",
        helper = "Use `as_AnnData(from, output_class = \"HDF5AnnData\", filename = <path>)` to provide the output file."
      )
    ),
    list(
      from = "AbstractAnnData",
      to = "SingleCellExperiment",
      handler = .make_convert_handler(
        converter = as_SingleCellExperiment,
        warn = warn_sce
      )
    )
  )

  .register_set_as_rules(single_cell_rules)
}

# Seurat coercion rules ------------------------------------------------------

if (rlang::is_installed("SeuratObject")) {
  warn_customise <- .format_control_recommendation("as_AnnData(...)")
  warn_seurat <- .format_control_recommendation("adata$as_Seurat(...)")

  seurat_rules <- list(
    list(
      from = "Seurat",
      to = "InMemoryAnnData",
      handler = .make_convert_handler(
        converter = function(from) {
          as_AnnData(from, output_class = "InMemoryAnnData")
        },
        warn = warn_customise
      )
    ),
    list(
      from = "Seurat",
      to = "ReticulateAnnData",
      handler = .make_convert_handler(
        converter = function(from) {
          as_AnnData(from, output_class = "ReticulateAnnData")
        },
        warn = warn_customise
      )
    ),
    list(
      from = "Seurat",
      to = "HDF5AnnData",
      handler = .make_abort_handler(
        from_class = "Seurat",
        to_class = "HDF5AnnData",
        helper = "Use `as_AnnData(from, output_class = \"HDF5AnnData\", filename = <path>)` to provide the output file."
      )
    ),
    list(
      from = "AbstractAnnData",
      to = "Seurat",
      handler = .make_convert_handler(
        converter = as_Seurat,
        warn = warn_seurat
      )
    )
  )

  .register_set_as_rules(seurat_rules)
}
