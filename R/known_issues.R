# This file contains the known issues that are currently present in the package.
# It can be used to generate documentation, but also throw warnings instead of errors
# in tests.
read_known_issues <- function() {
  check_requires("Reading known issues", "yaml")

  data <- yaml::read_yaml(system.file(
    "known_issues.yaml",
    package = "anndataR"
  ))

  map_dfr(
    data$known_issues,
    function(row) {
      expected_names <- c(
        "backend",
        "slot",
        "dtype",
        "process",
        "error_message",
        "description",
        "proposed_solution",
        "to_investigate",
        "to_fix"
      )
      if (!all(expected_names %in% names(row))) {
        cli_abort(c(
          "Unexpected columns in {.file known_issues.yaml}",
          "i" = "Expected columns: {.val {expected_names}}",
          "i" = "Actual columns: {.val {names(row)}}"
        ))
      }

      expand.grid(row)
    }
  )
}

is_known <- function(backend, slot, dtype, process, known_issues = NULL) {
  if (is.null(known_issues)) {
    known_issues <- read_known_issues()
  }

  filt <- rep(TRUE, nrow(known_issues))

  if (!is.null(backend)) {
    filt <- filt & known_issues$backend %in% backend
  }
  if (!is.null(slot)) {
    filt <- filt & known_issues$slot %in% slot
  }
  if (!is.null(dtype)) {
    filt <- filt & known_issues$dtype %in% dtype
  }
  if (!is.null(process)) {
    filt <- filt & known_issues$process %in% process
  }

  filt
}

message_if_known <- function(
  backend,
  slot,
  dtype,
  process,
  known_issues = NULL
) {
  if (is.null(known_issues)) {
    known_issues <- read_known_issues()
  }

  filt <- is_known(backend, slot, dtype, process, known_issues)

  if (any(filt)) {
    # take first
    row <- known_issues[which(filt)[[1]], ]

    paste0(
      "Known issue for backend '",
      row$backend,
      "', slot '",
      row$slot,
      "', dtype '",
      row$dtype,
      "', process '",
      row$process,
      "': ",
      row$description
    )
  } else {
    NULL
  }
}
