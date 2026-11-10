test_that("as_Seurat() fails gracefully", {
  expect_error(as_Seurat(), regexp = "adata.*is missing")
  expect_error(as_Seurat("foo"), regexp = "must be a <AbstractAnnData>")
})

skip_if_not_installed("Seurat")
library(Seurat)

known_issues <- read_known_issues()

ad <- generate_dataset(n_obs = 10L, n_vars = 20L, format = "AnnData")
ad$obsm[["X_pca"]] <- matrix(1:50, 10, 5)
ad$varm[["PCs"]] <- matrix(1:100, 20, 5)
ad$layers["counts"] <- ad$layers["numeric_csparse"]

# Remove column names to avoid warnings
for (name in c("integer", "numeric", "integer_with_nas", "numeric_with_nas")) {
  colnames(ad$obsm[[name]]) <- NULL
}

seu <- ad$as_Seurat()

test_that("as_Seurat retains number of observations and features", {
  expect_equal(nrow(seu), 20)
  expect_equal(ncol(seu), 10)

  # trackstatus: class=Seurat, feature=test_get_var_names, status=done
  expect_equal(rownames(seu), rownames(ad$var))
  # trackstatus: class=Seurat, feature=test_get_obs_names, status=done
  expect_equal(colnames(seu), rownames(ad$obs))
})

# trackstatus: class=Seurat, feature=test_get_obs, status=done
for (obs_key in colnames(ad$obs)) {
  test_that(paste0("as_Seurat retains obs key: ", obs_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("obs"),
      dtype = obs_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(obs_key %in% colnames(seu[[]]))
    expect_equal(
      seu[[]][[obs_key]],
      ad$obs[[obs_key]],
      info = paste0("obs_key: ", obs_key),
      ignore_attr = TRUE
    )
  })
}

# trackstatus: class=Seurat, feature=test_get_var, status=done
for (var_key in colnames(ad$var)) {
  test_that(paste0("as_Seurat retains var key: ", var_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("var"),
      dtype = var_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    active_assay <- seu[[DefaultAssay(seu)]]
    expect_true(var_key %in% colnames(active_assay[[]]))
    expect_equal(active_assay[[var_key]][[1]], ad$var[[var_key]])
  })
}


# trackstatus: class=Seurat, feature=test_get_layers, status=done
for (layer_key in names(ad$layers)) {
  test_that(paste0("as_Seurat retains layer: ", layer_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("layers"),
      dtype = layer_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    active_assay <- seu[[DefaultAssay(seu)]]
    expect_true(layer_key %in% active_assay[])
    expect_true(
      all.equal(
        as.matrix(t(active_assay[layer_key])),
        as.matrix(ad$layers[[layer_key]]),
        check.attributes = FALSE
      ),
      info = paste0("layer_key: ", layer_key)
    )
  })
}

# trackstatus: class=Seurat, feature=test_get_obsp, status=done
for (obsp_key in names(ad$obsp)) {
  test_that(paste0("as_Seurat retains obsp key: ", obsp_key), {
    graph_name <- paste0("RNA_", obsp_key)
    expect_true(graph_name %in% Graphs(seu))

    expect_equal(
      as.matrix(seu[[graph_name]]),
      as.matrix(ad$obsp[[obsp_key]]),
      ignore_attr = TRUE
    )
  })
}

# trackstatus: class=Seurat, feature=test_get_uns, status=done
for (uns_key in names(ad$uns)) {
  test_that(paste0("as_Seurat retains uns key: ", uns_key), {
    msg <- message_if_known(
      backend = "to_Seurat",
      slot = c("uns"),
      dtype = uns_key,
      process = "convert",
      known_issues = known_issues
    )
    skip_if(!is.null(msg), message = msg)

    expect_true(uns_key %in% names(Misc(seu)))
    expect_equal(Misc(seu)[[uns_key]], ad$uns[[uns_key]], ignore_attr = TRUE)
  })
}

test_that("as_Seurat retains pca dimred", {
  msg <- message_if_known(
    backend = "to_Seurat",
    slot = c("obsm"),
    dtype = "pca",
    process = "convert",
    known_issues = known_issues
  )

  skip_if(!is.null(msg), message = msg)

  # trackstatus: class=Seurat, feature=test_get_obsm, status=wip
  expect_true("X_pca" %in% Reductions(seu))
  expect_equal(
    Embeddings(seu, reduction = "X_pca"),
    ad$obsm[["X_pca"]],
    ignore_attr = TRUE
  )

  # trackstatus: class=Seurat, feature=test_get_varm, status=wip
  expect_equal(
    Loadings(seu, reduction = "X_pca"),
    ad$varm[["PCs"]],
    ignore_attr = TRUE
  )
})

test_that("as_Seurat works with list mappings", {
  expect_no_error(
    ad$as_Seurat(
      object_metadata_mapping = as.list(.as_Seurat_guess_object_metadata(ad)),
      layers_mapping = as.list(.as_Seurat_guess_layers(ad)),
      assay_metadata_mapping = as.list(.as_Seurat_guess_assay_metadata(ad)),
      reduction_mapping = as.list(.as_Seurat_guess_reductions(ad)),
      graph_mapping = as.list(.as_Seurat_guess_graphs(ad)),
      misc_mapping = as.list(.as_Seurat_guess_misc(ad))
    )
  )

  expect_error(
    ad$as_Seurat(
      reduction_mapping = list(numeric = "numeric_matrix")
    )
  )
})

test_that("as_Seurat works with a vector reduction_mapping", {
  expect_no_error(
    ad$as_Seurat(reduction_mapping = c(numeric = "numeric_matrix"))
  )
})

test_that("as_Seurat works with unnamed mappings", {
  expect_no_error(
    ad$as_Seurat(
      object_metadata_mapping = unname(.as_Seurat_guess_object_metadata(ad)),
      layers_mapping = c(
        na.omit(unname(.as_Seurat_guess_layers(ad))),
        data = NA
      ),
      assay_metadata_mapping = unname(.as_Seurat_guess_assay_metadata(ad)),
      graph_mapping = unname(.as_Seurat_guess_graphs(ad)),
      misc_mapping = unname(.as_Seurat_guess_misc(ad))
    )
  )
})

test_that("as_Seurat works with a 0 dim AnnData", {
  adata <- AnnData(X = as(matrix(nrow = 10, ncol = 0), "dgCMatrix"))

  expect_error(
    adata$as_Seurat(layers_mapping = c(counts = NA)),
    "must have non-zero dimensions"
  )
})

# Seurat replaces underscores in feature names with dashes when creating an
# object, see https://github.com/scverse/anndataR/issues/476
underscore_cases <- list(
  "a counts layer" = list(layers = "counts", x = FALSE),
  "a data layer" = list(layers = "data", x = FALSE),
  "counts and data layers" = list(layers = c("counts", "data"), x = FALSE),
  "X and a counts layer" = list(layers = "counts", x = TRUE)
)

for (case_name in names(underscore_cases)) {
  test_that(
    paste("as_Seurat keeps underscores in feature names with", case_name),
    {
      case <- underscore_cases[[case_name]]
      mat <- matrix(as.numeric(seq_len(200)), nrow = 10, ncol = 20)
      layers <- rep(list(mat), length(case$layers))
      names(layers) <- case$layers
      adata <- AnnData(
        X = if (case$x) mat else NULL,
        obs = data.frame(row.names = paste0("cell_", 1:10)),
        var = data.frame(row.names = paste0("gene_", 1:20), score = 1:20),
        layers = layers
      )

      warnings <- capture_warnings(seu <- adata$as_Seurat())

      # warn that the names are kept, instead of that they are replaced
      expect_match(warnings, "underscores", all = FALSE)
      expect_no_match(warnings, "replacing with dashes")

      expect_identical(rownames(seu), adata$var_names)
      expect_identical(colnames(seu), adata$obs_names)
      for (layer in c(case$layers, if (case$x) "X")) {
        expect_equal(
          as.matrix(LayerData(seu, layer = layer)),
          t(mat),
          ignore_attr = TRUE
        )
        expect_identical(
          rownames(LayerData(seu, layer = layer)),
          adata$var_names
        )
      }
      expect_identical(seu[["RNA"]][[]]$score, 1:20)

      # the names survive the round trip
      expect_identical(as_AnnData(seu)$var_names, adata$var_names)
    }
  )
}
