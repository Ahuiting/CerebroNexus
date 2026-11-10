## convertSeuratToCerebro() follows the exported assay for its per-cell totals
## when nUMI/nGene are left at their RNA defaults (Xenium, Atera, Visium, ...).

test_that(".seuratAssayCountColumn picks the assay's column only when needed", {
  meta <- data.frame(
    nCount_Atera = 1:3,
    nFeature_Atera = 1:3,
    nCount_RNA2 = 1:3
  )
  expect_identical(
    .seuratAssayCountColumn(meta, "nCount", "Atera", "nCount_RNA"),
    "nCount_Atera"
  )
  expect_identical(
    .seuratAssayCountColumn(meta, "nFeature", "Atera", "nFeature_RNA"),
    "nFeature_Atera"
  )
  ## RNA assay, absent column, or a present default: unchanged
  expect_identical(
    .seuratAssayCountColumn(meta, "nCount", "RNA", "nCount_RNA"),
    "nCount_RNA"
  )
  expect_identical(
    .seuratAssayCountColumn(meta, "nCount", "Xenium", "nCount_RNA"),
    "nCount_RNA"
  )
  meta$nCount_RNA <- 1:3
  expect_identical(
    .seuratAssayCountColumn(meta, "nCount", "Atera", "nCount_RNA"),
    "nCount_RNA"
  )
})

test_that("convertSeuratToCerebro exports a non-RNA assay without nUMI/nGene", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  set.seed(3)
  counts <- matrix(
    stats::rpois(60 * 40, 3),
    60,
    40,
    dimnames = list(paste0("G", 1:60), paste0("c", 1:40))
  )
  obj <- SeuratObject::CreateSeuratObject(
    methods::as(counts, "dgCMatrix"),
    assay = "Atera"
  )
  obj$grp <- rep(c("a", "b"), each = 20)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      stats::rnorm(80),
      40,
      2,
      dimnames = list(colnames(obj), c("umap_1", "umap_2"))
    ),
    key = "umap_",
    assay = "Atera"
  )
  out <- withr::local_tempdir()
  expect_no_error(convertSeuratToCerebro(
    obj,
    out,
    assay = "Atera",
    groups = "grp",
    add_most_expressed_genes = FALSE,
    verbose = FALSE
  ))
  crb <- list.files(out, "\\.crb$", full.names = TRUE)
  expect_length(crb, 1L)
  meta <- readCerebro(crb)$getMetaData()
  expect_true(all(c("nUMI", "nGene") %in% colnames(meta)))
  expect_equal(unname(meta$nUMI), unname(obj$nCount_Atera))

  ## An explicit nUMI is never replaced
  expect_error(
    convertSeuratToCerebro(
      obj,
      withr::local_tempdir(),
      assay = "Atera",
      groups = "grp",
      nUMI = "nCount_RNA",
      add_most_expressed_genes = FALSE,
      verbose = FALSE
    ),
    "nCount_RNA"
  )
})
