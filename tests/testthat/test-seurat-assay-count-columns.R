## convertSeuratToCerebro() picks its per-cell UMI / gene count columns by the
## Builder's rule when nUMI/nGene are not given (Xenium, Atera, SCT, scanpy...).

qc_meta <- function(columns, n = 4L) {
  meta <- data.frame(row.names = paste0("c", seq_len(n)))
  for (column in columns) {
    meta[[column]] <- seq_len(n) + 1
  }
  meta
}

test_that(".seuratQcColumn follows the Builder's order of preference", {
  pick <- function(columns, kind, assay = "RNA") {
    .seuratQcColumn(qc_meta(columns), kind, assay)
  }
  ## RNA: nCount_RNA as before
  expect_identical(
    pick(c("nCount_RNA", "nFeature_RNA"), "umi"),
    "nCount_RNA"
  )
  expect_identical(
    pick(c("nCount_RNA", "nFeature_RNA"), "gene"),
    "nFeature_RNA"
  )
  ## The exported assay's own column comes first
  expect_identical(
    pick(c("nCount_RNA", "nCount_SCT"), "umi", "SCT"),
    "nCount_SCT"
  )
  expect_identical(pick("nCount_Atera", "umi", "Atera"), "nCount_Atera")
  expect_identical(pick("nFeature_Xenium", "gene", "Xenium"), "nFeature_Xenium")
  ## scanpy names
  expect_identical(
    pick(c("total_counts", "n_counts"), "umi"),
    "total_counts"
  )
  expect_identical(
    pick(c("n_genes", "n_genes_by_counts"), "gene"),
    "n_genes_by_counts"
  )
  ## Nothing that looks like a count
  expect_null(pick(c("percent_mt", "doublet_prob"), "umi"))
  ## A non-numeric column is never chosen
  meta <- qc_meta("nCount_RNA")
  meta$nCount_Atera <- letters[1:4]
  expect_identical(.seuratQcColumn(meta, "umi", "Atera"), "nCount_RNA")
})

test_that("the export picks the same nUMI/nGene as the Builder", {
  skip_if_not_installed("SeuratObject")
  builder_profile_source_runtime()
  builder_repo_source("inspect.R")
  counts <- methods::as(
    matrix(
      stats::rpois(30 * 6, 3),
      30,
      6,
      dimnames = list(paste0("G", 1:30), paste0("c", 1:6))
    ),
    "dgCMatrix"
  )
  cases <- list(
    list(assay = "RNA", columns = c("nCount_RNA", "nFeature_RNA", "nUMI")),
    list(
      assay = "SCT",
      columns = c("nCount_RNA", "nCount_SCT", "nFeature_SCT")
    ),
    list(assay = "Atera", columns = c("nCount_Atera", "nFeature_Atera")),
    list(assay = "RNA", columns = c("total_counts", "n_genes_by_counts")),
    list(assay = "RNA", columns = c("UMI", "nCount_ADT", "nFeature_ADT")),
    list(assay = "Xenium", columns = c("nUMI", "nGene"))
  )
  for (case in cases) {
    object <- suppressWarnings(SeuratObject::CreateSeuratObject(
      counts,
      assay = case$assay
    ))
    object@meta.data <- qc_meta(case$columns, n = 6L)
    rownames(object@meta.data) <- colnames(object)
    builder <- describe_seurat(object)$assay_profiles[[case$assay]]
    label <- paste(case$assay, paste(case$columns, collapse = ","))
    expect_identical(
      .seuratQcColumn(object@meta.data, "umi", case$assay) %||% NA_character_,
      builder$nUMI,
      label = label
    )
    expect_identical(
      .seuratQcColumn(object@meta.data, "gene", case$assay) %||% NA_character_,
      builder$nGene,
      label = label
    )
  }
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
  ## Stop here rather than report a follow-on error for a missing file.
  expect_length(crb, 1L)
  skip_if(length(crb) != 1L, "conversion failed; see the failure above")
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
