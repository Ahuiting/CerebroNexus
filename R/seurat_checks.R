## Argument checks shared by the functions that read a Seurat object. Their
## messages are user-facing; keep them unchanged.

## Seurat is a Suggests dependency and must be at least v5.0.
.checkSeuratAvailable <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(
      "The 'Seurat' package is needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  if (utils::packageVersion("Seurat") < "5.0.0") {
    stop(
      paste0(
        "The installed Seurat package is of version `",
        utils::packageVersion("Seurat"),
        "`, but at least v5.0 is required."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.checkSeuratObject <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop(
      paste0(
        "Provided object is of class `",
        paste(class(object), collapse = ", "),
        "` but must be of class 'Seurat'."
      ),
      call. = FALSE
    )
  }
  if (object@version < "3") {
    stop(
      paste0(
        "Provided Seurat object has version `",
        object@version,
        "` but must be at least 3.0."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.checkSeuratAssay <- function(object, assay) {
  if (assay %in% names(object@assays) == FALSE) {
    stop(
      paste0(
        'Specified assay slot `',
        assay,
        '` could not be found in provided Seurat object.'
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.checkSeuratGroups <- function(object, groups) {
  missing_groups <- groups[!groups %in% colnames(object@meta.data)]
  if (length(missing_groups)) {
    stop(
      paste0(
        "Group(s) `",
        paste0(missing_groups, collapse = '`, `'),
        "` were not ",
        "found in meta data of provided Seurat object. Only grouping variables ",
        "that are present in the meta data can be used."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## Per-cell UMI / gene count column for an assay, chosen by the same rule as
## the Builder's default (`qc_choices()` in inst/builder/inspect.R, which is
## sourced without the package and so keeps its own copy; a test compares
## the two): the first numeric column among the preferred names, otherwise the
## first column whose name looks like a count. NULL when there is none.
.seuratQcColumn <- function(meta, kind = c("umi", "gene"), assay) {
  kind <- match.arg(kind)
  umi_pattern <- paste0(
    "^(nCount|nUMI)|^total[._]?(count|UMI)s?$|",
    "^n[._]?counts?$|^(n[._])?UMIs?([._]counts?)?$"
  )
  gene_pattern <- paste0(
    "^(nFeature|nGene)|^total[._]?genes?$|",
    "^n[._]?genes?([._]by[._]counts)?$|^n[._]?features?$"
  )
  numeric_values <- function(values) {
    if (is.numeric(values)) {
      return(TRUE)
    }
    if (!(is.factor(values) || is.character(values))) {
      return(FALSE)
    }
    text <- trimws(as.character(values))
    missing <- is.na(values) | !nzchar(text)
    numeric <- suppressWarnings(as.numeric(text))
    any(!missing) && !any(!missing & !is.finite(numeric))
  }
  columns <- colnames(meta)
  qc_columns <- columns[grepl(
    paste0(umi_pattern, "|", gene_pattern),
    columns,
    ignore.case = TRUE
  )]
  numeric_cols <- qc_columns[vapply(
    qc_columns,
    function(column) numeric_values(meta[[column]]),
    logical(1)
  )]
  if (identical(kind, "umi")) {
    preferred <- c(
      paste0("nCount_", assay),
      "nUMI",
      "nCount_RNA",
      "total_counts",
      "n_counts",
      "total_UMIs",
      "UMI"
    )
    pattern <- umi_pattern
  } else {
    preferred <- c(
      paste0("nFeature_", assay),
      "nGene",
      "nFeature_RNA",
      "n_genes_by_counts",
      "n_genes",
      "total_genes"
    )
    pattern <- gene_pattern
  }
  relevant <- numeric_cols[grepl(pattern, numeric_cols, ignore.case = TRUE)]
  choices <- unique(c(intersect(preferred, numeric_cols), relevant))
  if (length(choices)) choices[[1L]] else NULL
}

## Levels of a grouping column: the factor levels, otherwise the unique values
## without NA, with a warning about the cells that have none. `na_text` is the
## warning text between the cell count and the group name.
.seuratGroupLevels <- function(
  values,
  group,
  na_text = ' cell(s) without group assignment (NA) for `'
) {
  if (is.factor(values)) {
    return(levels(values))
  }
  levels <- unique(values)
  if (anyNA(levels)) {
    levels <- stats::na.omit(levels)
    warning(
      'Found ',
      sum(is.na(values)),
      na_text,
      group,
      '`. These cells will be ignored during the analysis.',
      call. = FALSE
    )
  }
  levels
}
