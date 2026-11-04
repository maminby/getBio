#' getBio: Helper functions for bulk and single-cell RNA-seq analysis
#'
#' A collection of convenience wrappers for expression-matrix preprocessing,
#' limma differential expression, GO/KEGG/GSEA/GSVA, common plots (heatmap,
#' PCA, volcano, violin, Venn, scatter) and Seurat single-cell workflows.
#'
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom methods new
"_PACKAGE"

# fallback for R < 4.4 (identical semantics to base::`%||%`)
`%||%` <- function(x, y) if (is.null(x)) y else x

# internal: read an assay matrix on both Seurat v4 (slot) and v5 (layer)
.get_assay_data <- function(object, what = "data", assay = "RNA") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::GetAssayData(object, assay = assay, layer = what)
  } else {
    SeuratObject::GetAssayData(object, assay = assay, slot = what)
  }
}
