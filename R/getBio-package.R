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
