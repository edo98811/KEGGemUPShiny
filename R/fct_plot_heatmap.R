#' Plot heatmap of expression values for selected features
#'
#' #' This function generates a heatmap of expression values for a set of selected features
#' across samples in a \code{SummarizedExperiment} object. The heatmap
#' can be customized with options for clustering, scaling, and winsorization.
#'
#' @param se A \code{SummarizedExperiment} object containing the expression data.
#' @param features A character vector of feature identifiers to include in the heatmap.
#' @param row_column An optional character string specifying the column in \code{rowData
#'
#' (se)} to match features against. If \code{NULL}, row names of \code{se} are used.
#' @param pathway An optional character string specifying the pathway name for title generation.
#' @param cluster_rows Logical. If \code{TRUE}, cluster the rows of the heatmap. Default is \code{TRUE}.
#' @param cluster_columns Logical. If \code{TRUE}, cluster the columns of the heatmap. Default is \code{FALSE}.
#' @param center_mean Logical. If \code{TRUE}, center the expression values by subtracting
#' the row means. Default is \code{TRUE}.
#' @param scale_row Logical. If \code{TRUE}, scale the rows to have unit variance. Default is \code{FALSE}.
#' @param winsorize_threshold An optional numeric value specifying the threshold for winsorization.
#' Values beyond this threshold will be capped. If \code{NULL}, no winsorization is applied. Default is \code{NULL}.
#' @param plot_title An optional character string specifying a custom title for the heatmap.
#' If \code{NULL}, a default title based on the pathway name is used.
#' @param export_data Logical. If \code{TRUE}, returns the processed data matrix used for the heatmap instead of generating the plot. Default is \code{FALSE}.
#' @param ... Additional arguments passed to the heatmap function.
#' @return A heatmap plot of expression values or a data matrix if \code{export_data = TRUE}.
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @noRd
plot_heatmap <- function(features, se,
                         row_column = NULL,
                         pathway = NULL,
                         cluster_rows = TRUE,
                         cluster_columns = FALSE,
                         center_mean = TRUE,
                         scale_row = FALSE,
                         winsorize_threshold = NULL,
                         plot_title = NULL,
                         export_data = FALSE,
                         ...) {
  # parameters check
  if (is.null(features) || length(features) == 0) {
    stop("features cannot be empty or NULL")
  }

  if (!is.null(winsorize_threshold)) {
    stopifnot(is.numeric(winsorize_threshold))
    stopifnot(winsorize_threshold >= 0)
  }

  # Select features to plot and check their presence
  if (!is.null(row_column)) {
    availablefeatures <- features[features %in% rowData(se)[[row_column]]]
    feature_indices <- match(availablefeatures, rowData(se)[[row_column]])
    heatmap_data <- assay(se)[feature_indices, , drop = FALSE]
    rownames(heatmap_data) <- features
  } else {
    availablefeatures <- features[features %in% rownames(se)]
    heatmap_data <- assay(se)[availablefeatures, , drop = FALSE]
  }

  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(heatmap_data, 1, var) == 0
  heatmap_data <- heatmap_data[!to_remove, , drop = FALSE]

  if (nrow(heatmap_data) < 2) {
    warning("Creating a heatmap with only one gene...")
  }

  hm_name <- "Expression \nvalues"

  # Handle pllotting of other data
  if (center_mean) {
    heatmap_data <- heatmap_data - rowMeans(heatmap_data)
    hm_name <- "Expression \nvalues"
  }

  if (scale_row) {
    heatmap_data <- t(scale(t(heatmap_data)))
    hm_name <- "Z-scores \nExpression \nvalues"
  }

  # If cutting extreme values
  if (!is.null(winsorize_threshold)) {
    # do the winsoring
    heatmap_data[heatmap_data < -winsorize_threshold] <- -winsorize_threshold
    heatmap_data[heatmap_data > winsorize_threshold] <- winsorize_threshold
  }

  # Generate title
  if (is.null(plot_title)) {
    title <- paste0("Signature heatmap ", pathway)
  } else {
    title <- plot_title
  }

  # Export data if requested
  if (export_data) {
    return(heatmap_data)
  }

  if (nrow(heatmap_data) == 0 || ncol(heatmap_data) == 0) {
    warning("No data available for the heatmap. Please check your input features.")
    return(NULL)
  }

  # Generate heatmap using base R graphics
  # heatmap(
  #   x = as.matrix(heatmap_data),
  #   main = title,
  #   col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  #   Rowv = NA,
  #   Colv = NA,
  #   dendrogram = "none",
  #   labRow = rownames(heatmap_data),
  #   margins = c(5, 10),
  #   cexRow = 0.8,
  #   cexCol = 0.8,
  #   ...
  # )
  # # Generate heatmap
  # ch <- ComplexHeatmap::Heatmap(
  #   matrix = heatmap_data,
  #   column_title = title,
  #   name = hm_name,
  #   col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  #   rect_gp = grid::gpar(col = "white", lwd = 0.5),
  #   cluster_rows = cluster_rows,
  #   cluster_columns = cluster_columns,
  #   row_labels = rownames(heatmap_data),
  #   show_heatmap_legend = FALSE,
  #   # row_labels = ifelse(use_symbol && !is.null(rowData(se)$SYMBOL),
  #   #   rowData(se)$SYMBOL[match(rownames(heatmap_data), rownames(rowData(se)))],
  #   #   rownames(heatmap_data)
  #   # ),
  #   ...
  # )

ch <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(heatmap_data),
  column_title = title,
  name = hm_name,

  # Color scale
  col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),

  # Cell borders
  rect_gp = grid::gpar(col = "white", lwd = 0.5),

  # Disable clustering and dendrograms
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,

  # Labels
  row_labels = rownames(heatmap_data),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),

  # Titles
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),

  # Legend
  show_heatmap_legend = FALSE,

  ...
)

  # return(ComplexHeatmap::draw(ch, merge_legend = TRUE))
  return(ch)
}

# iDEA -> mettere solo la funzione per creare la heatamp, le altre operazioni nel overview.
