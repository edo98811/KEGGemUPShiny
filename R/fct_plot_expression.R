#' Plot Expression for a Single feature in a SummarizedExperiment Object
#'
#' This function generates a plot of expression values for a specified feature (e.g., gene, protein, or metabolite)
#' across samples in a \code{SummarizedExperiment} object. The plot includes group information and uses different
#' background colors to indicate the data type (proteomics, transcriptomics, or metabolomics) based on the feature identifier.
#' Optionally, the function can export the plotting data instead of generating a plot.
#'
#' @param feature Character. The identifier of the feature to plot (e.g., gene symbol, UniProt ID, InChI key).
#' @param se A \code{SummarizedExperiment} object containing the count matrix and sample metadata.
#' @param export_data Logical. If \code{TRUE}, returns the plotting data as a data frame instead of a plot. Default is \code{FALSE}.
#' @param data_type Character. The type of data ("proteomics", "transcriptomics", "metabolomics", or "unknown"). Used for display purposes. Default is "unknown".
#' @param group_var Character. The column name in \code{colData(se)} to use for grouping samples. Default is "group".
#'
#' @return A \code{ggplot} object showing the expression values of the specified feature across groups, or a data frame if \code{export_data = TRUE}.
#'
#' @details
#' The function checks if the feature exists in the count matrix. If not, it returns a plot with an error message.
#' The plot includes a sina plot (via \code{ggforce::geom_sina}), boxplot, and sample labels (via \code{ggrepel::geom_text_repel}).
#' The background color of the plot is determined by the feature type, using helper functions \code{check_uniprot}, \code{check_ensembl}, and \code{check_inchi}.
#'
#' @import ggplot2
#' @import ggforce
#' @import ggrepel
#' @importFrom SummarizedExperiment assays colData
#' @noRd
plot_expression <- function(
    feature,
    se,
    row_column = NULL,
    export_data = FALSE,
    data_type = "unknown",
    group_var = "group",
    selected_assay = "counts") {
  # Use count matrix from SummarizedExperiment if not provided

  # Check if feature exists in count matrix
  if (is.null(se)) {
    return(ggplot() +
      annotate("text",
        x = 0.5, y = 0.5,
        label = paste0("Select a subset"),
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      theme_void())
  }

  count_matrix <- assays(se)$counts

  # GEt assay of feature to plot
  if (!is.null(row_column)) {
    se <- se[!is.na(rowData(se)[[row_column]])]
    heatmap_data <- assay(se[rowData(se)[[row_column]] == feature, ], selected_assay)
  } else {
    heatmap_data <- assay(se[rownames(rowData(se)) == feature, ], selected_assay)
  }

  # Check if feature exists in count matrix
  if (nrow(heatmap_data) == 0) {
    return(ggplot() +
      annotate("text",
        x = 0.5, y = 0.5,
        label = paste0("feature '", feature, "' not found in data"),
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      theme_void())
  }

  # GEt assay of feature to plot
  if (!is.null(row_column)) {
    plotting_data_raw <- assay(se)[rowData(se)[[row_column]] == feature, , drop = FALSE]
  } else {
    plotting_data_raw <- assay(se)[rowData(se)[[row_column]] == feature, , drop = FALSE]
  }

  # Create plotting data
  # plotting_data <- t(count_matrix[feature, , drop = FALSE])
  plotting_data <- t(plotting_data_raw)
  rownames(plotting_data) <- colnames(count_matrix)
  colnames(plotting_data) <- "Value"

  # Merge with group information
  plotting_data <- merge(plotting_data, colData(se)[c(group_var)], by = "row.names")
  plotting_data$group <- as.factor(plotting_data[[group_var]])

  # Set display name
  display_name <- feature

  # plotting_data <- plotting_data[order(plotting_data$group), ]

  # Export the plotting data if requested
  if (export_data) {
    return(plotting_data)
  }

  # Return plot
  return(plotting_data %>%
    ggplot(aes(x = group, y = Value, fill = group)) +
    ggforce::geom_sina(aes(color = group), size = 1.5) +
    ggrepel::geom_text_repel(aes(label = Row.names), size = 2) +
    geom_boxplot(alpha = 0.3) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 21),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    ggtitle(display_name) +
    xlab("") +
    ylab(display_name))
}
