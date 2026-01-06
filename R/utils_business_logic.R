get_features_from_selected_ids <- function(dde, ids, col) {
  rd <- rowData(dde)

  all_ids <- unlist(lapply(ids, function(id) {
    strsplit(id, ";", fixed = TRUE)[[1]]
  }))
  
  if ("KEGG" %in% colnames(rd)) {
    out <- rd[[col]][rd$KEGG %in% all_ids]
  } else if ("ENTREZID" %in% colnames(rd)) {
    out <- rd[[col]][rd$ENTREZID %in% all_ids]
  } else {
    stop("Neither KEGG nor ENTREZID columns are present in rowData(dde).")
  }

  unique(out)
}

#' Helper functions to interact with KEGGemUP data
#' @importFrom KEGGemUP get_kegg_db
#' @noRd
get_pathways_list <- function() {
  KEGGemUP::get_kegg_db("pathway/mmu")[, 1]
}

#' @importFrom DeeDeeExperiment getDEANames
#' @noRd
get_dea_names <- function(dde) {
  tryCatch(
    {
      getDEANames(dde)
    },
    error = function(e) {
      message(paste0("Error retrieving DE analysis names: ", e$message))
      NULL
    }
  )
}

#' @importFrom DeeDeeExperiment getDEA
#' @noRd
get_dea <- function(dde, name) {
  tryCatch(
    {
      getDEA(dde, name, format = "original")
    },
    error = function(e) {
      message(paste0("Error retrieving DE analysis names: ", e$message))
      NULL
    }
  )
}
