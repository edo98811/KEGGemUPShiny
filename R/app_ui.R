#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bslib page_fillable page_sidebar sidebar layout_columns
#' @importFrom golem add_resource_path
#' @importFrom visNetwork visNetworkOutput
#'
#' @noRd
app_ui <- function(request) {
  page_fillable(
    padding = 0,
    page_sidebar(
      title = "KEGGemUPShiny",
      sidebar = sidebar(
        width = 400,
        title = "Feature Plot",
        selectInput("group_selector", "Select rowData col:", choices = NULL),
        selectInput("id_selector", "Select ID:", choices = NULL),
        selectInput("assay_selector", "Select Assay:", choices = NULL),
        plotOutput("feature_plot")
      ),
      layout_columns(
        selectInput("pathway_selector", "Select Pathway:", choices = NULL),
        selectInput("de_selector", "Select DE:", choices = NULL)
      ),
      visNetworkOutput("graph_vis", height = "600px")
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @importFrom KEGGemUP kegg_to_graph map_results_to_graph
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "KEGGemUPShiny"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}

# app_ui <- function(request) {

#   fluidPage(
#     sidebarLayout(
#       sidebarPanel(
#         h3("Selected Plot"),
#         plotOutput("feature_plot")
#       ),
#       mainPanel(
#         fluidRow(
#           column(
#             6,
#             selectInput("pathway_selector", "Select Pathway:", choices = NULL)
#           ),
#           column(
#             6,
#             selectInput("de_selector", "Select DE:", choices = NULL)
#           )
#         ),
#         visNetworkOutput("graph_vis", height = "600px")
#       )
#     )
#   )
# }
