#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @importFrom visNetwork renderVisNetwork visNetworkProxy visUpdateNodes 
#' @importFrom KEGGemUP kegg_to_graph map_results_to_graph
#' @importFrom golem get_golem_options
#' @noRd
app_server <- function(input, output, session) {
  data <- golem::get_golem_options("data")

  # Prepare DE and pathway lists
  de_list <- get_dea_names(data)
  pathway_list <- get_pathways_list()

  # Update selectors
  updateSelectInput(session, "de_selector", choices = de_list)
  updateSelectInput(session, "pathway_selector", choices = pathway_list)

  # Track selected nodes in the graph
  selected_ids <- reactiveVal(character())

  # Reactive DE based on selector
  selected_de <- reactive({
    req(input$de_selector)
    get_dea(data, input$de_selector)
  })

  # Reactive graph based on pathway and DE
  graph_data <- reactive({
    req(input$pathway_selector)
    selected_ids(character())
    KEGGemUP::kegg_to_graph(input$pathway_selector)
  })

  # Map DE to visualization graph
  visgraph <- reactive({
    req(graph_data(), input$de_selector)
    KEGGemUP::map_results_to_graph(graph_data(), data.frame(selected_de()), value_column = "log2FoldChange", feature_column = "ENTREZID")
  })

  # Render interactive graph
  output$graph_vis <- visNetwork::renderVisNetwork({
    req(visgraph())
    visgraph()
  })

  last_click_id <- reactiveVal(NULL)

  observeEvent(input$graph_click, {
    sel <- input$graph_click # or however you extract ID

    if (identical(sel, last_click_id())) {
      return()
    }

    last_click_id(sel)

    current <- selected_ids()

    if (is.null(current) || length(current) == 0) {
      updated <- sel
    } else if (sel %in% current) {
      updated <- setdiff(current, sel)
    } else {
      updated <- union(current, sel)
    }

    selected_ids(updated)
    
    g <- graph_data()
    node_ids <- V(g)$name
    
    if (length(node_ids) > 0) {
      nodes_to_update <- data.frame(
        id = node_ids,
        color.border = ifelse(node_ids %in% updated, "red", "black"),
        borderWidth = ifelse(node_ids %in% updated, 4, 1)
      )
      
      visNetwork::visNetworkProxy("graph_vis") %>%
        visNetwork::visUpdateNodes(nodes = nodes_to_update)
    }
    
    message("Selected IDs: ", paste(updated, collapse = ", "))
  })

  # Render selected features
  selected_features <- reactive({
    req(selected_ids())
    g <- graph_data()
    # get node IDs from igraph
    node_ids <- V(g)$name # assuming your node ids are stored in 'name'
    kegg_attr <- V(g)$KEGG # KEGG attributes
    # Filter KEGG attributes for selected nodes
    kegg_ids <- kegg_attr[node_ids %in% selected_ids(), drop = FALSE]
    get_features_from_selected_ids(data, kegg_ids)
  })

  # Render plot based on selected features
  output$feature_plot <- renderPlot({
    features <- selected_features()
    req(features)
    if (length(features) == 1) {
      plot_expression(unlist(features), data, row_column = "SYMBOL")
    } else if (length(features) > 1) {
      plot_heatmap(features, data, row_column = "SYMBOL")
    }
  })
}
