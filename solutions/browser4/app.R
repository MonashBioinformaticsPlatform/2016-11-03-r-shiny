
source("browser.R")

gene_df <- data.frame(
    gene=names(genes(txdb)),                 # txdb defined in browser.R
    location=as.character(genes(txdb)),
    row.names=NULL, stringsAsFactors=FALSE)

ui <- fluidPage(
    titlePanel("Gene table"),
    DT::dataTableOutput("genes"),
    h3("PAT-seq reads"),
    browser_ui("browser"))

server <- function(input,output,session) {
    output$genes <- DT::renderDataTable(
            server = TRUE,
            selection = "single",
            options = list(pageLength=10), {
        gene_df
    })

    observeEvent(input$genes_rows_selected, {
        if (length(input$genes_rows_selected) != 1)
            return()

        loc <- gene_df$location[input$genes_rows_selected]
        updateTextInput(session, "browser-location_str", value=loc)
    })

    callModule(browser_server, "browser")
}

shinyApp(ui, server)
