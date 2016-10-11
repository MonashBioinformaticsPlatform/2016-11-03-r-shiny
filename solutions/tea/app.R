
ui <- fluidPage(
    titlePanel("Ronald's exact test"),
    numericInput("tea", "Tea first", 3),
    numericInput("milk", "Milk first", 3),
    numericInput("tea_correct", "Tea first correctly called", 2),
    numericInput("milk_correct", "Milk first correctly called", 2),
    textOutput("p_text"))

server <- function(input, output, server) {
    x <- reactive( c(rep(0,input$tea), rep(1,input$milk)) )

    y <- reactive( c(rep(0,input$tea_correct), rep(1,input$tea-input$tea_correct), 
                     rep(0,input$milk-input$milk_correct), rep(1,input$milk_correct)) )

    statistic <- reactive( sum(x() == y()) )

    x_perms <- reactive( permutations(x()) ) # <- this is slow

    distribution <- reactive( colSums(x_perms() == y()) )

    p <- reactive( mean(distribution() >= statistic()) )

    output$p_text <- renderText( withProgress(message="Computing p", {
        paste0("p-value is ",p())
    }))
}

shinyApp(ui, server)
