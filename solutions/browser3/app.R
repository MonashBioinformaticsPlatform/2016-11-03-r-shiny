
source("browser.R")

ui_browsermod <- fluidPage(
    titlePanel("Using a genome browser module"),
    browser_ui("browser"))

server_browsermod <- function(input,output,session) {
    callModule(browser_server, "browser")
}

shinyApp(ui_browsermod, server_browsermod)
