
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

axis_track <- GenomeAxisTrack()
seq_track <- SequenceTrack(genome)
gene_track <- GeneRegionTrack(
    txdb, genome=genome, name="Genes", showId=TRUE, shape="arrow")

plot_genome <- function(location) {
    # Load data, at a reasonable level of detail for width(location)
    n <- min(width(location), 1000)
    d1 <- rtracklayer::summary(
        BigWigFile("forward.bw"), location, n, "max")[[1]]
    d2 <- rtracklayer::summary(
        BigWigFile("reverse.bw"), location, n, "max")[[1]]
    data_track <- DataTrack(
        d1, data=rbind(d1$score,-d2$score), groups=c(1,2), 
        name="PAT-seq", type=c("l"), col="#000000")
    
    plotTracks(
        list(axis_track, seq_track, gene_track, data_track),
        chromosome=as.character(seqnames(location)), 
        from=start(location), to=end(location))
}


# ==================================================================


ui <- fluidPage(
    titlePanel("Genome browser v2"),
    fluidRow(
        column(3, 
            textInput("location_str", "Genomic location", "chrI:140000-180000")),
        column(9,
            actionButton("go_left", "<<"),
            actionButton("go_right", ">>"),
            actionButton("zoom_in", "+"),
            actionButton("zoom_out", "-"))),
    plotOutput("genome_plot"))

server <- function(input, output, session) {
    #To observe input (eg exactly what brushing provides):
    #observe({ print(as.list(input)) })
    
    location <- reactive({ 
        GRanges(input$location_str, seqinfo=seqinfo(genome)) 
    })
    
    observeEvent(input$go_left, {
        amount <- max(1, width(location()) %/% 4) 
        new_location <- shift(location(), -amount)
        new_location <- trim(new_location)
        updateTextInput(session, "location_str", value=as.character(new_location))
    })
    
    observeEvent(input$go_right, {
        amount <- max(1, width(location()) %/% 4) 
        new_location <- shift(location(), amount)
        new_location <- trim(new_location)
        updateTextInput(session, "location_str", value=as.character(new_location))
    })
    
    observeEvent(input$zoom_in, {
        amount <- width(location()) %/% 4 
        new_location <- location()
        start(new_location) <- start(new_location) + amount
        end(new_location) <- end(new_location) - amount
        updateTextInput(session, "location_str", value=as.character(new_location))
    })

    observeEvent(input$zoom_out, {
        amount <- width(location()) %/% 2 
        new_location <- location()
        start(new_location) <- start(new_location) - amount
        end(new_location) <- end(new_location) + amount
        new_location <- trim(new_location)
        updateTextInput(session, "location_str", value=as.character(new_location))
    })

    output$genome_plot <- renderPlot({
        plot_genome( location() )
    })
}

shinyApp(ui, server)

