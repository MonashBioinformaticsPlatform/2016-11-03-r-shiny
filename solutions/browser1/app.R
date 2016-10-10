
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

#
# The actual plotting code can go in a function,
# so that it can be used from R as well.
#
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
    titlePanel("Genome browser v1"),
    textInput("location_str", "Genomic location", "chrI:140000-180000"),
    plotOutput("genome_plot"))

server <- function(input, output, session) {
    location <- reactive({ 
        GRanges(input$location_str, seqinfo=seqinfo(genome)) 
    })
    
    output$genome_plot <- renderPlot({
        plot_genome( location() )
    })
}

shinyApp(ui, server)

