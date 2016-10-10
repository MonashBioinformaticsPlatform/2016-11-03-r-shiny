# This file is generated from the corresponding .Rmd file



# =====
# Setup
# =====

install.packages("shiny")

source("https://bioconductor.org/biocLite.R")
biocLite(c(
    "BSgenome.Scerevisiae.UCSC.sacCer3",
    "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
    "GenomicRanges",
    "rtracklayer",
    "Gviz"))




# ============
# Hello, world
# ============

library(shiny)

ui1 <- fluidPage("Hello, world.")

server1 <- function(input, output, session) { }

app1 <- shinyApp(ui1, server1)
app1


class(ui1)
as.character(ui1)

class(app1)


# Various ways of running an app
runApp(app1)
runGadget(app1)
runGadget(app1, viewer=dialogViewer("App 1"))
print(app1)
app1


# ================
# input and output
# ================

# =========
# DataTable
# =========

## ----------------
## Upload challenge
## ----------------
# 
# 
# 
#
# =========================
# Reactive values save time
# =========================

# ==============================================
# What you can't see doesn't need to be computed
# ==============================================

## ------------------
## Reactive challenge
## ------------------
# 
# Use what you have just learned to make this app more responsive.
# 
# ...
# 
# 
#
## --------------------------------
## Genome browser challenge, part 1
## --------------------------------
# 
# The following code produces a diagram for a region of a genome. Your
# collaborator is asking you to make diagrams for a whole lot of
# different locations in the genome. Create a Shiny app to create these
# diagrams for them.
# 

library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

# We want to be able to interactively specify this location:
location <- as("chrI:140000-160000", "GRanges")

# You don't need to understand the rest of this code, just wrap it in Shiny.
genome <- BSgenome.Scerevisiae.UCSC.sacCer3
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

axis_track <- GenomeAxisTrack()
seq_track <- SequenceTrack(genome)
gene_track <- GeneRegionTrack(txdb, genome=genome, name="Genes", showId=TRUE, shape="arrow")

d1 <- rtracklayer::summary(BigWigFile("forward.bw"), location, 500, "max")[[1]]
d2 <- rtracklayer::summary(BigWigFile("reverse.bw"), location, 500, "max")[[1]]
d2$score <- -d2$score
data_track <- DataTrack(GRangesList(d1,d2), name="PAT-seq", type="l")

plotTracks(
    list(axis_track, seq_track, gene_track, data_track),
    chromosome=as.character(seqnames(location)),
    from=start(location),
    to=end(location))

# 
# 
#
# =================================
# Escaping from the reactive system
# =================================

# ===================
# Pushing back inputs
# ===================

## --------------------------------
## Genome browser challenge, part 2
## --------------------------------
# 
# Add buttons to your genome browser to navigate left and right, and/or
# zoom in and out.
# 
# 
# 
#
# =======
# Modules
# =======

## ----------------
## Module challenge
## ----------------
# 
# Adapt your genome browser to be a module.
# 
# There is a list of genes of particular interest. Make an app that
# shows a table of these genes. Use your genome browser module to
# display a selected gene.
# 
# 
# 
# ...solutions...
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
