library(shiny)
library(DT)
library(tidyverse)
library(BiocManager)
library(Biostrings)
library(stringr)
library(dplyr)
library(openxlsx)

ui <- fluidPage(
  
  # App title
  titlePanel("eBlock Design"),
  
  # Horizontal line
  hr(),
  
  fluidRow(
    # Left column (occupies the left half)
    column(width = 6,
      helpText("Paste your FASTA sequences below:"),
      textAreaInput(inputId = "fasta_input", value =
">Seq1
ATGCTAGCTAGCTAGGTCTCTAGCTAGCTAGCTAGCGGTCTCTAGCTAGCTAGC
>Seq2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
>Seq3
TACGTACGTACGTACGAGACCTTCGTACGTACGTACGTACGTACGTACGTACGTACGTACG", 
label = "Sequence(s)", placeholder = "Paste FASTA sequences here...", width = '100%', height = '20vh'),
      

      hr(),
      HTML("<strong><u>Output</u></strong>"),    
      br(),
      br(),
      textAreaInput("mod_seq", "Suggested sequence without BsaI sites:", width = "100%"),
downloadButton("downloadCSV_wo_bsai", "Download CSV"),
hr(),
HTML("<strong>Fragments</strong>"),
br(),
DTOutput("frag_table")
    ),
    
    # Right column (occupies the right half)
    column(width = 6,
      HTML("<strong><u>Extracted Sequences</u></strong>"),
      br(),
      DTOutput("fasta_table"),
    )
  ),


)