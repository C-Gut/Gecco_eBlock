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
ATGCTAGCTAGCTAGGTCTCTAGCTAGCTAGCTAGCGGTCTCTAGCTAGCTAGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>Seq2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>Seq3
TACGTACGTACGTACGAGACCTTCGTACGTACGTACGTACGTACGTACGTACGTACGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
      label = "Sequence(s)", placeholder = "Paste FASTA sequences here...", width = '100%', height = '20vh'),
      hr(),
      HTML("<strong><u>Output</u></strong>"),
      br(),
      br(),
      textAreaInput("mod_seq", "Suggested sequence without BsaI sites:", width = "100%"),
      downloadButton("downloadCSV_wo_bsai", "Download Excel File"),
      hr(),
      HTML("<strong>Fragments</strong>"),
      br(), 
      # # Create a text input for the DNA sequence to be introduced by the user
      textInput('frag_len', 'eblock fragment max. length', value = 35),
      DTOutput("frag_table"),
      #downloadButton("downloadXLS_fragm", "Download Excel File"),
    ),
    
    # Right column (occupies the right half)
    column(width = 6,
      HTML("<strong><u>Extracted Sequences</u></strong>"),
      br(),
      DTOutput("fasta_table"),
    )
  ),


)