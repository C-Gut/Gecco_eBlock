library(shiny)
library(DT)
library(tidyverse)
library(BiocManager)
library(Biostrings)
library(stringr)
library(dplyr)

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
      
      br(),
      hr(),
      HTML("<strong><u>Output</u></strong>"),    
      br(),
      br(),
      textAreaInput("mod_seq", "Suggested sequence without BsaI sites:", width = "100%"),
      hr(),
      DTOutput("frag_table"),     
    ),
    
    # Right column (occupies the right half)
    column(width = 6,
      HTML("<strong><u>Extracted Sequences</u></strong>"),
      br(),
      br(),
      DTOutput("fasta_table"),

    )
  )
)
