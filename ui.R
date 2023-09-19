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

  # Show text input box to paste fasta sequence and display the parsed sequence as a table
  sidebarLayout(
    sidebarPanel(
      helpText("Paste your FASTA sequences below:"),
      textAreaInput(inputId = "fasta_input", value =
#Seq1 has 2 GGTCTC in 15 and 37, Seq3 has 1 gagacc in 15
">Seq1
ATGCTAGCTAGCTAGGTCTCTAGCTAGCTAGCTAGCGGTCTCTAGCTAGCTAGC
>Seq2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
>Seq3
TACGTACGTACGTACGAGACCTTCGTACGTACGTACGTACGTACGTACGTACGTACGTACG", label = "Sequence(s)", placeholder = "Paste FASTA sequences here..."),
      actionButton("parse_button", "Parse")
    ),
    
    mainPanel(
      tableOutput("fasta_table")
    )
  ),
  
  # # Create a text input for the DNA sequence to be introduced by the user
  # textInput('ex_len', 'example seq length', value = 41),
  # textInput('frag_len', 'eblock fragment max. length', value = 870),
  # textInput(inputId = "sequence", value = 'AAGGTCTCAAAAGGTCTCAAAAAGAGACCAAAAGAGACCAA', width = '100%', label = "Sequence"),
  # Horizontal line
  hr(),
  
  # Title for the outputs section
  HTML("<strong>Output</strong>"),
  
  # Line break
  br(),hr(),
  
  htmlOutput("total_len"),
  # Line break
  hr(),
  
  DTOutput("bsai_table"),
  br(),
  hr(),

textInput("mod_seq", "Suggested sequence without BsaI sites:", width = "100%"),
  
  br(),
  hr(),
  DTOutput("frag_table")
)
