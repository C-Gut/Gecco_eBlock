library(shiny)
library(DT)
library(tidyverse)

ui <- fluidPage(
  
  # App title
  titlePanel("eBlock Design"),
  
  # Horizontal line
  hr(),
  

  # Create a text input for the DNA sequence to be introduced by the user
  textInput('ex_len', 'example seq length', value = 41),
  textInput('frag_len', 'eblock fragment max. length', value = 870),
  textInput(inputId = "sequence", value = 'AAGGTCTCAAAAGGTCTCAAAAAGAGACCAAAAGAGACCAA', width = '100%', label = "Sequence"),
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
  DTOutput("frag_table")
)
