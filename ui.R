library(shiny)
library(DT)
library(tidyverse)

ui <- fluidPage(
  
  # App title
  titlePanel("eBlock Design"),
  
  # Horizontal line
  hr(),
  
  # Create a text input for the DNA sequence to be introduced by the user
  textInput(inputId = "sequence", value = "AAGGTCTCAAAAGGTCTCAAAAAGAGACCAAAAGAGACCAA",
            label = "Sequence"
  ),
  
  # Horizontal line
  hr(),
  
  # Title for the outputs section
  HTML("<strong>Output</strong>"),
  
  # Line break
  br(),
  
  textOutput("bsaI_search"),
  
  DTOutput("table")
)

