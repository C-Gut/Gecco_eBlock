library(shiny)
library(DT)
library(tidyverse)

ui <- fluidPage(
  
  # App title
  titlePanel("eBlock Design"),
  
  # Horizontal line
  hr(),
  
  # Create a text input for the DNA sequence to be introduced by the user
  textInput(inputId = "sequence", value = "AAGGTCTCAAAAGGTCTCAAAAAGAGACCAAAAAGAGACCAA",
            label = "Sequence"
  ),
  
  # Horizontal line
  hr(),
  
  # Title for the outputs section
  HTML("<strong>Output</strong>"),
  
  # Line break
  br(),
  DTOutput("table")
)

