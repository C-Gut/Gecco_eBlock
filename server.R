

server <- function(input, output, session) {
  
  processed_input <- reactive({

    
    #Create 2 data frames containing the positions where bsaI1 site "GGTCTC" and bsaI2 site "GAGACC" are found in the given input DNA sequence.
      #[[1]] is to make a data frame instead of a matrix
        # to call the input use the id given in the UI after the $ sign. 
    
    bsaI1.df <- str_locate_all(input$sequence, "GGTCTC")[[1]]
    bsaI2.df <- str_locate_all(input$sequence, "GAGACC")[[1]]
    
     if (is.na(bsaI1.df[1]) == TRUE && is.na(bsaI2.df[1]) == TRUE ){
       bsaI_search <- "No BsaI sites were found"
       print(bsaI_search)
     } else {

     bsaI1.df <- as.data.frame(bsaI1.df)
     bsaI2.df <- as.data.frame(bsaI2.df)
    
    #Create extra columns to indicate whether fw or rv BsaI sites have been found
    bsaI1.df$GGTCTC <- TRUE
    bsaI1.df$GAGACC <- FALSE 
    bsaI2.df$GGTCTC <- FALSE
    bsaI2.df$GAGACC <- TRUE

    #Merge both data frames  
    bsaI.df <- rbind(bsaI1.df, bsaI2.df) %>% as.data.frame()
    
    #Create an extra column to indicate the position of the BsaI site in the codon.
      #%% 3 gives you the rest after dividing by 3. When the rest is 1 is in position 1, 
      #when it´s 2 in position 2, and when it´s 0 in position 3. 
        #Substitute 0s by 3s
    bsaI.df$position <- bsaI.df$start %% 3
    bsaI.df$position[bsaI.df$position == 0] <- 3
    
    #If there are not BsaI sites a message saying: 
    #"No BsaI sites found" appears and the rest of steps are done with the same input sequence
    
    #If there are BsaI sites, the table with the sites and positions found appear and a message appear asking to provide 
    #a sequence without the sites and an alternative sequence
    bsaI.df  
        }
## 1: GGT = GGC
## 2: GTC = GTG
## 3: TCT = TCC
    
    
     
  })
  
  output$table <- renderDT(processed_input())
}
