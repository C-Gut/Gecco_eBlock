

server <- function(input, output, session) {
  
  processed_input <- reactive({
    #Create 2 data frames containing the positions where bsaI1 site "GGTCTC" and bsaI2 site "GAGACC" are found in the given input DNA sequence.
    #[[1]] is to make a data frame instead of a matrix
    # to call the input use the id given in the UI after the $ sign.
    
    # bsaI1.df <- str_locate_all(input$sequence, "GGTCTC")[[1]] %>% as.data.frame() # output of str_locate is matrix, need to convert to df
    # bsaI2.df <- str_locate_all(input$sequence, "GAGACC")[[1]] %>% as.data.frame()
    # 
    # if (nrow(bsaI1.df) == 0 && nrow(bsaI2.df) == 0) {
    #   output$bsaI_search <- renderText("No BsaI sites were found")
    #   return(NULL)  # Return NULL when no BsaI sites are found
    # } else if (nrow(bsaI1.df) > 0 && nrow(bsaI2.df) == 0) {
    #   output$bsaI_search <- renderText("BsaI1 sites were found")
    #   bsaI1.df$GGTCTC <- TRUE
    #   bsaI1.df$GAGACC <- FALSE
    #   bsaI1.df$suggested_sequence <-
    #     sub("GGT", "GGC", input$sequence)  # Create the suggested sequence for BsaI1
    #   return(bsaI1.df)  # Return the data frame for BsaI1 sites when only BsaI1 sites are found
    #   
    # } else if (nrow(bsaI1.df) == 0 && nrow(bsaI2.df) > 0) {
    #   output$bsaI_search <- renderText("BsaI2 sites were found")
    #   bsaI2.df$GGTCTC <- FALSE
    #   bsaI2.df$GAGACC <- TRUE
    #   bsaI2.df$suggested_sequence <-
    #     sub("GAGACC", "GAGACC", input$sequence)  # Create the suggested sequence for BsaI2
    #   return(bsaI2.df)  # Return the data frame for BsaI2 sites when only BsaI2 sites are found
    #   
    # } else {
    #   output$bsaI_search <-
    #     renderText("Both BsaI1 and BsaI2 sites were found")
    #   # Both data frames have BsaI sites, so merge them
    #   #Create extra columns to indicate whether fw or rv BsaI sites have been found
    #   bsaI1.df$GGTCTC <- TRUE
    #   bsaI1.df$GAGACC <- FALSE
    #   bsaI2.df$GGTCTC <- FALSE
    #   bsaI2.df$GAGACC <- TRUE
    #   
    #   #Merge both data frames
    #   bsaI.df <- rbind(bsaI1.df, bsaI2.df) %>% as.data.frame()
    #   
    #   #Create an extra column to indicate the position of the BsaI site in the codon.
    #   #%% 3 gives you the rest after dividing by 3. When the rest is 1 is in position 1,
    #   #when it´s 2 in position 2, and when it´s 0 in position 3.
    #   #Substitute 0s by 3s
    #   bsaI.df$position <- bsaI.df$start %% 3
    #   bsaI.df$position[bsaI.df$position == 0] <- 3
    #   
    #   # Create suggested sequences for both BsaI1 and BsaI2
    #   bsaI.df$suggested_sequence <- ifelse(
    #     bsaI.df$position == 1,
    #     sub("GGT", "GGC", input$sequence),
    #     sub("GAG", "GAA", input$sequence)
    #   )
    #   
    #   return(bsaI.df)  # Return the merged data frame with BsaI sites
    # }
    ## 1: GGT = GGC
    ## 2: GTC = GTG
    ## 3: TCT = TCC
    
    seq <- toupper(input$sequence)
    
    df <- rbind(
      seq %>% str_locate_all(., "GGTCTC") %>% as.data.frame() %>% mutate(dir = 'fw'),
      seq %>% str_locate_all(., "GAGACC") %>% as.data.frame() %>% mutate(dir = 'rv')
    ) %>% mutate(cod_pos = ((start + 2) %%3) + 1) 
    
    df$first_cod_to_change <- ifelse(df$dir == 'fw', 
                                    c('GTC','GGT','TCT')[(df$cod_pos %% 3) + 1], 
                                    c('AGA', 'GAG', 'GAC')[(df$cod_pos %% 3) + 1])
    df$change_cod_to <- ifelse(df$dir == 'fw', 
                               c('GTG','GGC','TCC')[(df$cod_pos %% 3) + 1], 
                               c('AGG', 'GAA', 'GAT')[(df$cod_pos %% 3) + 1])
    
    df[order(df$start),]
  
  })
  
  output$table <- renderDT({
    if (!is.null(processed_input())) {
      datatable(processed_input(), options = list(dom = 't'))  # Display the data table
    } else {
      data.frame()  # Return an empty data frame when no BsaI sites are found
    }
  })
}
