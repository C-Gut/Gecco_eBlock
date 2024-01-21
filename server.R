################ 1. Set fixed variables ################

BsaSTART <- "CTGGTCTCGTGGT"
BsaSTOPCTTG <- "TAACTTGAGAGACCTG"
BsaMid2 <- "AGAGACCTG"
BsaMid1 <- "CTGGTCTCG"
BBseq1 <- "ATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGGGCAGCAGCCATCATCATCATCATCACGGCAGCGGCCTGGTGCCGCGCGGCAGCGCTGGT"
BBseq2 <- "CTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGACTCCAGCTTGGCTGTTTTGGCGGATGAGAGAAGATTTTCAGCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTGTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAACGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTTTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGAGTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATCGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGGTCATGGCTGCGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGGCAGCAGATCAATTCGCGCGCGAAGGCGAAGCGGCATGCATAATGTGCCTGTCAAATGGACGAAGCAGGGATTCTGCAAACCCTATGCTACTCCGTCAAGCCGTCAATTGTCTGATTCGTTACCAATTATGACAACTTGACGGCTACATCATTCACTTTTTCTTCACAACCGGCACGGAACTCGCTCGGGCTGGCCCCGGTGCATTTTTTAAATACCCGCGAGAAATAGAGTTGATCGTCAAAACCAACATTGCGACCGACGGTGGCGATAGGCATCCGGGTGGTGCTCAAAAGCAGCTTCGCCTGGCTGATACGTTGGTCCTCGCGCCAGCTTAAGACGCTAATCCCTAACTGCTGGCGGAAAAGATGTGACAGACGCGACGGCGACAAGCAAACATGCTGTGCGACGCTGGCGATATCAAAATTGCTGTCTGCCAGGTGATCGCTGATGTACTGACAAGCCTCGCGTACCCGATTATCCATCGGTGGATGGAGCGACTCGTTAATCGCTTCCATGCGCCGCAGTAACAATTGCTCAAGCAGATTTATCGCCAGCAGCTCCGAATAGCGCCCTTCCCCTTGCCCGGCGTTAATGATTTGCCCAAACAGGTCGCTGAAATGCGGCTGGTGCGCTTCATCCGGGCGAAAGAACCCCGTATTGGCAAATATTGACGGCCAGTTAAGCCATTCATGCCAGTAGGCGCGCGGACGAAAGTAAACCCACTGGTGATACCATTCGCGAGCCTCCGGATGACGACCGTAGTGATGAATCTCTCCTGGCGGGAACAGCAAAATATCACCCGGTCGGCAAACAAATTCTCGTCCCTGATTTTTCACCACCCCCTGACCGCGAATGGTGAGATTGAGAATATAACCTTTCATTCCCAGCGGTCGGTCGATAAAAAAATCGAGATAACCGTTGGCCTCAATCGGCGTTAAACCCGCCACCAGATGGGCATTAAACGAGTATCCCGGCAGCAGGGGATCATTTTGCGCTTCAGCCATACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCT"
sumoBB1 <- "ATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGGGCAGCAGCCATCATCATCATCATCACGGCAGCGGCCTGGTGCCGCGCGGCAGCGCTAGCATGTCGGACTCAGAAGTCAATCAAGAAGCTAAGCCAGAGGTCAAGCCAGAAGTCAAGCCTGAGACTCACATCAATTTAAAGGTGTCCGATGGATCTTCAGAGATCTTCTTCAAGATCAAAAAGACCACTCCTTTAAGAAGGCTGATGGAAGCGTTCGCTAAAAGACAGGGTAAGGAAATGGACTCCTTAAGATTCTTGTACGACGGTATTAGAATTCAAGCTGATCAGACCCCTGAAGATTTGGACATGGAGGATAACGATATTATTGAGGCTCACAGAGAACAGATTGGTGGT"
sumoBB2 <- "GGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGACTCCAGCTTGGCTGTTTTGGCGGATGAGAGAAGATTTTCAGCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTGTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAACGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTTTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGAGTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATCGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGGTCATGGCTGCGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGGCAGCAGATCAATTCGCGCGCGAAGGCGAAGCGGCATGCATAATGTGCCTGTCAAATGGACGAAGCAGGGATTCTGCAAACCCTATGCTACTCCGTCAAGCCGTCAATTGTCTGATTCGTTACCAATTATGACAACTTGACGGCTACATCATTCACTTTTTCTTCACAACCGGCACGGAACTCGCTCGGGCTGGCCCCGGTGCATTTTTTAAATACCCGCGAGAAATAGAGTTGATCGTCAAAACCAACATTGCGACCGACGGTGGCGATAGGCATCCGGGTGGTGCTCAAAAGCAGCTTCGCCTGGCTGATACGTTGGTCCTCGCGCCAGCTTAAGACGCTAATCCCTAACTGCTGGCGGAAAAGATGTGACAGACGCGACGGCGACAAGCAAACATGCTGTGCGACGCTGGCGATATCAAAATTGCTGTCTGCCAGGTGATCGCTGATGTACTGACAAGCCTCGCGTACCCGATTATCCATCGGTGGATGGAGCGACTCGTTAATCGCTTCCATGCGCCGCAGTAACAATTGCTCAAGCAGATTTATCGCCAGCAGCTCCGAATAGCGCCCTTCCCCTTGCCCGGCGTTAATGATTTGCCCAAACAGGTCGCTGAAATGCGGCTGGTGCGCTTCATCCGGGCGAAAGAACCCCGTATTGGCAAATATTGACGGCCAGTTAAGCCATTCATGCCAGTAGGCGCGCGGACGAAAGTAAACCCACTGGTGATACCATTCGCGAGCCTCCGGATGACGACCGTAGTGATGAATCTCTCCTGGCGGGAACAGCAAAATATCACCCGGTCGGCAAACAAATTCTCGTCCCTGATTTTTCACCACCCCCTGACCGCGAATGGTGAGATTGAGAATATAACCTTTCATTCCCAGCGGTCGGTCGATAAAAAAATCGAGATAACCGTTGGCCTCAATCGGCGTTAAACCCGCCACCAGATGGGCATTAAACGAGTATCCCGGCAGCAGGGGATCATTTTGCGCTTCAGCCATACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCT"

################ 2. Define costum functions that will be called later on inside the server ################

### function to make a random DNA sequence of length n, for testing
# example_seq <- function(n) {
#   sample(c('A', 'C', 'T', 'G'), n, replace = TRUE) %>% paste0(collapse = '')
# }


### Function to clean the FASTA input and create a data frame and a list. It accepts both fasta format and tab-separated field as input
clean_fasta <- function(input_text) {
  # Initialize empty vectors to store sequence names and sequences
  seq_names <- character(0)
  seq_data <- character(0)

  # Split the input text into lines, remove spaces
  lines <- unlist(strsplit(input_text, "\n"))
  lines <- gsub(" ", "", lines)

  # Define what will be the name and the sequence data
  # the nchar(lines) > 0 condition checks if the line is not empty (has a length greater than 0) after removing leading and trailing whitespace using trimws. Lines that are empty or contain only whitespace characters will be ignored when extracting seq_data.
  seq_names <- lines[str_detect(lines, "^>")]
  seq_data <- toupper(lines[!str_detect(lines, "^>") & nchar(lines) > 0])
    
  
  if (length(seq_names) != length(seq_data)) {
    print("problem with sequence")
    return()
   }
  
  # Create a data frame from the vectors
  fasta_df <- data.frame(Name = seq_names, Sequence = seq_data)
  # Create a list from the vectors
  fasta_list <- as.list(seq_data)
  names(fasta_list) <- seq_names
  
  return(list(fasta_list, fasta_df))
  }

### Bsai_locate function uses values of the fasta_list as input, finds BsaI sites, 
### and returns a data frame with info about when BsaI starts, ends, if fw or rv, and the position in the codon

bsai_locate <- function(seq) {
  fw_matches <- str_locate_all(seq, "GGTCTC")[[1]]
  rv_matches <- str_locate_all(seq, "GAGACC")[[1]]
  
  if (is.null(fw_matches) && is.null(rv_matches)) {
    # If no matches were found, create a row with NA values
    bsaI_locate_df <-
      data.frame(
        start = NA,
        end = NA,
        dir = NA,
        cod_pos = NA,
        seq_name = "Seq"
      )
  } else {
    if (length(fw_matches) == 0) {
      fw_df <- data.frame(start = NA,
                          end = NA,
                          dir = 'fw')
    } else {
      fw_df <- data.frame(start = fw_matches[, 1],
                          end = fw_matches[, 2],
                          dir = 'fw')
    }
    
    if (length(rv_matches) == 0) {
      rv_df <- data.frame(start = NA,
                          end = NA,
                          dir = 'rv')
    } else {
      rv_df <-
        data.frame(start = rv_matches[, 1],
                   end = rv_matches[, 2],
                   dir = 'rv')
    }
    
    bsaI_locate_df <- bind_rows(fw_df, rv_df) %>%
      mutate(cod_pos = ifelse(!is.na(start), ((start + 2) %% 3) + 1, NA),
             seq_name = "Seq")
  }
  
  bsaI_locate_df <- bsaI_locate_df[order(bsaI_locate_df$start), ]
  
  bsaI_locate_df <- bsaI_locate_df[!is.na(bsaI_locate_df$start), ]
  return(bsaI_locate_df)
  
}

### This function changes the necessary codons in the original sequence to remove BsaI sites found
### returns a text with the modified sequence
remove_bsai <- function(df, seq) {
  # Replace codons in the input sequence
  modified_sequence <- seq
  for (i in seq_len(nrow(df))) {
    start_pos <- df[i, "start"]
    end_pos <- df[i, "end"]
    replacement <- df[i, "change_cod_to"]
    pos_in_seq <- df[i, "calc_pos"]
    
    modified_sequence <-
      paste0(
        substring(modified_sequence, 1, pos_in_seq - 1),
        replacement,
        substring(modified_sequence, pos_in_seq + 3)
      )
  }
  modified_sequence
}

### func. to determine how many fragments are required given a sequence and a max. frag. length
# The first and last fragments need to have 4 and 3 nucleotides less than the rest, 
# Adding 8 nucleotides to the length of the sequences allows to calculate the amount of fragments by dividing by max_len

calc_n_fragments <- function(seq, max_len) {
  ceiling((nchar(seq) + 13) / max_len) # ceiling() rounds up to next integer
}

# ### function to calculate the break points in a sequence, given max. frag. length
# 
# calc_break_points <- function(seq, max_len) {
#   (calc_n_fragments(seq, max_len) - 1) * 2
# }


split_vector = function(seq, n_chunks, x = 0) {
  # Ensure valid input
  if (n_chunks > seq || n_chunks < 1) {
    stop("Invalid number of parts")
  }
  
  # Basic split
  base_size = seq %/% n_chunks
  remainder = seq %% n_chunks
  
  # Initial distribution
  parts = rep(base_size, n_chunks)
  parts[1:remainder] = parts[1:remainder] + 1
  
  # Adjust the first and last part
  if (n_chunks > 1) {
    parts[1] = max(1, parts[1] - 4)
    parts[n_chunks] = max(1, parts[n_chunks] - 3)  
  }

  # Ensure no part is empty
  parts = pmax(parts, 1)
  
 # # Recalculate middle parts if necessary
 #  if (sum(parts) != seq) {
 #    diff = seq - sum(parts)
 #    mid_indices = 2:(n_chunks-1)
 #    mid_parts = length(mid_indices)
 #    parts[mid_indices] = rep(base_size + diff %/% mid_parts, mid_parts)
 #    parts[mid_indices[1:(diff %% mid_parts)]] = parts[mid_indices[1:(diff %% mid_parts)]] + 1
 #  }  
  return(parts)
}
# Existing split_vector function here (as defined earlier)

seq_chunks = function(seq, n_chunks, x_first = 0, x_last = 0) {
  # Calculate the lengths for splitting
  split_lengths = split_vector(nchar(seq), n_chunks, x_first)
  print("201")
  
  # Split the string according to the lengths
  start = 1
  end = 0
  result = c()
  for (len in split_lengths) {
    end = start + len - 1
    result = c(result, substr(seq, start, end))
    start = end + 1
  }
  return(result)
}

### func. to create df with seq. fragments, given sequence and max. frag. length
split_seq_in_chunks <- function(seq, max_len) {
  if (seq == '' || !is.numeric(max_len) || is.na(max_len)) {
    return(data.frame())
  }
  chunks <- seq_chunks(seq = seq, n_chunks = calc_n_fragments(seq, max_len)) %>% lapply(., function(x)
    paste(x, collapse = ''))
  data.frame(
    length = lapply(chunks, function(x)
      str_length(x)) %>% unlist(),
    fragments = unlist(chunks)
  ) # convert list into data frame, can add more info (columns) about fragments
}

## Function to check the fragments that do not pass all the tests and fix them
# This function needs a data frame as an input
fix_checks <- function(fragm.df) {
  print('fixing')
  # repeat {
    # Check for false values
    false_values <- fragm.df$test == FALSE
    empty <- FALSE
    # If there are false values, manipulate the strings
    if (any(false_values) ||  empty) {
      # Loop through false values and update strings
      for (i in which(false_values)) {
        print(i)
        if (i > 1) {
          # Remove the first character from the current row
          current_string <- substr(fragm.df[i, 'fragments'], 2, nchar(fragm.df[i, 'fragments']))
          print(current_string)
          if (current_string == '') {
            print('empty')
            break()
          }
          # Add the removed character to the end of the previous row's string
          fragm.df$fragments[i - 1] <- paste0(fragm.df[i -1, 'fragments'], substr(fragm.df[i, 'fragments'], 1, 1))

          # Update the current row's string
          fragm.df$fragments[i] <- current_string
        }
      }
    } else {
      # If no false values, break the loop
      break
    }
  return(fragm.df)
}


################ 3. Server where the functions are called ################ 

server <- function(input, output, session) {

#######################################################################################
##################################### PARSE INPUT #####################################
#######################################################################################
  
 ### Use the clean_fasta function to remove spaces and empty lines from text input
 ### Accepts multiple fasta sequences as input and gives a table as output 
  
  clean_multiple_fasta <- reactive({
    
    input_text <- input$fasta_input
    
    # Here it checks if the input is a tab separated value (copied from excel),
      #if so, converts it into fasta format and then uses the clean_fasta function to process it

    if (str_detect(input_text, "\t")) {
      # Convert the input string to a tibble
      tsv_data <- read_tsv(input_text, col_names = FALSE)
      # Create FASTA format using paste
      input_text <- paste(">", tsv_data$X1, "\n", tsv_data$X2, sep = "")
      # Apply the clean_fasta function
      clean_fasta(input_text)

    } else {
      clean_fasta(input_text)
       }
  })
  
  output$fasta_table <- renderDT({
    clean_multiple_fasta()[[2]]
  })
  
  # for testing; create random sequences of length specified in input$ex_len
  # observeEvent(input$ex_len, {
  #   if (str_detect(input$ex_len, '^[0-9]+$')) {
  #     updateTextInput(session = getDefaultReactiveDomain(),
  #                     inputId = 'sequence',
  #                     value = example_seq(n = as.numeric(input$ex_len)))
  #   }
  # }, ignoreInit = TRUE)
  
#######################################################################################
################################## Remove BsaI sites ##################################
####################################################################################### 

  ### Use the bsai_locate function. Process the sequence: upper case, detect bsai sites, find position within codon, suggest swapped codon
  
  bsai_table <- function(seq) {
    df <- bsai_locate(seq = toupper(seq))
    df$first_cod_to_change <- ifelse(df$dir == 'fw',
                                     c('GTC', 'GGT', 'TCT')[(df$cod_pos %% 3) + 1],
                                     c('AGA', 'GAG', 'GAC')[(df$cod_pos %% 3) + 1])
    df$change_cod_to <- ifelse(df$dir == 'fw',
                               c('GTG', 'GGC', 'TCC')[(df$cod_pos %% 3) + 1],
                               c('AGG', 'GAA', 'GAT')[(df$cod_pos %% 3) + 1])
    
    # Create a new column with the position on the sequence where the codon that has to be changed starts.
    for (i in seq_len(nrow(df))) {
      if (df$cod_pos[i] == 3)
        df$calc_pos[i] <- df$start[i] + 1
      
      if (df$cod_pos[i] == 2)
        df$calc_pos[i] <- df$start[i] + 2
      
      if (df$cod_pos[i] == 1)
        df$calc_pos[i] <- df$start[i]
    }
    df
  }
  
  ### take the cleaned up input with multiple fasta sequences and create a df and a list of seq without bsai sites
  processed_input <- reactive({
    
    seq_list <- clean_multiple_fasta()[[1]]
    list_wo_bsai <- lapply(seq_list, function(x) {
      df <- bsai_table(x)
      modified_sequence <- remove_bsai(df, x)
      list(df, modified_sequence)
    })
    
    seqs <- lapply(list_wo_bsai, function(x)
      x[[2]])

    # create text output with input sequences without bsai sites (it is actually an input so it can be modified if the user wants)    
    fasta_out <-
      c(rbind(str_replace(names(seqs), '$', '_noBsai'), as.character(unlist(seqs)))) %>% paste(collapse = '\n')
    # c(rbind(str_replace(gsub(">", "", names(seqs)), '$', '_noBsai'), as.character(unlist(seqs)))) %>% paste(collapse = '\n')
    updateTextInput(session = getDefaultReactiveDomain(),
                    inputId = 'mod_seq',
                    value = fasta_out)
    fasta_out_df <- as.data.frame(fasta_out, row.names = NULL)
    list(list_wo_bsai, fasta_out)
  })
  
  # Download button with table without BsaI sites 
  output$downloadCSV_wo_bsai<- downloadHandler(
    filename = function() {
      paste("woBsa-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      #Use the clean_fasta function to create a data frame using fasta_out as input (which is the text output containing sequences without bsai sites)
        # Extract the names and sequences from processed_input()
        df_wo_bsai <- clean_fasta(processed_input()[[2]])[[2]]
        # Remove the > symbol from the data frame that can be downloaded as csv
        df_wo_bsai$Name <- sub(">", "", df_wo_bsai$Name)
        # Write the data frame to a CSV file
        write.xlsx(df_wo_bsai, file, sep = ";", rowNames = FALSE, quote = FALSE)
    }
  )
  observe(processed_input())
  
# create a table with the whole plasmid sequences  
 output$whole_seq.df <- renderDT({
    whole_seq <- clean_fasta(processed_input()[[2]])[[2]]
    if (input$plasmid == "pBAD"){
      whole_seq$Sequence <- paste0(BBseq1, whole_seq$Sequence, BBseq2)
      whole_seq <- whole_seq %>%
        mutate(Name = gsub("^>(.*)_noBsai$", "pBAD-\\1", Name))
      # Change column names
      colnames(whole_seq) <- c("Vector Name", "Vector Sequence")
    }
    if (input$plasmid == "pBAD SUMO"){
      whole_seq$Sequence <- paste0(sumoBB1, whole_seq$Sequence, sumoBB2)
      whole_seq <- whole_seq %>%
        mutate(Name = gsub("^>(.*)_noBsai$", "pBAD-SUMO-\\1", Name))
      # Change column names
      colnames(whole_seq) <- c("Vector Name", "Vector Sequence")
    }
    datatable(whole_seq, 
              options = list(
                scrollX = TRUE)
              )
    })
 
 # # Download button with vector 
 # output$downloadCSV_vector<- downloadHandler(
 #   filename = function() {
 #     paste("vector-", Sys.Date(), ".xlsx", sep="")
 #   }, 
 #   content = function(file) {
 #     whole_seq <- output$whole_seq.df
 #    # Write the data frame to a CSV file
 #    write.xlsx(whole_seq, file, sep = ";", rowNames = FALSE, quote = FALSE)
 #    }
  #)
 
  # # output table for bsai sites
  # # filling a new text input field with the new sequence
  # output$bsai_table <- renderDT({
  #   updateTextInput(session = getDefaultReactiveDomain(),
  #   inputId = 'mod_seq',
  #                   value = processed_input())
  #  datatable(processed_input()[[1]], options = list(dom = 't'))  # option removes (here) pointless search field
  # })
  
#######################################################################################
###################################### FRAGMENTS ######################################
####################################################################################### 

 ##Clean multiple fasta sequences coming from the modified sequence without bsaI sites
  clean_multiple_mod_fasta <- reactive({
      clean_fasta(input$mod_seq)
  })

  process_frags <- function(fragm.df) {
    fragm.df$p5_Bsa <- BsaMid1
    fragm.df$p3_Bsa <- BsaMid2
    fragm.df[1, "p5_Bsa"] <- BsaSTART
    fragm.df[nrow(fragm.df), "p3_Bsa"] <- BsaSTOPCTTG
    fragm.df$p5_overhang <- substr(fragm.df$fragments, 1, 4)

  ###### checks overhangs #######
    
    ##1## are all overhangs unique?
    
    # Find duplicated values and report row numbers
    
    # Define the columns you want to compare for uniqueness
    column1_to_check <- fragm.df$p5_overhang
    #column2_to_check <- fragm.df$p3_overhang
    
    # Combine both columns into a single vector for comparison
    #combined_column <- c(column1_to_check, column2_to_check)
    
    # Find duplicated values and report row numbers
    duplicated_rows <-
      which(duplicated(column1_to_check) |
              duplicated(column1_to_check, fromLast = TRUE))
    
    # Check if there are any duplicated values
    fragm.df$p5_overhang_check_unique <- TRUE
    if (length(duplicated_rows) > 0) {
      cat("Duplicate values found in the following rows:\n")
      for (row_num in duplicated_rows) {
        cat("Row", row_num, ":", column1_to_check[row_num], "\n")
        fragm.df[row_num, "p5_overhang_check_unique"] <- FALSE
      }
    } else {
      cat("No duplicate values found between the two columns.\n")
      fragm.df$p5_overhang_check_unique <- TRUE
    }
    
    ##2## is any overhang palindromic?
    
    # Function to check if a DNAString is palindromic
    is_palindromic <- function(sequence) {
      complement_sequence <- reverseComplement(DNAString(sequence))
      identical(sequence, as.character(complement_sequence))
    }
    
    # Check if each sequence in the 'Sequence' column is palindromic
    fragm.df$p5_overhang_check_palindrome <-
      sapply(fragm.df$p5_overhang, function(seq)
        ! is_palindromic(seq))
    
    ##3## does any overhang have more than 2 repeats?
    
    # Function to check for repeated characters more than twice
    has_no_repeated_characters <- function(text) {
      !any(rle(strsplit(text, "")[[1]])$lengths > 2)
    }
    # Check if each value in the 'Text' column has repeated characters more than twice
    fragm.df$p5_overhang_check_repeats <-
      sapply(fragm.df$p5_overhang, has_no_repeated_characters)
  
    ##4## pass all checks?

    # Create a new column 'test' based on the conditions
    fragm.df$test <- apply(fragm.df[,c("p5_overhang_check_unique", "p5_overhang_check_palindrome", "p5_overhang_check_repeats")], 
                           1, all)
    # If all the checks of the over hangs don´t pass
     cat("Do all fragments pass the OH checks?", all(fragm.df$test))

     # # Find rows that don´t pass the checks and modify the fragment sequence by removing the las nucleotide
     # # Show that in a different column, for now
     # fragm.df$new_fragm <- fragm.df$fragments
     # fragm.df$new_fragm[fragm.df$test == FALSE] <- substring(fragm.df$new_fragm[fragm.df$test == FALSE], 1, nchar(fragm.df$new_fragm[fragm.df$test == FALSE]) - 1)

    # Create a new column with the overhangs added to each fragment
    fragm.df$OH5prev <- c(fragm.df$p5_overhang[-1], NA)
    
    fragm.df$fragm_OH <- paste0(fragm.df$p5_Bsa, fragm.df$fragments, fragm.df$OH5prev, fragm.df$p3_Bsa)
    # New column for first fragmen with bsai 5' site
    fragm.df[1, "fragm_OH"] <- paste0(BsaSTART, fragm.df[1, "fragments"], fragm.df[1, "OH5prev"], BsaMid2)   
    fragm.df[nrow(fragm.df), "fragm_OH"] <- paste0(BsaMid1, fragm.df[nrow(fragm.df), "fragments"], BsaSTOPCTTG)

    # Create a new column with the length of the final fragments
    fragm.df$length_final_fragm <- nchar(fragm.df$fragm_OH)
    
    # # If there are false values, call the function to fix the fragments that do not pass all the tests
    # fixed_fragments.df <- data.frame()
    #   if (any(fragm.df$test) == FALSE){
    #     print("fixing false")
    #     fixed_fragments.df <- fix_checks(fragm.df)  
    #   }
  fragm.df
  }

 ## output table for fragments
  output$frag_table <- renderDT({
    seqs_wo_bsa.l <- clean_multiple_mod_fasta()[[1]]
    
      # The fragment length given by the user will refer to the final fragments with the addition of over hangs. 
      # Because the frag_len is used to split initial fragments without OH, we subtract 22, which is the length of OH added later
    
    if (as.numeric(input$frag_len) > 24){
      frag_len <- as.numeric(input$frag_len) -22      
    } else {
      break
    }

    
    if (length(seqs_wo_bsa.l) > 0) {
      split_fragm.l <-
        lapply(seqs_wo_bsa.l, function(x)
          split_seq_in_chunks(x, frag_len))
      processed_frag.l <-
        lapply(split_fragm.l, function(x)
          process_frags(x))
      fragments.df <- bind_rows(processed_frag.l, .id = "seq")
      as.data.frame(fragments.df)
      fragments.df$new_names <- 
      row.names(fragments.df) <- NULL
      colnames(fragments.df) <-
        c(
          "Name",
          "Length",
          "Fragments",
          "5'BsaI",
          "3'BsaI",
          "5'OH",
          "5'OH unique",
          "5'OH no palindrome",
          "5'OH no repeats",
          "Pass all checks",
          "5´OH prev",
          "Fragm OH", 
          "Length final fragm"
        )
      
      # Download button with table with fragments 
      output$downloadXLS_fragm<- downloadHandler(
        filename = function() {
          paste("fragm-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
          # Write the data frame to a XLSX file
          write.xlsx(fragments.df, file, sep = ";", rowNames = FALSE, quote = FALSE)
        }
      )
      
      # Change the order of the columns
      fragments.df <- fragments.df %>%
        select(1, "Fragm OH", "Length final fragm", "Pass all checks", "5'OH unique", "5'OH no palindrome", "5'OH no repeats", "Fragments", "Length", "5'BsaI", "3'BsaI", "5'OH", "5´OH prev")
      
      ##Change the fragment names
      # Split the "Name" column into parts
      name_parts <- strsplit(fragments.df$Name, "_")
      seq_identifiers <- sapply(name_parts, function(x) x[1])
      chunk_identifiers <- toupper(letters[sequence(table(seq_identifiers))])
      
      # Generate new row names
      new_row_names <- paste0(seq_identifiers, "_f", chunk_identifiers)
      
      # Update the values in the "Name" column
      fragments.df$Name <- new_row_names

      datatable(fragments.df, 
                options = list(
                  # autoWidth = TRUE,
                  # columnDefs = list(
                  #   list(width="100px", targets = c(2))
                  # )
                   scrollX = TRUE
                  # columnDefs = list(list(targets = c(2), className = "dt-left"))  
                )) %>% formatStyle(columns = c("Pass all checks"), target = "row", backgroundColor = styleEqual(c(TRUE, FALSE), c("white", "orange")))
    }
  })
}

# # return some info about fragments as text
# output$total_len <- renderText({
#   frag_input_seq <- input$mod_seq
#   mid_sites_to_add <- calc_break_points(seq = frag_input_seq, max_len = as.numeric(input$frag_len))
#   total_len <- nchar(frag_input_seq) + nchar(BsaSTART) + nchar(BsaSTOPCTTG) + (mid_sites_to_add * nchar(BsaMid1))
#   HTML(paste0("Input sequence length: &nbsp", as.character(nchar(frag_input_seq)), "<br>",
#               "Number of fragments &nbsp: &nbsp", as.character(ceiling(nchar(frag_input_seq)/(max_len = as.numeric(input$frag_len)))), "<br>",
#               "Final number of bases: &nbsp", as.character(total_len)))
