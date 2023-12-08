################ 1. Set fixed variables ################

BsaSTART <- "CTGGTCTCGTGGT"
BsaSTOPCTTG <- "TAACTTGAGAGACCTG"
BsaMid2 <- "AGAGACCTG"
BsaMid1 <- "CTGGTCTCG"
BBseq1 <- "ATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGGGCAGCAGCCATCATCATCATCATCACGGCAGCGGCCTGGTGCCGCGCGGCAGCGCTGGT"
BBseq2 <- "CTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGACTCCAGCTTGGCTGTTTTGGCGGATGAGAGAAGATTTTCAGCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTGTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAACGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTTTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGAGTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATCGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGGTCATGGCTGCGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGGCAGCAGATCAATTCGCGCGCGAAGGCGAAGCGGCATGCATAATGTGCCTGTCAAATGGACGAAGCAGGGATTCTGCAAACCCTATGCTACTCCGTCAAGCCGTCAATTGTCTGATTCGTTACCAATTATGACAACTTGACGGCTACATCATTCACTTTTTCTTCACAACCGGCACGGAACTCGCTCGGGCTGGCCCCGGTGCATTTTTTAAATACCCGCGAGAAATAGAGTTGATCGTCAAAACCAACATTGCGACCGACGGTGGCGATAGGCATCCGGGTGGTGCTCAAAAGCAGCTTCGCCTGGCTGATACGTTGGTCCTCGCGCCAGCTTAAGACGCTAATCCCTAACTGCTGGCGGAAAAGATGTGACAGACGCGACGGCGACAAGCAAACATGCTGTGCGACGCTGGCGATATCAAAATTGCTGTCTGCCAGGTGATCGCTGATGTACTGACAAGCCTCGCGTACCCGATTATCCATCGGTGGATGGAGCGACTCGTTAATCGCTTCCATGCGCCGCAGTAACAATTGCTCAAGCAGATTTATCGCCAGCAGCTCCGAATAGCGCCCTTCCCCTTGCCCGGCGTTAATGATTTGCCCAAACAGGTCGCTGAAATGCGGCTGGTGCGCTTCATCCGGGCGAAAGAACCCCGTATTGGCAAATATTGACGGCCAGTTAAGCCATTCATGCCAGTAGGCGCGCGGACGAAAGTAAACCCACTGGTGATACCATTCGCGAGCCTCCGGATGACGACCGTAGTGATGAATCTCTCCTGGCGGGAACAGCAAAATATCACCCGGTCGGCAAACAAATTCTCGTCCCTGATTTTTCACCACCCCCTGACCGCGAATGGTGAGATTGAGAATATAACCTTTCATTCCCAGCGGTCGGTCGATAAAAAAATCGAGATAACCGTTGGCCTCAATCGGCGTTAAACCCGCCACCAGATGGGCATTAAACGAGTATCCCGGCAGCAGGGGATCATTTTGCGCTTCAGCCATACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCT"
sumoBB1 <- "ATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGGGCAGCAGCCATCATCATCATCATCACGGCAGCGGCCTGGTGCCGCGCGGCAGCGCTAGCATGTCGGACTCAGAAGTCAATCAAGAAGCTAAGCCAGAGGTCAAGCCAGAAGTCAAGCCTGAGACTCACATCAATTTAAAGGTGTCCGATGGATCTTCAGAGATCTTCTTCAAGATCAAAAAGACCACTCCTTTAAGAAGGCTGATGGAAGCGTTCGCTAAAAGACAGGGTAAGGAAATGGACTCCTTAAGATTCTTGTACGACGGTATTAGAATTCAAGCTGATCAGACCCCTGAAGATTTGGACATGGAGGATAACGATATTATTGAGGCTCACAGAGAACAGATTGGTGGT"
sumoBB2 <-"GGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGACTCCAGCTTGGCTGTTTTGGCGGATGAGAGAAGATTTTCAGCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTGTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAACGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTTTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGAGTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATCGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGGTCATGGCTGCGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGGCAGCAGATCAATTCGCGCGCGAAGGCGAAGCGGCATGCATAATGTGCCTGTCAAATGGACGAAGCAGGGATTCTGCAAACCCTATGCTACTCCGTCAAGCCGTCAATTGTCTGATTCGTTACCAATTATGACAACTTGACGGCTACATCATTCACTTTTTCTTCACAACCGGCACGGAACTCGCTCGGGCTGGCCCCGGTGCATTTTTTAAATACCCGCGAGAAATAGAGTTGATCGTCAAAACCAACATTGCGACCGACGGTGGCGATAGGCATCCGGGTGGTGCTCAAAAGCAGCTTCGCCTGGCTGATACGTTGGTCCTCGCGCCAGCTTAAGACGCTAATCCCTAACTGCTGGCGGAAAAGATGTGACAGACGCGACGGCGACAAGCAAACATGCTGTGCGACGCTGGCGATATCAAAATTGCTGTCTGCCAGGTGATCGCTGATGTACTGACAAGCCTCGCGTACCCGATTATCCATCGGTGGATGGAGCGACTCGTTAATCGCTTCCATGCGCCGCAGTAACAATTGCTCAAGCAGATTTATCGCCAGCAGCTCCGAATAGCGCCCTTCCCCTTGCCCGGCGTTAATGATTTGCCCAAACAGGTCGCTGAAATGCGGCTGGTGCGCTTCATCCGGGCGAAAGAACCCCGTATTGGCAAATATTGACGGCCAGTTAAGCCATTCATGCCAGTAGGCGCGCGGACGAAAGTAAACCCACTGGTGATACCATTCGCGAGCCTCCGGATGACGACCGTAGTGATGAATCTCTCCTGGCGGGAACAGCAAAATATCACCCGGTCGGCAAACAAATTCTCGTCCCTGATTTTTCACCACCCCCTGACCGCGAATGGTGAGATTGAGAATATAACCTTTCATTCCCAGCGGTCGGTCGATAAAAAAATCGAGATAACCGTTGGCCTCAATCGGCGTTAAACCCGCCACCAGATGGGCATTAAACGAGTATCCCGGCAGCAGGGGATCATTTTGCGCTTCAGCCATACTTTTCATACTCCCGCCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCT"

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
  seq_data <- lines[!str_detect(lines, "^>") & nchar(lines) > 0]
    
  
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

### func. to create n_chunks sequence fragments, returns vector of strings
# seq_chunks <- function(seq, n_chunks, max_len) {
#   # if only 1 chunk is requested, return input and stop
#   if (n_chunks == 1) {
#     return(seq)
#   } 
#   # create vector from string 
# 
#   seq <- str_split(seq, '')[[1]] 
#   print("seq:")
#   print(seq)
#   print("cut()")
#   subtr
#   print(cut(seq_along(seq), n_chunks, labels = FALSE))
#   print("split:")
#   print(split(seq, cut(seq_along(seq), n_chunks, labels = FALSE)))
#   # split(seq, cut(seq_along(seq), n_chunks, labels = FALSE)) # returns list where each item is vector chunk
#   print("145")
#   # Calculate the length of the first and last fragments
#   first_last_length <- max_len - 4
#   print("148")
#   # Calculate number of fragments except first and last
#   full_fragments <- function(seq, max_len) {
#        (calc_n_fragments(seq, max_len) - 2)
#   }
#   print("153")
#   print("full_fragments")
  
#########
split_vector = function(seq, n_chunks, x = 0) {
  # Ensure valid input
  if (n_chunks > seq || n_chunks < 1) {
    stop("Invalid number of parts")
  }
  
  # Basic split
  base_size = seq %/% n_chunks
  print("base size")
  print(base_size)
  remainder = seq %% n_chunks
  print(remainder)
  
  # Initial distribution
  parts = rep(base_size, n_chunks)
  parts[1:remainder] = parts[1:remainder] + 1
  
  # Adjust the first and last part
  if (n_chunks > 1) {
    parts[1] = max(1, parts[1] - 4)
    parts[n_chunks] = max(1, parts[n_chunks] - 3)  
  }
print(parts)
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
  print("212")
  print(result)
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
    # fasta_out_df <- as.data.frame(fasta_out, row.names = NULL)
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
    #fragm.df$p3_overhang <- substr(fragm.df$fragments, nchar(fragm.df$fragments) - 3, nchar(fragm.df$fragments))
    print(fragm.df[, "p5_overhang"])
    ### checks overhangs
    
    #1# are all overhangs unique?
    
    # Find duplicated values and report row numbers
    
    # Define the columns you want to compare for uniqueness
    column1_to_check <- fragm.df$p5_overhang
    column2_to_check <- fragm.df$p3_overhang
    
    # Combine both columns into a single vector for comparison
    combined_column <- c(column1_to_check, column2_to_check)
    
    # Find duplicated values and report row numbers
    duplicated_rows <-
      which(duplicated(column1_to_check) |
              duplicated(column1_to_check, fromLast = TRUE))
    
    # Check if there are any duplicated values
    fragm.df$p5_overhang_check_unique <- TRUE
    if (length(duplicated_rows) > 0) {
      cat("Duplicate values found in the following rows:\n")
      for (row_num in duplicated_rows) {
        cat("Row", row_num, ":", combined_column[row_num], "\n")
        fragm.df[row_num, "p5_overhang_check_unique"] <- FALSE
      }
    } else {
      cat("No duplicate values found between the two columns.\n")
      fragm.df$p5_overhang_check_unique <- TRUE
    }
    
    #2# is any overhang palindromic?
    
    # Function to check if a DNAString is palindromic
    is_palindromic <- function(sequence) {
      complement_sequence <- reverseComplement(DNAString(sequence))
      identical(sequence, as.character(complement_sequence))
    }
    
    # Check if each sequence in the 'Sequence' column is palindromic
    fragm.df$p5_overhang_check_palindrome <-
      sapply(fragm.df$p5_overhang, function(seq)
        ! is_palindromic(seq))
    
    
    #3# does any overhang have more than 2 repeats?
    
    # Function to check for repeated characters more than twice
    has_no_repeated_characters <- function(text) {
      !any(rle(strsplit(text, "")[[1]])$lengths > 2)
    }
    # Check if each value in the 'Text' column has repeated characters more than twice
    fragm.df$p5_overhang_check_repeats <-
      sapply(fragm.df$p5_overhang, has_no_repeated_characters)

    #4# pass all checks?

    # Create a new column 'test' based on the conditions
    fragm.df$test <- apply(fragm.df[,c("p5_overhang_check_unique", "p5_overhang_check_palindrome", "p5_overhang_check_repeats")], 
                           1, all)
    # If all the checks of the over hangs don´t pass
     cat("Do all fragments pass the OH checks?", all(fragm.df$test))
     
     if (!all(fragm.df$test)) {
       print("not all tests are passed")
     }
     
     # Find rows that don´t pass the checks and modify the fragment sequence by removing the las nucleotide
     # Show that in a different column, for now
     fragm.df$new_fragm <- fragm.df$fragments
     fragm.df$new_fragm[fragm.df$test == FALSE] <- substring(fragm.df$new_fragm[fragm.df$test == FALSE], 1, nchar(fragm.df$new_fragm[fragm.df$test == FALSE]) - 1)

    # Create a new column with the overhangs added to each fragment
    fragm.df$OH5prev <- c(fragm.df$p5_overhang[-1], NA)
    
    fragm.df$fragm_OH <- paste0(fragm.df$p5_Bsa, fragm.df$fragments, fragm.df$OH5prev, fragm.df$p3_Bsa)
    # New column for first fragmen with bsai 5' site
    fragm.df[1, "fragm_OH"] <- paste0(BsaSTART, fragm.df[1, "fragments"], fragm.df[1, "OH5prev"], BsaMid2)   
    fragm.df[nrow(fragm.df), "fragm_OH"] <- paste0(BsaMid1, fragm.df[nrow(fragm.df), "fragments"], BsaSTOPCTTG)
    print("fragm.df$fragm_OH")
    
    ## Create a new column with the length of the final fragments
    fragm.df$length_final_fragm <- nchar(fragm.df$fragm_OH)
    print(fragm.df)

    # fragm.df$test <-
    #   ifelse(
    #     fragm.df$p5_overhang_check_unique &
    #       fragm.df$p5_overhang_check_palindrome &
    #       fragm.df$p5_overhang_check_repeats,
    #     TRUE,
    #     FALSE
    #   )
    
    ### this only if all checks are true, otherwise change fragments


    fragm.df
  }
 ## output table for fragments
  output$frag_table <- renderDT({
    seqs_wo_bsa.l <- clean_multiple_mod_fasta()[[1]]
    print(input$frag_len)
    if (length(seqs_wo_bsa.l) > 0) {
      split_fragm.l <-
        lapply(seqs_wo_bsa.l, function(x)
        
          split_seq_in_chunks(x, as.numeric(input$frag_len)))
      processed_frag.l <-
        lapply(split_fragm.l, function(x)
          process_frags(x))
      fragments.df <- bind_rows(processed_frag.l, .id = "seq")
      as.data.frame(fragments.df)
      
      row.names(fragments.df) <- NULL
      colnames(fragments.df) <-
        c(
          "Name",
          "Length",
          "Fragments",
          "5'BasI",
          "3'BsaI",
          "5'OH",
          "5'OH unique",
          "5'OH palindrome",
          "5' repeats",
          "Pass all checks",
          "New fragm",
          "5´OH prev",
          "Fragm OH", 
          "Length final fragm"
        )
      
      fragments.df
    }
    
  })
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
}


# # return some info about fragments as text
# output$total_len <- renderText({
#   frag_input_seq <- input$mod_seq
#   mid_sites_to_add <- calc_break_points(seq = frag_input_seq, max_len = as.numeric(input$frag_len))
#   total_len <- nchar(frag_input_seq) + nchar(BsaSTART) + nchar(BsaSTOPCTTG) + (mid_sites_to_add * nchar(BsaMid1))
#   HTML(paste0("Input sequence length: &nbsp", as.character(nchar(frag_input_seq)), "<br>",
#               "Number of fragments &nbsp: &nbsp", as.character(ceiling(nchar(frag_input_seq)/(max_len = as.numeric(input$frag_len)))), "<br>",
#               "Final number of bases: &nbsp", as.character(total_len)))
