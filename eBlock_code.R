
########### VERSION 1: FINDS FIRST PATTERN IN SEQUENCE

# Function to find "GGTCTC" within a DNA sequence
# findPattern <- function(sequence) {
#   pattern <- "GGTCTC"  
#   # Pattern to search for
#   result <- regexpr(pattern, sequence)  
#   # Find the position of the pattern
#   
#   if (result == -1) {
#     cat("Pattern not found in the sequence.")
#   } else {
#     cat("Pattern found at position", result)
#   }
# }
# 
# # Accept DNA sequence from user
# cat("Enter the DNA sequence: ")
# sequence <- readline()
# 
# # Call the findPattern function
# findPattern(sequence)
# 
# 

########### VERSION 2: FINDS ALL OCCURRENCES OF GGTCTC PATTERN IN SEQUENCE

# # Function to find all occurrences of "GGTCTC" within a DNA sequence
# findAllPatterns <- function(sequence) {
#   pattern <- "GGTCTC"  
#   # Pattern to search for
#   #positions <- regexpr(pattern, sequence)  
#   # Find all positions of the pattern
#   matches <- gregexpr(pattern, sequence)  
#   # Alternative approach to find all positions
#   
#   # Extract the positions of all occurrences of the pattern
#   #all_positions <- unlist(positions[positions != -1])
#   all_matches <- unlist(matches)
#   
#   if (length(all_matches) == 0) {
#     cat("Pattern not found in the sequence.")
#   } else {
#     cat("Pattern found at positions:", paste(all_matches, collapse = ", "))
#     cat("\n")
#     
#     # Alternatively, you can use the following line to print the matches with start and end positions
#     # cat("Pattern found at matches:\n", as.data.frame(cbind(Start = all_matches, End = all_matches + nchar(pattern) - 1)))
#   }
# }
# 
# # Accept DNA sequence from user
# cat("Enter the DNA sequence: ")
# sequence <- readline()
# 
# # Call the findAllPatterns function
# findAllPatterns(sequence)



########### VERSION 3: FINDS ALL OCCURRENCES OF GGTCTC AND GAGACC PATTERNS IN SEQUENCE

# # Function to find all occurrences of multiple patterns within a DNA sequence
# 
# findAllPatterns <- function(sequence) {
#  
#   # Patterns to search for
#   patterns <- c("GGTCTC", "GAGACC") 
#   
#   all_matches <- vector("list", length(patterns))
#   
#   for (i in seq_along(patterns)) {
#     # Find all positions of the pattern
#     matches <- gregexpr(patterns[i], sequence)
#     print(matches)
#     all_matches[[i]] <- unlist(matches)
#   }
#   
#   any_pattern_found <- FALSE
#   
#   for (i in seq_along(all_matches)) {
#     if (length(all_matches[[i]]) > 0) {
#       any_pattern_found <- TRUE
#       cat("Pattern '", patterns[i], "' found at positions:", paste(all_matches[[i]], collapse = ", "), "\n")
#     }
#   }
#   
#   if (!any_pattern_found) {
#     cat("No patterns found in the sequence.")
#   }
# }
# 
# # Accept DNA sequence from user
# cat("Enter the DNA sequence: ")
# sequence <- readline()
# 
# # Call the findAllPatterns function
# findAllPatterns(sequence)

############## VERSION 4: MORE SIMPLE

#Create 2 data frames containing the positions where bsaI1 site "GGTCTC" and bsaI2 site "GAGACC" are found in the given input DNA sequence.
  #[[1]] is to make a data frame instead of a matrix

bsaI1.df <- str_locate_all(sequence, "GGTCTC")[[1]]
bsaI2.df <- str_locate_all(sequence, "GGTCTC")[[1]]


