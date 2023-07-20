# Function to find "GGTCTC" within a DNA sequence
findPattern <- function(sequence) {
  pattern <- "GGTCTC"  # Pattern to search for
  result <- regexpr(pattern, sequence)  # Find the position of the pattern
  
  if (result == -1) {
    cat("Pattern not found in the sequence.")
  } else {
    cat("Pattern found at position", result)
  }
}

# Accept DNA sequence from user
cat("Enter the DNA sequence: ")
sequence <- readline()

# Call the findPattern function
findPattern(sequence)
