################################################################################
#                           libraries
################################################################################

library(readr)

library(gplots)



################################################################################
#                           functions
################################################################################



discretize_data <- function(data, num_bins = NULL) {
  
  # Input validation
  # add more
  if (is.null(num_bins)) {
    num_bins <- nclass.Sturges(data)
  }
  
  # so that all bins have the same width
  # Use consistent breaks across all variables
  range_data <- range(data, na.rm = TRUE)
  breaks <- seq(from=range_data[1], to=range_data[2], length.out = num_bins + 1)
  
  # Ensure rightmost bin includes maximum value
  breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps # to avoid issues with rightmost bin not including max value
  discretized <- cut(data, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  
  return(discretized)
}




calculate_probabilities <- function(data) {
  
  # Consider adding bias corrections for small sample sizes! 
  
  # Input validation
  if (length(data) == 0 || !is.integer(data)) {
    stop("Data must be a non-empty integer vector")
  }

  
  # frequency table (counts of each unique value)
  freq_table <- table(data)
  n <- sum(freq_table) # total number of observations
  
  probabilities <- freq_table / n

  return(probabilities)
  
}





calculate_entropy <- function(probabilities) {
  
  # calculate entropy
  # Use ifelse to set log2(p) = 0 when p = 0
  entropy <- -sum(probabilities * ifelse(probabilities > 0, log2(probabilities), 0)) 
  
  return(entropy)
  
}





calculate_joint_probabilities <- function(data_A, data_B) {
  
  # Input validation
  if (length(data_A) != length(data_B)) {
    stop("Data A and B must have the same length")
  }
  if (!is.integer(data_A) || !is.integer(data_B)) {
    stop("Data must be discrete integers. Use discretize_data() first.")
  }

  
  # Create joint frequency table
  joint_table <- table(data_A, data_B)
  n <- sum(joint_table)

  joint_probabilities <- joint_table/ n
  
  return(joint_probabilities)
  
}



calculate_joint_entropy <- function(joint_probabilities) {
  
  # calculate joint entropy
  # H(A,B) = -∑∑p(a,b)log2(p(a,b))
  # Use ifelse to set log2(p) = 0 when p = 0
  H_AB <- -sum(joint_probabilities * ifelse(joint_probabilities > 0, log2(joint_probabilities), 0))
  
  return(H_AB)
  
}




calculate_mutual_information <- function(probabilities_A, probabilities_B, joint_probabilities) {
  
  H_A <- calculate_entropy(probabilities_A)
  H_B <- calculate_entropy(probabilities_B)

  H_AB <- calculate_joint_entropy(joint_probabilities)
  
  # calculate mutual information
  MI <- H_A + H_B - H_AB

  return(MI)
  
}





################################################################################
#                           main
################################################################################

bins = 10


mean_value_data <- readr::read_tsv("GSE128816.top.table.tsv")
single_value_data <- read.csv("GSE128816_Processed_Data_File.csv")




# example usage:
gene_A_expression <- as.numeric(single_value_data[1, 3:13])
gene_B_expression <- as.numeric(single_value_data[2, 3:13])

probs_A <- calculate_probabilities(gene_A_expression, num_bins = bins)
probs_B <- calculate_probabilities(gene_B_expression, num_bins = bins)
joint_probs <- calculate_joint_probabilities(gene_A_expression, gene_B_expression, num_bins = bins)
MI <- calculate_mutual_information(probs_A, probs_B, joint_probs)



# verifying the result 
# ...
