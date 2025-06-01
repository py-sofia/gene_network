################################################################################
#                           libraries
################################################################################

library(readr)

library(gplots)



################################################################################
#                           functions
################################################################################



calculate_probabilities <- function(data, num_bins=NULL) {
  
  # Input validation
  if (length(data) == 0 || !is.numeric(data)) {
    stop("Data must be a non-empty numeric vector")
  }
  if (is.null(num_bins)) {
    # Sturges' formula
    num_bins <- nclass.Sturges(data)
  }
  if (num_bins <= 0) {
    stop("num_bins must be positive")
  }

  

  # create histogram
  range_data <- range(data)
  breaks_vec <- seq(range_data[1], range_data[2], length.out = num_bins + 1)
  hist_result <- hist(data, breaks = breaks_vec, plot = FALSE)
  
  # hist_result$counts = vector with number of values per bin
  probabilities <- hist_result$counts / sum(hist_result$counts)
  

  return(probabilities)
  
}





calculate_entropy <- function(probabilities) {
  
  # calculate entropy
  # Use ifelse to set log2(p) = 0 when p = 0
  entropy <- -sum(probabilities * ifelse(probabilities > 0, log2(probabilities), 0)) 
  
  return(entropy)
  
}





calculate_joint_probabilities <- function(data_A, data_B, num_bins=NULL) {
  
  # Input validation
  if (length(data_A) == 0 || !is.numeric(data_A)) {
    stop("Data A must be a non-empty numeric vector")
  }
  if (length(data_B) == 0 || !is.numeric(data_B)) {
    stop("Data B must be a non-empty numeric vector")
  }
  if (length(data_A) != length(data_B)) {
    stop("Data A and B must have the same length")
  }
  if (is.null(num_bins)) {
    # Sturges' formula
    num_bins <- nclass.Sturges(data_A) ## or data_B???
  }
  if (num_bins <= 0) {
    stop("num_bins must be positive")
  }

  
  
  # create joint histogram
  joint_hist <- hist2d(data_A, data_B, nbins = c(num_bins, num_bins), plot = TRUE)
  
  # joint probabilities
  joint_probabilities <- joint_hist$counts / sum(joint_hist$counts)

  
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
