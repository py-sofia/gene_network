
########################### number of valid genes ##############################
# Load your data
single_value_data <- read.csv("GSE128816_test.csv")

# Identify indices of genes with zero expression
zero_genes <- which(rowSums(single_value_data[, 3:12]) == 0)

# Get the total number of genes
n_total_genes <- nrow(single_value_data)

# Calculate the number of valid (non-zero) genes
n_valid_genes <- n_total_genes - length(zero_genes)

print(n_valid_genes)


########################## number of expected pairs ############################

# optimized
expected_pairs_opt = (n_valid_genes*(n_valid_genes-1))/2
print(expected_pairs_opt)

# original
expected_pairs_orig = (n_valid_genes-1)^2
print(expected_pairs_orig)



############################### manual calculation #############################

# Number of bins for discretization
bins <- 3

# Extract expression data for the first two genes
expr_A <- as.numeric(single_value_data[1, 3:12])
expr_B <- as.numeric(single_value_data[2, 3:12])

# Run the calculation steps
discrete_A <- discretize_data(expr_A, num_bins = bins)
discrete_B <- discretize_data(expr_B, num_bins = bins)

probs_A <- calculate_probabilities(discrete_A)
probs_B <- calculate_probabilities(discrete_B)

joint_probs <- calculate_joint_probabilities(discrete_A, discrete_B)

MI <- calculate_mutual_information(probs_A, probs_B, joint_probs)

print(paste("Manual MI for Gene 1 and 2:", MI))







#################### compare with bioconductor package #########################

# bioconductor calculation

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minet")

if (!require("infotheo", quietly = TRUE))
  install.packages("infotheo")

library(minet)
library(infotheo)

data(syn.data) # sample dataset
# syn.data has samples as rows and genes as columns

syn.data_transposed <- t(syn.data)
syn.data_discr <- discretize(syn.data_transposed, disc="equalwidth", nbins=3)

mim <- minet(syn.data_discr, method = "mrnet", estimator = "mi.empirical")



# implement_network_parallel calculation

source("implement_network_parallel_clean.R")

syn.data_transposed <- t(syn.data)

discretized <- as.data.frame(t(apply(syn.data_transposed, 1, discretize_data, num_bins = 3)))


run_mi_analysis(
  data = discretized,
  sample_cols = 2:ncol(discretized),
  bins = 3,
  chunk_size = 500000,
  out_file = "mi_results_ascp_clean.csv"
)