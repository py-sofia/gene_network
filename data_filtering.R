# Gene Network Mutual Information Threshold Analysis
# For filtering large gene interaction networks

library(dplyr)
library(ggplot2)
library(data.table)
library(readr)

# Function to analyze thresholds and suggest optimal cutoffs
analyze_network_thresholds <- function(file_path, sample_size = 100000) {
  
  cat("Loading CSV file...\n")
  
  # For very large files, use data.table for better memory efficiency
  tryCatch({
    # Try with data.table for speed and memory efficiency
    data <- fread(file_path, header = TRUE)
  }, error = function(e) {
    cat("Trying alternative CSV reading method...\n")
    # Fallback to readr
    data <- read_csv(file_path, show_col_types = FALSE)
  })
  
  # Rename columns for consistency
  colnames(data) <- c("gene_A", "gene_B", "mutual_info")
  
  cat(sprintf("Loaded %d gene pairs\n", nrow(data)))
  
  # Sample data if too large for analysis
  if (nrow(data) > sample_size) {
    cat(sprintf("Sampling %d rows for threshold analysis...\n", sample_size))
    data_sample <- data[sample(nrow(data), sample_size), ]
  } else {
    data_sample <- data
  }
  
  # Basic statistics
  cat("\n=== MUTUAL INFORMATION STATISTICS ===\n")
  mi_stats <- summary(data_sample$mutual_info)
  print(mi_stats)
  
  cat(sprintf("Standard deviation: %.6f\n", sd(data_sample$mutual_info, na.rm = TRUE)))
  cat(sprintf("Total unique genes: %d\n", 
              length(unique(c(data_sample$gene_A, data_sample$gene_B)))))
  
  # Calculate percentile-based thresholds
  percentiles <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999)
  thresholds <- quantile(data_sample$mutual_info, percentiles, na.rm = TRUE)
  
  cat("\n=== SUGGESTED THRESHOLDS ===\n")
  threshold_analysis <- data.frame(
    Percentile = paste0(percentiles * 100, "%"),
    Threshold = round(thresholds, 6),
    Edges_Remaining = sapply(thresholds, function(t) {
      round(sum(data_sample$mutual_info >= t, na.rm = TRUE) * nrow(data) / nrow(data_sample))
    }),
    Percent_Remaining = round(sapply(thresholds, function(t) {
      mean(data_sample$mutual_info >= t, na.rm = TRUE) * 100
    }), 2)
  )
  
  print(threshold_analysis)
  
  # Network size recommendations
  cat("\n=== NETWORK SIZE RECOMMENDATIONS ===\n")
  cat("For different network analysis purposes:\n")
  cat("- Small network (< 10K edges): Use 99.5% threshold or higher\n")
  cat("- Medium network (10K-100K edges): Use 95-99% threshold\n")
  cat("- Large network (100K-1M edges): Use 90-95% threshold\n")
  cat("- Very large network (> 1M edges): Use 80-90% threshold\n")
  
  # Create visualization
  create_threshold_plots(data_sample, thresholds)
  
  return(list(
    data = data,
    statistics = mi_stats,
    thresholds = threshold_analysis,
    sample_data = data_sample
  ))
}

# Function to create visualization plots
create_threshold_plots <- function(data_sample, thresholds) {
  
  # Distribution plot
  p1 <- ggplot(data_sample, aes(x = mutual_info)) +
    geom_histogram(bins = 100, fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = thresholds, color = "red", linetype = "dashed", alpha = 0.7) +
    scale_x_continuous(trans = "log10") +
    labs(title = "Distribution of Mutual Information Values",
         subtitle = "Red lines show potential thresholds",
         x = "Mutual Information (log scale)",
         y = "Frequency") +
    theme_minimal()
  
  print(p1)
  
  # Cumulative distribution
  p2 <- ggplot(data_sample, aes(x = mutual_info)) +
    stat_ecdf(color = "blue", size = 1) +
    geom_vline(xintercept = thresholds, color = "red", linetype = "dashed", alpha = 0.7) +
    scale_x_continuous(trans = "log10") +
    labs(title = "Cumulative Distribution of Mutual Information",
         subtitle = "Shows what fraction of edges remain at each threshold",
         x = "Mutual Information (log scale)",
         y = "Cumulative Probability") +
    theme_minimal()
  
  print(p2)
}

# Function to apply threshold and create filtered dataset
filter_network <- function(data, threshold, output_file = NULL) {
  
  cat(sprintf("Applying threshold: %.6f\n", threshold))
  
  # Filter data
  filtered_data <- data %>%
    filter(mutual_info >= threshold)
  
  cat(sprintf("Original edges: %d\n", nrow(data)))
  cat(sprintf("Filtered edges: %d\n", nrow(filtered_data)))
  cat(sprintf("Reduction: %.1f%%\n", 
              (1 - nrow(filtered_data)/nrow(data)) * 100))
  
  # Network statistics
  unique_genes <- length(unique(c(filtered_data$gene_A, filtered_data$gene_B)))
  avg_degree <- (2 * nrow(filtered_data)) / unique_genes
  
  cat(sprintf("Unique genes in network: %d\n", unique_genes))
  cat(sprintf("Average degree: %.2f\n", avg_degree))
  
  # Save filtered data if requested
  if (!is.null(output_file)) {
    write.csv(filtered_data, output_file, row.names = FALSE)
    cat(sprintf("Filtered data saved to: %s\n", output_file))
  }
  
  return(filtered_data)
}

# Function to find optimal threshold based on network size target
find_optimal_threshold <- function(data, target_edges = 50000, sample_size = 100000) {
  
  # Sample for faster computation
  if (nrow(data) > sample_size) {
    data_sample <- data[sample(nrow(data), sample_size), ]
    scale_factor <- nrow(data) / sample_size
  } else {
    data_sample <- data
    scale_factor <- 1
  }
  
  # Binary search for optimal threshold
  mi_values <- sort(data_sample$mutual_info, decreasing = TRUE)
  target_index <- round(target_edges / scale_factor)
  
  if (target_index > length(mi_values)) {
    cat("Target edge count is larger than available edges.\n")
    return(min(data_sample$mutual_info, na.rm = TRUE))
  }
  
  optimal_threshold <- mi_values[target_index]
  
  cat(sprintf("For approximately %d edges, use threshold: %.6f\n", 
              target_edges, optimal_threshold))
  
  return(optimal_threshold)
}

# Main execution example
# Replace 'your_file.xlsx' with your actual file path
main_analysis <- function(file_path) {
  
  # Step 1: Analyze thresholds
  results <- analyze_network_thresholds(file_path)
  
  # Step 2: Interactive threshold selection
  cat("\n=== THRESHOLD SELECTION ===\n")
  cat("Based on the analysis above, choose a threshold or use the suggestions below:\n")
  
  # Example: Filter for a medium-sized network (~50K edges)
  optimal_threshold <- find_optimal_threshold(results$data, target_edges = 50000)
  
  # Step 3: Apply filter
  filtered_network <- filter_network(results$data, 
                                     optimal_threshold, 
                                     "filtered_gene_network.csv")
  
  return(list(
    original_data = results$data,
    filtered_data = filtered_network,
    threshold_used = optimal_threshold
  ))
}

# Usage:
# Uncomment and modify the file path below
# results <- main_analysis("mi_results_ascp_clean.csv")

# Alternative: Quick threshold analysis only
# quick_analysis <- analyze_network_thresholds("your_file.csv")

cat("Script loaded. Use main_analysis('your_file_path.csv') to start analysis.\n")