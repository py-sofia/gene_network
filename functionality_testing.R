
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

syn.data_discr <- discretize(syn.data, disc="equalfreq", nbins=3)

mim <- minet(syn.data_discr, method = "mrnet", estimator = "mi.empirical")



# implement_network_parallel calculation

library(readr) 
library(gplots)      
library(doSNOW)     
library(foreach)       
library(data.table)  
library(minet)

discretize_data <- function(data, num_bins = NULL) {
  if (is.null(num_bins)) num_bins <- nclass.Sturges(data)
  rng    <- range(data, na.rm = TRUE)
  breaks <- seq(rng[1], rng[2], length.out = num_bins + 1)
  breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps 
  cut(data, breaks = breaks, include.lowest = TRUE, labels = FALSE)
}

calculate_probabilities <- function(x) {
  if (length(x) == 0 || !is.integer(x))
    stop("x must be a non-empty integer vector")
  tbl <- table(x); tbl / sum(tbl)
}

calculate_joint_probabilities <- function(a, b) {
  if (length(a) != length(b)) stop("vectors must have same length")
  if (!is.integer(a) || !is.integer(b))
    stop("vectors must be integer bins (use discretize_data first)")
  jt <- table(a, b); jt / sum(jt)
}

calculate_entropy <- function(p) {
  v <- if (is.table(p) || is.array(p)) as.numeric(p) else p
  -sum(v * ifelse(v > 0, log2(v), 0))
}

calculate_mutual_information <- function(pA, pB, pAB) {
  calculate_entropy(pA) + calculate_entropy(pB) - calculate_entropy(pAB)
}

bins          <- 3                   
chunk_size    <- 500000              
out_file      <- "mi_results_ascp.csv"

single_value_data <- data(syn.data)
sample_cols <- 1:50
n_genes     <- nrow(single_value_data)
zero_genes  <- which(rowSums(single_value_data[, sample_cols]) == 0) 
gene_idx    <- setdiff(seq_len(n_genes), zero_genes)
total_pairs <- choose(length(gene_idx), 2)



# -------- Parallel-Cluster ----------------------------------------------------

# select all possible cores minus one for system stability
n_cores <- max(1, parallel::detectCores() - 1)

cl      <- makeCluster(n_cores, type = "SOCK")

# register the cluster for parallel processing
registerDoSNOW(cl)

# export variables and functions to the cluster
clusterExport(
  cl,
  c("single_value_data", "sample_cols", "bins", "gene_idx",
    "discretize_data", "calculate_probabilities",
    "calculate_joint_probabilities", "calculate_entropy",
    "calculate_mutual_information"),
  envir = environment()
)

# -------- Fortschrittsanzeige -------------------------------------------------
start_time <- Sys.time()
last_update_time <- start_time
last_pairs_done <- 0

# updates each 10'000 pairs
update_interval <- 10000  

progress <- function(n_done) {
  # cumulative number of pairs done
  pairs_done <- cumsum(length(gene_idx) - seq_along(gene_idx))[n_done]
  
  current_time <- Sys.time()
  
  if (pairs_done - last_pairs_done >= update_interval || pairs_done == total_pairs) {
    if (pairs_done > 0) { # avoid division by zero at the start
      elapsed_time <- as.numeric(current_time - start_time, units = "secs")
      rate <- pairs_done / elapsed_time
      remaining_pairs <- total_pairs - pairs_done
      # ETA: estimated time of arrival (task completion)
      eta_seconds <- remaining_pairs / rate  # Zeit bis zum Ende
      
      # format ETA
      if (eta_seconds > 3600) {
        eta_formatted <- sprintf("%.1fh", eta_seconds / 3600)
      } else if (eta_seconds > 60) {
        eta_formatted <- sprintf("%.1fm", eta_seconds / 60)
      } else {
        eta_formatted <- sprintf("%.0fs", eta_seconds)
      }
      
      # calculate percentage of completion
      progress_percent <- 100 * pairs_done / total_pairs
      
      # Erweiterte Anzeige mit ETA und Rate
      cat(sprintf("\rprogress: %.1f%% (%s/%s) | ETA: %s | rate: %.0f pairs/s", 
                  progress_percent,
                  format(pairs_done, big.mark = "'"),
                  format(total_pairs, big.mark = "'"),
                  eta_formatted,
                  rate))
      # avoid buffering
      flush.console()
    }
    
    # <<- to modify the global variable
    last_update_time <<- current_time
    last_pairs_done <<- pairs_done
  }
}

# -------- Buffer-System Initialisierung --------------------------------------
if (file.exists(out_file)) file.remove(out_file)

# Globale Buffer-Variablen - AUSSERHALB der foreach-Schleife!
mi_buffer <- data.frame(gene1 = character(), 
                        gene2 = character(), 
                        MI = numeric(), 
                        stringsAsFactors = FALSE)
header_written <- FALSE
total_written <- 0
buffer_count <- 0


write_buffer <- function(new_data = NULL, force_write = FALSE) {
  # add new data to buffer
  if (!is.null(new_data) && nrow(new_data) > 0) {
    mi_buffer <<- rbind(mi_buffer, new_data) # row-bound append
  }
  
  # write if buffer is full or if force_write = TRUE and there are new rows
  if (nrow(mi_buffer) >= chunk_size || (force_write && nrow(mi_buffer) > 0)) {
    # count the number of buffer writes
    buffer_count <<- buffer_count + 1
    
    cat(sprintf("\n[buffer %d: write %s MI-values...]", 
                buffer_count, 
                format(nrow(mi_buffer), big.mark = "'")))
    
    fwrite(mi_buffer, 
           file = out_file, 
           append = header_written, 
           col.names = !header_written)
    
    total_written <<- total_written + nrow(mi_buffer)
    header_written <<- TRUE
    
    cat(sprintf("(total: %s)\n", format(total_written, big.mark = "'")))
    
    # clear buffer
    mi_buffer <<- data.frame(gene1 = character(), 
                             gene2 = character(), 
                             MI = numeric(), 
                             stringsAsFactors = FALSE)
  }
}


# takes data from each parallel core and inputs them to the buffer
# crucial to combine the parallel workers!
combine_and_buffer <- function(...) {
  results_list <- list(...)
  
  for (result in results_list) {
    if (is.data.frame(result) && nrow(result) > 0) {
      write_buffer(result)
    }
  }
  
  return(NULL)  # very important, though why?
}

# -------- pairwise MI-calculation (parallel) ----------------------------------
cat("start MI calculation...\n")

opts <- list(progress = progress)

result <- foreach(i = gene_idx,
                  .combine = combine_and_buffer,
                  .options.snow = opts) %dopar% {
                    
                    local_results <- data.frame(gene1 = character(),
                                                gene2 = character(),
                                                MI = numeric(),
                                                stringsAsFactors = FALSE)
                    
                    for (j in gene_idx[gene_idx > i]) {
                      dA <- discretize_data(as.numeric(single_value_data[i, sample_cols]), bins)
                      dB <- discretize_data(as.numeric(single_value_data[j, sample_cols]), bins)
                      
                      MI <- calculate_mutual_information(
                        calculate_probabilities(dA),
                        calculate_probabilities(dB),
                        calculate_joint_probabilities(dA, dB))
                      
                      local_results[nrow(local_results) + 1, ] <- list(
                        single_value_data[i, 1, drop = TRUE],
                        single_value_data[j, 1, drop = TRUE],
                        MI
                      )
                    }
                    
                    local_results  # implicitly returned to .combine function
                  }

# write remaining data in buffer to file
cat(sprintf("\n\nwrite final buffer (%s MI-values)...\n", 
            format(nrow(mi_buffer), big.mark = "'")))
write_buffer(force_write = TRUE)

# -------- end -----------------------------------------------------------
total_elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
cat(sprintf("\n=== completed ===\n"))
cat(sprintf("total time: %.1f minutes\n", total_elapsed / 60))
cat(sprintf("buffer write count: %d\n", buffer_count))
cat(sprintf("total MI-values: %s\n", format(total_written, big.mark = "'")))

stopCluster(cl)

# -------- verification --------------------------------------------------------
if (file.exists(out_file)) {
  file_size_mb <- round(file.info(out_file)$size / 1024^2, 1)
  cat(sprintf("output-file: %s (%.1f MB)\n", out_file, file_size_mb))
  
  # first few rows
  sample_data <- fread(out_file, nrows = 3)
  cat("example data:\n")
  print(sample_data)
} else {
  cat("ERROR: Output file was not created!\n")
  
  # debug information
  cat("debug info:\n")
  cat("number of used genes:", length(gene_idx), "\n")
  cat("expected pairs", total_pairs, "\n")
  cat("buffer content", nrow(mi_buffer), "Zeilen\n")
}