
################################################################################
#                           libraries
################################################################################

library(readr)          
library(gplots)         
library(doSNOW)         # parallel processing
library(foreach)        # parallel loops
library(data.table)     # fwrite 

################################################################################
#                           functions
################################################################################

# -------- discretization -----------------------------------------------------
discretize_data <- function(data, num_bins = NULL) {
  if (is.null(num_bins)) num_bins <- nclass.Sturges(data)
  rng    <- range(data, na.rm = TRUE)
  breaks <- seq(rng[1], rng[2], length.out = num_bins + 1)
  breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps
  cut(data, breaks = breaks, include.lowest = TRUE, labels = FALSE)
}

# -------- probabilities / entropy / MI -------------------------------
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


# -------- progress bar -------------------------------------------------

progress_bar_set_up <- function() {
  
  start_time <- Sys.time()
  assign("start_time", start_time, envir = .GlobalEnv)
  
  last_update_time <- start_time
  assign("last_update_time", last_update_time, envir = .GlobalEnv)
  
  last_pairs_done <- 0
  assign("last_pairs_done", last_pairs_done, envir = .GlobalEnv)
  
  update_interval <- 10000  
  assign("update_interval", update_interval, envir = .GlobalEnv)
}


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


# -------- buffering --------------------------------------

buffer_set_up <- function() {
  # set up the output file
  if (file.exists(out_file)) file.remove(out_file)
  
  mi_buffer <- data.frame(gene1 = character(), 
                          gene2 = character(), 
                          MI = numeric(), 
                          stringsAsFactors = FALSE)
  assign("mi_buffer", mi_buffer, envir = .GlobalEnv)
  
  header_written <- FALSE
  assign("header_written", header_written, envir = .GlobalEnv)
  
  total_written <- 0
  assign("total_written", total_written, envir = .GlobalEnv)
  
  buffer_count <- 0
  assign("buffer_count", buffer_count, envir = .GlobalEnv)
  
}


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


set_up <- function(bins, chunk_size=500000, out_file="mi_results_ascp.csv",
                   data, sample_cols) {
  
  assign("bins", bins, envir = .GlobalEnv)
  assign("chunk_size", chunk_size, envir = .GlobalEnv)
  assign("out_file", out_file, envir = .GlobalEnv)
  assign("single_value_data", data, envir = .GlobalEnv)
  assign("sample_cols", sample_cols, envir = .GlobalEnv)
  assign("n_genes", nrow(data), envir = .GlobalEnv)

  
  # returns vector with indices of genes with zero expression across all samples
  zero_genes  <- which(rowSums(single_value_data[, sample_cols]) == 0) 
  assign("zero_genes", zero_genes, envir = .GlobalEnv)
  
  # indices of genes with at least some expression
  gene_idx    <- setdiff(seq_len(n_genes), zero_genes)
  assign("gene_idx", gene_idx, envir = .GlobalEnv)
  
  # total number of gene pairs for MI calculation
  total_pairs <- choose(length(gene_idx), 2)
  assign("total_pairs", total_pairs, envir = .GlobalEnv)
  
}


parallel_cluster <- function() {
  # select all possible cores minus one for system stability
  n_cores <- max(1, parallel::detectCores() - 1)
  
  cl      <- makeCluster(n_cores, type = "SOCK")
  assign("cl", cl, envir = .GlobalEnv)
  
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
}



pairwise_MI_calculation <- function() {
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
  
  
}


end <- function() {
  
  total_elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  stopCluster(cl)
  write_buffer(force_write = TRUE)
  
}


# wrapper function
run_mi_analysis <- function(data, sample_cols, bins = 3, chunk_size = 500000, 
                            out_file = "mi_results_ascp.csv") {
  set_up(bins, chunk_size, out_file, data, sample_cols)
  progress_bar_set_up()
  buffer_set_up()
  parallel_cluster()
  pairwise_MI_calculation()
  end()
}



################################################################################
#                           main
################################################################################

run_mi_analysis(
  data = read_csv("GSE128816_test.csv", show_col_types = FALSE),
  sample_cols = 3:12,
  bins = 3,
  chunk_size = 500000,
  out_file = "mi_results_ascp_clean.csv"
)

