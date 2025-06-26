# Auf der Console ausführen mit:  source("network_ascp.R", echo = FALSE)

################################################################################
#                           libraries
################################################################################

library(readr)          # schnelles CSV-/TSV-Einlesen
library(gplots)         # bleibt verfügbar
library(doSNOW)         # SOCK-Cluster + Fortschritts-Callback
library(foreach)        # Parallel-Schleifen
library(data.table)     # fwrite / rbindlist für schnellen I/O

################################################################################
#                           functions
################################################################################

# -------- Diskretisierung -----------------------------------------------------
discretize_data <- function(data, num_bins = NULL) {
  if (is.null(num_bins)) num_bins <- nclass.Sturges(data)
  rng    <- range(data, na.rm = TRUE)
  breaks <- seq(rng[1], rng[2], length.out = num_bins + 1)
  breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps
  cut(data, breaks = breaks, include.lowest = TRUE, labels = FALSE)
}

# -------- Wahrscheinlichkeiten / Entropie / MI -------------------------------
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

################################################################################
#                           main
################################################################################

# -------- Parameter -----------------------------------------------------------
bins          <- 3                   # Diskretisierungs-Bins
chunk_size    <- 500000              # 500'000 MI-Werte pro Schreibvorgang
out_file      <- "mi_results_ascp.csv"

# -------- Daten laden ---------------------------------------------------------
single_value_data <- read_csv("GSE128816.csv", show_col_types = FALSE)
sample_cols       <- 3:12
n_genes     <- nrow(single_value_data)
zero_genes  <- which(rowSums(single_value_data[, sample_cols]) == 0)
gene_idx    <- setdiff(seq_len(n_genes), zero_genes)
total_pairs <- choose(length(gene_idx), 2)

# -------- Parallel-Cluster ----------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 1)
cl      <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(cl)

clusterExport(
  cl,
  c("single_value_data", "sample_cols", "bins", "gene_idx",
    "discretize_data", "calculate_probabilities",
    "calculate_joint_probabilities", "calculate_entropy",
    "calculate_mutual_information"),
  envir = environment()
)

# -------- Anzeige der Parameter------------------------------------------------
cat("=====================\n")
cat("Genexpressionsanalyse\n")
cat("=====================\n")
cat("Geladene Gene:        ", format(n_genes, big.mark = "'"), "\n")
cat("Gene ohne Expression: ", format(length(zero_genes), big.mark = "'"), "\n")
cat("Gene mit Expression:  ", format(length(gene_idx), big.mark = "'"), "\n")
cat("Gen-Paare für MI:     ", format(total_pairs, big.mark = "'"), "\n")
cat("Buffer-Größe:         ", format(chunk_size, big.mark = "'"), " MI-Werte\n")
cat("---------------------\n")

# -------- Fortschrittsanzeige -------------------------------------------------
start_time <- Sys.time()
last_update_time <- start_time
last_pairs_done <- 0
update_interval <- 10000  # Häufigere Updates für kleine Tests

progress <- function(n_done) {
  pairs_done <- cumsum(length(gene_idx) - seq_along(gene_idx))[n_done]
  current_time <- Sys.time()
  
  if (pairs_done - last_pairs_done >= update_interval || pairs_done == total_pairs) {
    if (pairs_done > 0) {
      elapsed_time <- as.numeric(current_time - start_time, units = "secs")
      rate <- pairs_done / elapsed_time
      remaining_pairs <- total_pairs - pairs_done
      eta_seconds <- remaining_pairs / rate  # Zeit bis zum Ende
      
      # ETA formatieren (Zeit bis zum Ende, nicht Gesamtzeit)
      if (eta_seconds > 3600) {
        eta_formatted <- sprintf("%.1fh", eta_seconds / 3600)
      } else if (eta_seconds > 60) {
        eta_formatted <- sprintf("%.1fm", eta_seconds / 60)
      } else {
        eta_formatted <- sprintf("%.0fs", eta_seconds)
      }
      
      progress_percent <- 100 * pairs_done / total_pairs
      
      # Erweiterte Anzeige mit ETA und Rate
      cat(sprintf("\rFortschritt: %.1f%% (%s/%s) | ETA: %s | Rate: %.0f Paare/s", 
                  progress_percent,
                  format(pairs_done, big.mark = "'"),
                  format(total_pairs, big.mark = "'"),
                  eta_formatted,
                  rate))
      flush.console()
    }
    
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

# Buffer-Schreibfunktion
write_buffer <- function(new_data = NULL, force_write = FALSE) {
  # Neue Daten zum Buffer hinzufügen
  if (!is.null(new_data) && nrow(new_data) > 0) {
    mi_buffer <<- rbind(mi_buffer, new_data)
  }
  
  # Schreiben wenn Buffer voll oder force_write
  if (nrow(mi_buffer) >= chunk_size || (force_write && nrow(mi_buffer) > 0)) {
    buffer_count <<- buffer_count + 1
    
    cat(sprintf("\n[Buffer %d: Schreibe %s MI-Werte...]", 
                buffer_count, 
                format(nrow(mi_buffer), big.mark = "'")))
    
    fwrite(mi_buffer, 
           file = out_file, 
           append = header_written, 
           col.names = !header_written)
    
    total_written <<- total_written + nrow(mi_buffer)
    header_written <<- TRUE
    
    cat(sprintf(" ✓ (Gesamt: %s)\n", format(total_written, big.mark = "'")))
    
    # Buffer leeren
    mi_buffer <<- data.frame(gene1 = character(), 
                             gene2 = character(), 
                             MI = numeric(), 
                             stringsAsFactors = FALSE)
  }
}

# -------- Combine-Funktion für Buffer-Management -----------------------------
combine_and_buffer <- function(...) {
  results_list <- list(...)
  
  for (result in results_list) {
    if (is.data.frame(result) && nrow(result) > 0) {
      write_buffer(result)
    }
  }
  
  return(NULL)  # Nichts für die nächste Combine-Iteration
}

# -------- Paarweise MI-Berechnung (parallel) ----------------------------------
cat("Starte MI-Berechnung...\n")

# KORRIGIERTE foreach-Optionen
opts <- list(progress = progress)  # Entfernt chunkSize

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
    
    local_results  # Wird an combine_and_buffer weitergegeben
  }

# -------- Finaler Buffer-Schreibvorgang --------------------------------------
cat(sprintf("\n\nSchreibe finalen Buffer (%s MI-Werte)...\n", 
            format(nrow(mi_buffer), big.mark = "'")))
write_buffer(force_write = TRUE)

# -------- Abschluss -----------------------------------------------------------
total_elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
cat(sprintf("\n=== Abgeschlossen ===\n"))
cat(sprintf("Gesamtzeit: %.1f Minuten\n", total_elapsed / 60))
cat(sprintf("Buffer-Schreibvorgänge: %d\n", buffer_count))
cat(sprintf("Geschriebene MI-Werte: %s\n", format(total_written, big.mark = "'")))

stopCluster(cl)

# -------- Verifikation --------------------------------------------------------
if (file.exists(out_file)) {
  file_size_mb <- round(file.info(out_file)$size / 1024^2, 1)
  cat(sprintf("✓ Output-Datei: %s (%.1f MB)\n", out_file, file_size_mb))
  
  # Erste paar Zeilen zur Kontrolle
  sample_data <- fread(out_file, nrows = 3)
  cat("✓ Beispiel-Daten:\n")
  print(sample_data)
} else {
  cat("✗ FEHLER: Output-Datei wurde nicht erstellt!\n")
  
  # Debug-Information
  cat("Debug-Info:\n")
  cat("- Anzahl Gene mit Expression:", length(gene_idx), "\n")
  cat("- Erwartete Paare:", total_pairs, "\n")
  cat("- Buffer-Inhalt:", nrow(mi_buffer), "Zeilen\n")
}