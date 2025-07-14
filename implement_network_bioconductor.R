
library("infotheo")
library("minet")
library("Rgraphviz")



data = read_csv("GSE128816_test.csv", show_col_types = FALSE)

gene_symbols <- data$GeneSymbol

# Remove the GeneSymbol and Description columns
counts_data <- data %>%
  select(starts_with("NormalizedCounts"))

# columns to rows, rows to columns
transposed_data <- t(counts_data)
colnames(transposed_data) <- gene_symbols

# convert to data frame
minet_ready_data <- as.data.frame(transposed_data)

result <- minet(minet_ready_data, method="mrnet", estimator="pearson",
                disc="equalwidth", nbins=3)

plot(as(result, "graphNEL"))
