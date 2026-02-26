#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[[1]] else "data/gasch2000.txt"
output_file <- if (length(args) >= 2) args[[2]] else "gasch_top10_heatmap.png"
transform_mode <- if (length(args) >= 3) args[[3]] else "auto" # auto|none|log2p1

if (!transform_mode %in% c("auto", "none", "log2p1")) {
  stop("transform_mode must be one of: auto, none, log2p1")
}

if (!file.exists(input_file)) {
  dir.create(dirname(input_file), recursive = TRUE, showWarnings = FALSE)
  message("Input file not found; downloading to: ", input_file)
  download.file(
    url = "https://www.shackett.org/files/gasch2000.txt",
    destfile = input_file,
    mode = "wb"
  )
}

raw_df <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "#",
  quote = "\"",
  na.strings = c("", "NA")
)

if (ncol(raw_df) < 4) {
  stop("Unexpected table format: expected at least 4 columns.")
}

gene_col <- names(raw_df)[1]
expr_col_start <- 4
max_expression_columns <- 30
last_expr_col <- min(ncol(raw_df), expr_col_start + max_expression_columns - 1)
expr_col_names <- names(raw_df)[expr_col_start:last_expr_col]

expr_df <- raw_df[, expr_col_start:last_expr_col, drop = FALSE]
expr_df <- as.data.frame(lapply(expr_df, function(x) suppressWarnings(as.numeric(x))))

numeric_ok <- vapply(expr_df, function(x) is.numeric(x), logical(1))
expr_df <- expr_df[, numeric_ok, drop = FALSE]

if (ncol(expr_df) == 0) {
  stop("No numeric expression columns found after conversion.")
}

expr_mat <- as.matrix(expr_df)
rownames(expr_mat) <- raw_df[[gene_col]]

na_count <- sum(is.na(expr_mat))
has_negative <- any(expr_mat < 0, na.rm = TRUE)

if (transform_mode == "auto") {
  if (has_negative) {
    message("Detected negative values: treating input as already log scale (no transform).")
  } else {
    message("No negative values detected: applying log2(x + 1).")
    expr_mat <- log2(expr_mat + 1)
  }
} else if (transform_mode == "log2p1") {
  message("Applying forced log2(x + 1) transform.")
  expr_mat <- log2(expr_mat + 1)
} else {
  message("Transform mode 'none': using values as-is.")
}

gene_mean <- rowMeans(expr_mat, na.rm = TRUE)
gene_sd <- apply(expr_mat, 1, sd, na.rm = TRUE)
cv <- ifelse(gene_mean == 0, NA_real_, gene_sd / gene_mean)

stats_df <- data.frame(
  gene = rownames(expr_mat),
  mean_expr = gene_mean,
  sd_expr = gene_sd,
  cv = cv,
  stringsAsFactors = FALSE
)

stats_df <- stats_df %>%
  filter(!is.na(cv), is.finite(cv)) %>%
  arrange(desc(cv))

if (nrow(stats_df) < 10) {
  stop("Fewer than 10 genes with finite CV were found.")
}

top10_genes <- stats_df$gene[1:10]
top_expr <- expr_mat[top10_genes, , drop = FALSE]

plot_df <- as.data.frame(top_expr)
plot_df$gene <- rownames(plot_df)
plot_long <- pivot_longer(plot_df, cols = -gene, names_to = "condition", values_to = "expression")
plot_long$gene <- factor(plot_long$gene, levels = top10_genes)

p <- ggplot(plot_long, aes(x = gene, y = condition, fill = expression)) +
  geom_tile(color = "black", linewidth = 0.15) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(256),
    na.value = "white"
  ) +
  labs(
    title = "gene top 10",
    x = "gene",
    y = "conditions",
    fill = "expression"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.2)
  )

ggsave(
  filename = output_file,
  plot = p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(stats_df[1:10, ], "top10_cv_genes.csv", row.names = FALSE)

message("Rows: ", nrow(raw_df))
message("Columns: ", ncol(raw_df))
message("Gene ID column: ", gene_col)
message("Expression columns used: ", ncol(expr_df), " (starting from original column 4)")
message("Missing values (NA) in expression matrix: ", na_count)
message("Output heatmap: ", output_file)
message("Top-10 CV table: top10_cv_genes.csv")
