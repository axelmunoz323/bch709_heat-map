library(data.table)
library(pheatmap)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1) Read top200 gene list from Example 2
top200 <- fread("results/yeast_stress_cv_top200.tsv")
gene_list <- top200$gene_id

# 2) Read Gasch 2000 full dataset
dt <- fread("data/gasch2000.txt", header = TRUE)
gene_ids <- dt[[1]]

meta_cols <- c("UID", "NAME", "GWEIGHT")
desc_col  <- names(dt)[3]
skip_cols <- c(meta_cols, desc_col)
cond_cols <- setdiff(names(dt), skip_cols)

mat <- as.matrix(dt[, ..cond_cols])
mode(mat) <- "numeric"

# 3) Subset to the top200 genes (keep order from gene_list)
idx <- match(gene_list, gene_ids)
if (any(is.na(idx))) {
  missing <- gene_list[is.na(idx)]
  stop(paste("Missing gene_ids in gasch2000.txt:", paste(missing, collapse = ", ")))
}
submat <- mat[idx, , drop = FALSE]
rownames(submat) <- gene_list

# 4) Keep only first 30 condition columns (original order)
if (ncol(submat) < 30) stop("Not enough condition columns to take first 30.")
submat <- submat[, 1:30, drop = FALSE]

# 5) Row-wise Z-score: (value - row_mean) / row_sd
row_means <- rowMeans(submat, na.rm = TRUE)
row_sds   <- apply(submat, 1, sd, na.rm = TRUE)
if (any(row_sds == 0 | is.na(row_sds))) stop("Some rows have sd=0 or NA; cannot z-score.")

zmat <- sweep(submat, 1, row_means, "-")
zmat <- sweep(zmat, 1, row_sds, "/")

# 6) Hierarchical clustering (euclidean distance, ward.D2)
d  <- dist(zmat, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

# 7) Cut tree k=4
clusters <- cutree(hc, k = 4)

# 8) Save cluster assignment TSV (sorted by cluster ascending)
cluster_table <- data.table(
  gene_id = names(clusters),
  cluster = as.integer(clusters)
)
setorder(cluster_table, cluster, gene_id)
fwrite(cluster_table, "results/cluster_assignment.tsv", sep = "\t")

# 9) Save clustered heatmap as PDF (8 x 12 inches)
annotation_row <- data.frame(cluster = factor(clusters))
rownames(annotation_row) <- names(clusters)

pdf("results/cv_top200_cluster_heatmap.pdf", width = 8, height = 12)
pheatmap(
  zmat,
  cluster_rows = hc,           # use our ward.D2 clustering
  cluster_cols = FALSE,        # keep original condition order
  annotation_row = annotation_row,
  show_rownames = FALSE,
  fontsize_col = 6,
  main = "CV Top200 Genes: Row-wise Z-score + Ward.D2 (k=4)"
)
dev.off()

cat("Saved: results/cv_top200_cluster_heatmap.pdf\n")
cat("Saved: results/cluster_assignment.tsv\n")
