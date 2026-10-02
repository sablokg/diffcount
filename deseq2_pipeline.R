#!/usr/bin/env Rscript
# ==============================================================================
# DESeq2 Differential Expression Pipeline
# ------------------------------------------------------------------------------
# Input  : raw counts matrix (genes x samples) + sample metadata table
# Output : normalized counts, DE results table, PCA plot, sample-distance
#          heatmap, MA plot, volcano plot, top-gene heatmap
#
# Usage:
#   Rscript deseq2_pipeline.R
#   (edit the CONFIG section below, or source this script and adjust objects)
# ==============================================================================

# ---- 0. CONFIG ---------------------------------------------------------------

counts_file    <- "counts.csv"      # rows = genes, cols = samples (raw integer counts)
metadata_file  <- "metadata.csv"    # rows = samples, must include a "condition" column
sample_id_col  <- "sample_id"       # column in metadata that matches counts colnames
condition_col  <- "condition"       # the grouping variable to test (e.g. treatment vs control)
ref_level      <- "control"         # reference/baseline level of condition_col
contrast_group <- "treatment"       # level to compare against the reference
alpha_cutoff   <- 0.05              # FDR (padj) significance threshold
lfc_cutoff     <- 1                 # |log2FoldChange| threshold used for volcano/heatmap highlighting
top_n_heatmap  <- 30                # number of top DE genes to show in the heatmap
out_dir        <- "deseq2_results"  # output directory for all files

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. PACKAGES ---------------------------------------------------------------

required_pkgs <- c("DESeq2", "pheatmap", "RColorBrewer", "ggplot2",
                    "ggrepel", "apeglm", "EnhancedVolcano")

install_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    BiocManager::install(missing, update = FALSE, ask = FALSE)
  }
}
install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(apeglm)
  library(EnhancedVolcano)
})

set.seed(42)

# ---- 2. LOAD DATA ---------------------------------------------------------------

message(">> Reading counts and metadata...")

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
rownames(metadata) <- metadata[[sample_id_col]]

# Keep only samples present in both counts and metadata, in matching order
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No overlapping sample names between counts columns and metadata[[sample_id_col]].")
}
counts   <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

stopifnot(all(colnames(counts) == rownames(metadata)))

metadata[[condition_col]] <- factor(metadata[[condition_col]])
metadata[[condition_col]] <- relevel(metadata[[condition_col]], ref = ref_level)

message(sprintf("   %d genes x %d samples loaded.", nrow(counts), ncol(counts)))
print(table(metadata[[condition_col]]))

# ---- 3. BUILD DESeqDataSet & PRE-FILTER ------------------------------------------

message(">> Building DESeqDataSet...")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = as.formula(paste("~", condition_col))
)

# Pre-filter: keep genes with at least 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message(sprintf("   %d genes retained after low-count filtering.", nrow(dds)))

# ---- 4. RUN DESeq2 ---------------------------------------------------------------

message(">> Running DESeq2...")
dds <- DESeq(dds)

# Variance-stabilizing transform for visualization (PCA, heatmaps)
vsd <- vst(dds, blind = FALSE)

# ---- 5. EXTRACT RESULTS ---------------------------------------------------------

message(">> Extracting results...")

res <- results(dds,
                contrast = c(condition_col, contrast_group, ref_level),
                alpha = alpha_cutoff)

# Shrink log2FC for more reliable ranking/plotting (apeglm requires the coef name)
coef_name <- paste0(condition_col, "_", contrast_group, "_vs_", ref_level)
res_shrunk <- tryCatch(
  lfcShrink(dds, coef = coef_name, type = "apeglm"),
  error = function(e) {
    message("   apeglm shrinkage failed (", conditionMessage(e), "); using unshrunken results.")
    res
  }
)

res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

write.csv(res_df, file.path(out_dir, "DESeq2_results_full.csv"), row.names = FALSE)

sig_df <- subset(res_df, !is.na(padj) & padj < alpha_cutoff & abs(log2FoldChange) >= lfc_cutoff)
write.csv(sig_df, file.path(out_dir, "DESeq2_results_significant.csv"), row.names = FALSE)

message(sprintf("   %d significant genes (padj < %s, |log2FC| >= %s).",
                 nrow(sig_df), alpha_cutoff, lfc_cutoff))

# Normalized counts (for downstream use / GSEA / sharing)
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(out_dir, "normalized_counts.csv"))

# ---- 6. PCA PLOT ---------------------------------------------------------------

message(">> Generating PCA plot...")

pca_data <- plotPCA(vsd, intgroup = condition_col, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[condition_col]], label = name)) +
  geom_point(size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_bw(base_size = 13) +
  ggtitle("PCA of samples (VST counts)")

ggsave(file.path(out_dir, "PCA_plot.pdf"), p_pca, width = 7, height = 6)
ggsave(file.path(out_dir, "PCA_plot.png"), p_pca, width = 7, height = 6, dpi = 300)

# ---- 7. SAMPLE-DISTANCE HEATMAP -------------------------------------------------

message(">> Generating sample-distance heatmap...")

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(out_dir, "sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         annotation_col = metadata[, condition_col, drop = FALSE],
         main = "Sample-to-sample distance")
dev.off()

# ---- 8. MA PLOT ------------------------------------------------------------------

message(">> Generating MA plot...")

pdf(file.path(out_dir, "MA_plot.pdf"), width = 7, height = 6)
plotMA(res_shrunk, alpha = alpha_cutoff, ylim = c(-5, 5),
       main = paste0(contrast_group, " vs ", ref_level))
dev.off()

# ---- 9. VOLCANO PLOT --------------------------------------------------------------

message(">> Generating volcano plot...")

p_volcano <- EnhancedVolcano(res_df,
  lab = res_df$gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = alpha_cutoff,
  FCcutoff = lfc_cutoff,
  title = paste0(contrast_group, " vs ", ref_level),
  subtitle = NULL,
  caption = paste0("padj < ", alpha_cutoff, " & |log2FC| >= ", lfc_cutoff),
  labSize = 3.0,
  pointSize = 2.0,
  drawConnectors = TRUE,
  max.overlaps = 20
)

ggsave(file.path(out_dir, "volcano_plot.pdf"), p_volcano, width = 8, height = 8)
ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 8, height = 8, dpi = 300)

# ---- 10. TOP-GENE HEATMAP ---------------------------------------------------------

message(">> Generating top-gene expression heatmap...")

top_genes <- head(res_df$gene[!is.na(res_df$padj)], top_n_heatmap)

if (length(top_genes) > 1) {
  mat <- assay(vsd)[top_genes, , drop = FALSE]
  mat <- mat - rowMeans(mat)   # center per gene (z-score-like, row-centered)

  annotation_col <- metadata[, condition_col, drop = FALSE]

  pdf(file.path(out_dir, "top_genes_heatmap.pdf"), width = 8, height = 10)
  pheatmap(mat,
           annotation_col = annotation_col,
           show_rownames = TRUE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "none",
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
           main = paste0("Top ", length(top_genes), " DE genes (row-centered VST)"))
  dev.off()
} else {
  message("   Fewer than 2 genes with non-NA padj — skipping top-gene heatmap.")
}

# ---- 11. SUMMARY -------------------------------------------------------------------

sink(file.path(out_dir, "DESeq2_summary.txt"))
cat("DESeq2 pipeline summary\n")
cat("========================\n")
cat("Contrast:", contrast_group, "vs", ref_level, "\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("Significant (padj <", alpha_cutoff, ", |log2FC| >=", lfc_cutoff, "):", nrow(sig_df), "\n\n")
summary(res_shrunk)
sink()

message(">> Done. All outputs written to: ", normalizePath(out_dir))
message("   - DESeq2_results_full.csv / DESeq2_results_significant.csv")
message("   - normalized_counts.csv")
message("   - PCA_plot.pdf/png, sample_distance_heatmap.pdf")
message("   - MA_plot.pdf, volcano_plot.pdf/png, top_genes_heatmap.pdf")
message("   - DESeq2_summary.txt")
