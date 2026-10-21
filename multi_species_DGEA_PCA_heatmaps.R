# ===============================================
# MULTI-SPECIES DGEA + PCA + HEATMAPS (ONE STAGE)
# ===============================================

library(tidyverse)
library(DESeq2)
library(limma)
library(pheatmap)
library(factoextra)

# ===============================================
# PARAMETERS
# ===============================================
set_name <- "example_dataset"
merged_matrix_file <- "path/to/multi_species_matrix_file.tsv"
sample_info_file <- "path/to/sample_metadata.tsv"
output_dir <- "path/to/output_directory"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The input count matrix must contain one developmental stage only.
# Rows represent shared one-to-one orthologues; columns represent samples.
species_levels <- c("species1", "species2", "species3", "species4", "species5")
tissue_levels <- c("tissue1", "tissue2", "tissue3", "tissue4",
                   "tissue5", "tissue6", "tissue7", "tissue8")

padj_cutoff <- 0.05

# ===============================================
# LOAD AND FILTER DATA
# ===============================================
cat("Loading raw expression matrix and metadata...\n")
merged_mat <- read.table(merged_matrix_file, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE)
col_data_all <- read.table(sample_info_file, header = TRUE, row.names = 1,
                           sep = "\t", check.names = FALSE)

if (!all(c("species", "tissue") %in% colnames(col_data_all))) {
  stop("Metadata must contain 'species' and 'tissue' columns.")
}

col_data <- col_data_all[
  col_data_all$species %in% species_levels &
    col_data_all$tissue %in% tissue_levels, , drop = FALSE
]
valid_samples <- intersect(rownames(col_data), colnames(merged_mat))
if (!length(valid_samples)) stop("No matching samples between matrix and metadata.")
col_data <- col_data[valid_samples, , drop = FALSE]
mat <- as.matrix(merged_mat[, valid_samples, drop = FALSE])
storage.mode(mat) <- "numeric"
if (anyNA(mat) || any(mat < 0)) {
  stop("Count matrix must contain nonnegative counts without NA values.")
}

col_data$species <- droplevels(factor(col_data$species, levels = species_levels))
col_data$tissue <- droplevels(factor(col_data$tissue, levels = tissue_levels))
if (nlevels(col_data$species) < 2) stop("At least two species are required.")
if (nlevels(col_data$tissue) < 2) stop("At least two tissues are required.")

mat <- round(mat)
mat <- mat[rowSums(mat) > 0, , drop = FALSE]

# ===============================================
# COLORS AND SHAPES
# ===============================================
tissue_colors <- c(
  "tissue1" = "#E41A1C", "tissue2" = "#FF7F00",
  "tissue3" = "#FFD92F", "tissue4" = "#4DAF4A",
  "tissue5" = "#00BFC4", "tissue6" = "#377EB8",
  "tissue7" = "#984EA3", "tissue8" = "#F781BF"
)
species_shapes <- c(
  "species1" = 3, "species2" = 4, "species3" = 15,
  "species4" = 16, "species5" = 17
)
color_scale <- colorRampPalette(c("gray0", "gray", "ivory1"))(100)
breaks_vst <- seq(-2, 2, length.out = 101)

# ===============================================
# HELPERS
# ===============================================
center_genes <- function(x) sweep(x, 1, rowMeans(x), "-")

pca_of_matrix <- function(x) {
  keep <- apply(x, 1, sd) > 0
  if (sum(keep) < 2) stop("Too few variable genes for PCA.")
  prcomp(t(x[keep, , drop = FALSE]), scale. = TRUE)
}

plot_pca <- function(x, label, output_dir, sample_meta) {
  pca <- pca_of_matrix(x)
  scores <- as.data.frame(pca$x)
  scores$Tissue <- factor(sample_meta[rownames(scores), "tissue"],
                          levels = tissue_levels)
  scores$Species <- factor(sample_meta[rownames(scores), "species"],
                           levels = species_levels)
  fractions <- summary(pca)$importance[2, ] * 100

  for (i in 1:4) {
    if (ncol(pca$x) < i + 1) next
    xpc <- paste0("PC", i)
    ypc <- paste0("PC", i + 1)
    p <- ggplot(scores, aes(x = .data[[xpc]], y = .data[[ypc]],
                            color = Tissue, shape = Species)) +
      geom_point(size = 3) +
      scale_color_manual(values = tissue_colors, drop = FALSE) +
      scale_shape_manual(values = species_shapes, drop = FALSE) +
      theme_bw(base_size = 14) +
      labs(title = paste("Multi-species", label, xpc, "vs", ypc),
           x = sprintf("%s (%.1f%%)", xpc, fractions[i]),
           y = sprintf("%s (%.1f%%)", ypc, fractions[i + 1]))
    ggsave(file.path(output_dir, paste0(label, "_", xpc, "_", ypc, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }

  scree <- factoextra::fviz_eig(pca, main = paste("Multi-species", label, "Scree plot"))
  ggsave(file.path(output_dir, paste0("Scree_", label, ".png")),
         scree, width = 8, height = 6, dpi = 300)
  pca
}

plot_heatmap <- function(x, label, output_dir, sample_meta) {
  annotation_df <- data.frame(
    tissue = factor(sample_meta[colnames(x), "tissue"], levels = tissue_levels),
    species = factor(sample_meta[colnames(x), "species"], levels = species_levels),
    row.names = colnames(x)
  )
  # Add a species annotation with a distinct color for each species.
  species_colors <- setNames(grDevices::hcl.colors(length(species_levels), "Dark 3"),
                             species_levels)
  hm <- pheatmap::pheatmap(
    x, cluster_rows = TRUE, cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_col = annotation_df,
    annotation_colors = list(tissue = tissue_colors, species = species_colors),
    color = color_scale, breaks = breaks_vst,
    main = paste("Multi-species", label), silent = TRUE
  )
  ggsave(file.path(output_dir, paste0(label, ".png")),
         plot = hm$gtable, width = 10, height = 12, dpi = 300)
}

# ===============================================
# JOINT MULTI-SPECIES DESEQ2 ANALYSIS
# ===============================================

dds <- DESeqDataSetFromMatrix(countData = mat, colData = col_data,
                              design = ~ tissue + species)
dds <- DESeq(dds)

# ===============================================
# UNION OF SIGNIFICANT GENES ACROSS PAIRWISE CONTRASTS
# ===============================================
tissue_names <- levels(col_data$tissue)
tissue_pairs <- combn(tissue_names, 2, simplify = FALSE)
sig_genes <- character()

for (pair in tissue_pairs) {
  res <- results(dds, contrast = c("tissue", pair[1], pair[2]))
  res_sig <- res[!is.na(res$padj), ]
  sig_genes <- union(sig_genes,
                     rownames(res_sig[res_sig$padj <= padj_cutoff, ]))
}

if (!length(sig_genes)) stop("No significant DE genes; PCA/heatmaps cannot be produced.")

# ===============================================
# VST + BATCH CORRECTION (FOR VISUALIZATION ONLY)
# ===============================================
vst_counts <- assay(vst(dds, blind = FALSE))
vst_uncorrected <- vst_counts[sig_genes, , drop = FALSE]

has_pool <- "pool" %in% colnames(col_data) &&
            length(unique(stats::na.omit(col_data$pool))) > 1
if (has_pool) {
  if (anyNA(col_data$pool)) stop("Missing pool assignments cannot be batch corrected.")
  corrected <- limma::removeBatchEffect(
    vst_uncorrected, batch = factor(col_data$pool)
  )
} else {
  corrected <- vst_uncorrected
  message("No variable 'pool' column; batch correction skipped.")
}

vst_centered_uncorrected <- center_genes(vst_uncorrected)
vst_centered_limma <- center_genes(corrected)

# ===============================================
# PCA PLOTS AND HEATMAPS
# ===============================================
pca_uncorrected <- plot_pca(vst_centered_uncorrected, "Uncorrected_PCA",
                            output_dir, col_data)
pca_corrected <- plot_pca(vst_centered_limma, "Batch_corrected_PCA",
                          output_dir, col_data)
plot_heatmap(vst_centered_uncorrected, "Heatmap_uncorrected",
             output_dir, col_data)
plot_heatmap(vst_centered_limma, "Heatmap_batch_corrected",
             output_dir, col_data)

# ===============================================
# SAVE PCA AND HEATMAP DATA AND RESULTS
# ===============================================
saveRDS(
  list(
    dds = dds,
    sig_genes = sig_genes,
    vst_counts = vst_counts,
    vst_centered_uncorrected = vst_centered_uncorrected,
    vst_centered_limma = vst_centered_limma,
    pca_uncorrected = pca_uncorrected,
    pca_corrected = pca_corrected,
    col_data = col_data
  ),
  file.path(output_dir, paste0("PCA_heatmap_state_", set_name, ".rds"))
)

cat("Multi-species PCA and heatmap results saved at: ", output_dir, "\n")
