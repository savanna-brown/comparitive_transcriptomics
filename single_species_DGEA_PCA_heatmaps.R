# ===============================================
# SINGLE-SPECIES DGEA + PCA + HEATMAPS (ONE STAGE)
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

# Input count matrix should contain samples from one developmental stage only.
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

# Analyze each species separately (no stage grouping).
sample_groups <- split(rownames(col_data),
                       factor(col_data$species, levels = species_levels),
                       drop = TRUE)

# ===============================================
# COLORS
# ===============================================
tissue_colors <- c(
  "tissue1" = "#E41A1C", "tissue2" = "#FF7F00",
  "tissue3" = "#FFD92F", "tissue4" = "#4DAF4A",
  "tissue5" = "#00BFC4", "tissue6" = "#377EB8",
  "tissue7" = "#984EA3", "tissue8" = "#F781BF"
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

plot_pca <- function(x, label, out_dir, group_label, sample_meta) {
  pca <- pca_of_matrix(x)
  scores <- as.data.frame(pca$x)
  scores$Tissue <- factor(sample_meta[rownames(scores), "tissue"],
                          levels = tissue_levels)
  fractions <- summary(pca)$importance[2, ] * 100

  for (i in 1:4) {
    if (ncol(pca$x) < i + 1) next
    xpc <- paste0("PC", i)
    ypc <- paste0("PC", i + 1)
    p <- ggplot(scores, aes(x = .data[[xpc]], y = .data[[ypc]], color = Tissue)) +
      geom_point(size = 3) +
      scale_color_manual(values = tissue_colors, drop = FALSE) +
      theme_bw(base_size = 14) +
      labs(title = paste(group_label, label, xpc, "vs", ypc),
           x = sprintf("%s (%.1f%%)", xpc, fractions[i]),
           y = sprintf("%s (%.1f%%)", ypc, fractions[i + 1]))
    ggsave(file.path(out_dir, paste0(label, "_", xpc, "_", ypc, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }

  scree <- factoextra::fviz_eig(pca, main = paste(group_label, label, "Scree plot"))
  ggsave(file.path(out_dir, paste0("Scree_", label, ".png")),
         scree, width = 8, height = 6, dpi = 300)
  pca
}

plot_heatmap <- function(x, label, out_dir, group_label, sample_meta) {
  annotation_df <- data.frame(
    tissue = factor(sample_meta[colnames(x), "tissue"], levels = tissue_levels),
    row.names = colnames(x)
  )
  present_tissues <- unique(as.character(annotation_df$tissue))
  annotation_colors <- list(tissue = tissue_colors[present_tissues])
  hm <- pheatmap::pheatmap(
    x, cluster_rows = TRUE, cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = color_scale, breaks = breaks_vst,
    main = paste(group_label, label), silent = TRUE
  )
  ggsave(file.path(out_dir, paste0(label, ".png")),
         plot = hm$gtable, width = 10, height = 12, dpi = 300)
}

# ===============================================
# SINGLE-SPECIES DESEQ2 ANALYSIS
# ===============================================
for (group_label in names(sample_groups)) {
  samples <- sample_groups[[group_label]]
  sample_meta <- col_data[samples, , drop = FALSE]
  sample_meta$tissue <- droplevels(factor(sample_meta$tissue,
                                          levels = tissue_levels))
  if (nlevels(sample_meta$tissue) < 2) {
    warning("Skipping ", group_label, ": fewer than two tissues.")
    next
  }

  group_mat <- round(mat[, samples, drop = FALSE])
  group_mat <- group_mat[rowSums(group_mat) > 0, , drop = FALSE]
  out_dir <- file.path(output_dir, group_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nProcessing ", group_label, " (", ncol(group_mat), " samples)...\n", sep = "")

  
  dds <- DESeqDataSetFromMatrix(countData = group_mat,
                                colData = sample_meta, design = ~ tissue)
  dds <- DESeq(dds)

  # ===============================================
  # UNION OF SIGNIFICANT GENES ACROSS PAIRWISE CONTRASTS
  # ===============================================
  tissue_names <- levels(sample_meta$tissue)
  tissue_pairs <- combn(tissue_names, 2, simplify = FALSE)
  sig_genes <- character()

  for (pair in tissue_pairs) {
    res <- results(dds, contrast = c("tissue", pair[1], pair[2]))
    res_sig <- res[!is.na(res$padj), ]
    sig_genes <- union(sig_genes,
                       rownames(res_sig[res_sig$padj <= padj_cutoff, ]))
  }

  if (!length(sig_genes)) {
    warning("No significant DE genes in ", group_label, "; skipping PCA/heatmaps.")
    next
  }

  # ===============================================
  # VST + BATCH CORRECTION (FOR VISUALIZATION ONLY)
  # ===============================================
  vst_counts <- assay(vst(dds, blind = FALSE))
  vst_uncorrected <- vst_counts[sig_genes, , drop = FALSE]

  has_pool <- "pool" %in% colnames(sample_meta) &&
              length(unique(stats::na.omit(sample_meta$pool))) > 1
  if (has_pool) {
    corrected <- limma::removeBatchEffect(
      vst_uncorrected, batch = factor(sample_meta$pool)
    )
  } else {
    corrected <- vst_uncorrected
    message("No variable 'pool' column in ", group_label,
            "; batch correction skipped.")
  }

  vst_centered_uncorrected <- center_genes(vst_uncorrected)
  vst_centered_limma <- center_genes(corrected)

  # ===============================================
  # PCA PLOTS AND HEATMAPS
  # ===============================================
  pca_uncorrected <- plot_pca(vst_centered_uncorrected, "Uncorrected_PCA",
                               out_dir, group_label, sample_meta)
  pca_corrected <- plot_pca(vst_centered_limma, "Batch_corrected_PCA",
                             out_dir, group_label, sample_meta)
  plot_heatmap(vst_centered_uncorrected, "Heatmap_uncorrected",
               out_dir, group_label, sample_meta)
  plot_heatmap(vst_centered_limma, "Heatmap_batch_corrected",
               out_dir, group_label, sample_meta)

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
      col_data = sample_meta
    ),
    file.path(out_dir, paste0("PCA_heatmap_state_", group_label, ".rds"))
  )
  cat("PCA and heatmap results saved at: ", out_dir, "\n")
}
