# ===============================================
# PCA -> LDA : ALL SAMPLES
# ===============================================

library(tidyverse)
library(DESeq2)
library(MASS)

# ===============================================
# PARAMETERS
# ===============================================

set_name <- "example_dataset"

merged_matrix_file <- "path/to/multi_species_matrix_file.tsv"

sample_info_file <- "path/to/sample_metadata.tsv"

output_dir <- "path/to/output_directory"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

max_pcs_to_test <- XX

top_pcs <- XX

# ===============================================
# LOAD AND FILTER DATA
# ===============================================

merged_mat <- read.table(
  merged_matrix_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

col_data_all <- read.table(
  sample_info_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

#species and tissues to include in analysis, in order to appear on figs
species_levels <- c("species1", "species2", "species3", "species4", "species5")
tissue_levels <- c("tissue1", "tissue2", "tissue3", "tissue4", "tissue5", "tissue6", "tissue7", "tissue8")


col_data <- col_data_all[
  col_data_all$tissue %in% tissue_levels &
    col_data_all$species %in% species_levels, ]

valid_samples <- intersect(rownames(col_data), colnames(merged_mat))

col_data <- col_data[valid_samples, , drop = FALSE]

mat <- merged_mat[, valid_samples, drop = FALSE]


# Final alignment
samples_in_matrix <- intersect(colnames(mat), rownames(col_data))
col_data <- col_data[samples_in_matrix, , drop = FALSE]
mat <- mat[, samples_in_matrix, drop = FALSE]

col_data$species <- factor(col_data$species, levels = species_levels)
col_data$tissue <- factor(col_data$tissue, levels = tissue_levels)

#vst normalization
dds <- DESeqDataSetFromMatrix(round(mat), col_data, design = ~ tissue + species)
dds <- DESeq(dds)
vst_counts <- assay(vst(dds, blind = FALSE))

# ===============================================
# COLORS AND SHAPES
# ===============================================

tissue_colors <- c(
  "tissue1" = "#E41A1C",  # Red
  "tissue2" = "#FF7F00",  # Orange
  "tissue3" = "#FFD92F",  # Yellow
  "tissue4" = "#4DAF4A",  # Green
  "tissue5" = "#00BFC4",  # Cyan
  "tissue6" = "#377EB8",  # Blue
  "tissue7" = "#984EA3",  # Purple
  "tissue8" = "#F781BF"   # Pink
)


species_shapes <- c(
  "species1" = 3,
  "species2" = 4,
  "species3" = 15,
  "species4" = 16,
  "species5" = 17
)

# ===============================================
# PCA ON ALL SAMPLES
# ===============================================

train_samples <- rownames(col_data)
train_mat <- vst_counts[, train_samples, drop = FALSE]
pca_train <- prcomp(t(train_mat), scale. = TRUE)

# PCA visualization
pca_df <- as.data.frame(pca_train$x)
pca_df$Species <- factor(col_data[rownames(pca_df), "species"], levels = species_levels)
pca_df$Tissue <- factor(col_data[rownames(pca_df), "tissue"], levels = tissue_levels)
pca_var <- summary(pca_train)$importance[2, ] * 100

plot_pca <- function(xpc, ypc) {
  x_idx <- as.numeric(sub("PC", "", xpc))
  y_idx <- as.numeric(sub("PC", "", ypc))
  
  ggplot(pca_df, aes_string(x = xpc, y = ypc, color = "Tissue", shape = "Species")) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = tissue_colors, drop = FALSE) +
    scale_shape_manual(values = species_shapes, drop = FALSE) +
    labs(
      title = paste0("PCA: all samples (", xpc, " vs ", ypc, ")"),
      x = paste0(xpc, " (", round(pca_var[x_idx], 1), "% variance)"),
      y = paste0(ypc, " (", round(pca_var[y_idx], 1), "% variance)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

ggsave(
  file.path(output_dir, paste0("PCA_PC1_PC2_", set_name, ".png")),
  plot = plot_pca("PC1", "PC2"),
  width = 8, height = 6, dpi = 300
)

if (ncol(pca_train$x) >= 4) {
  ggsave(
    file.path(output_dir, paste0("PCA_PC3_PC4_", set_name, ".png")),
    plot = plot_pca("PC3", "PC4"),
    width = 8, height = 6, dpi = 300
  )
}

# ===============================================
# PCA-LDA LEAVE-ONE-OUT CROSS-VALIDATION
# ===============================================

max_pcs_to_test <- 50

cat("Running PCA-LDA LOOCV...\n")

expression_data <- t(train_mat)
tissue_labels <- droplevels(col_data[rownames(expression_data), "tissue"])
pc_range <- 2:max_pcs_to_test

predictions <- matrix(
  NA_character_,
  nrow = nrow(expression_data),
  ncol = length(pc_range),
  dimnames = list(rownames(expression_data), as.character(pc_range))
)

for (i in seq_len(nrow(expression_data))) {
  
  train_x <- expression_data[-i, , drop = FALSE]
  test_x <- expression_data[i, , drop = FALSE]
  train_y <- droplevels(tissue_labels[-i])
  
  if (nlevels(train_y) < nlevels(tissue_labels)) next
  
  keep_genes <- apply(train_x, 2, sd) > 0
  train_x <- train_x[, keep_genes, drop = FALSE]
  test_x <- test_x[, keep_genes, drop = FALSE]
  
  pca_fold <- prcomp(train_x, center = TRUE, scale. = TRUE)
  test_pcs <- predict(pca_fold, newdata = test_x)
  
  for (n in pc_range) {
    
    if (n > ncol(pca_fold$x) ||
        n >= nrow(train_x) - nlevels(train_y)) next
    
    train_pcs <- pca_fold$x[, 1:n, drop = FALSE]
    test_pc <- test_pcs[, 1:n, drop = FALSE]
    
    predictions[i, as.character(n)] <- tryCatch({
      model <- MASS::lda(x = train_pcs, grouping = train_y)
      as.character(predict(model, newdata = test_pc)$class)
    }, error = function(e) NA_character_)
  }
  
  if (i %% 10 == 0) {
    cat("Completed", i, "of", nrow(expression_data), "samples\n")
  }
}

cv_results <- map_df(pc_range, function(n) {
  pred <- predictions[, as.character(n)]
  
  tibble(
    nPC = n,
    Accuracy = if (anyNA(pred)) NA_real_
    else mean(pred == as.character(tissue_labels))
  )
})

cv_results <- filter(cv_results, !is.na(Accuracy))
if (!nrow(cv_results)) stop("No PC counts could be evaluated by LOOCV")

write.csv(
  cv_results,
  file.path(output_dir, paste0("LDA_CV_results_", set_name, ".csv")),
  row.names = FALSE
)

p_cv <- ggplot(cv_results, aes(nPC, Accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCA-LDA LOOCV Accuracy",
    x = "Number of PCs",
    y = "Classification accuracy"
  )

ggsave(
  file.path(output_dir, paste0("LDA_CV_accuracy_vs_nPC_", set_name, ".png")),
  p_cv,
  width = 8, height = 6, dpi = 300
)

# ===============================================
# LDA ON PCs OF ALL SAMPLES
# ===============================================

if (top_pcs > ncol(pca_train$x)) {
  stop("top_pcs exceeds available PCs.")
}

train_df <- as.data.frame(pca_train$x[, 1:top_pcs, drop = FALSE])
train_df$tissue <- col_data$tissue

lda_model <- lda(tissue ~ ., data = train_df)


# LDA visualization

lda_pred <- predict(lda_model)

lda_df <- as.data.frame(lda_pred$x)
lda_df$Species <- factor(col_data[rownames(lda_df), "species"],
                         levels = species_levels)
lda_df$Tissue <- factor(col_data[rownames(lda_df), "tissue"],
                        levels = tissue_levels)

lda_var <- (lda_model$svd^2) / sum(lda_model$svd^2) * 100

plot_lda <- function(xld, yld) {
  
  x_idx <- as.numeric(sub("LD", "", xld))
  y_idx <- as.numeric(sub("LD", "", yld))
  
  ggplot(
    lda_df,
    aes_string(x = xld, y = yld,
               color = "Tissue",
               shape = "Species")
  ) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = tissue_colors, drop = FALSE) +
    scale_shape_manual(values = species_shapes, drop = FALSE) +
    labs(
      title = paste0("LDA on PCs: all samples (", xld, " vs ", yld, ")"),
      subtitle = paste0("trained on first ", top_pcs, " PCs"),
      x = paste0(xld, " (", round(lda_var[x_idx], 1), "% discriminant variance)"),
      y = paste0(yld, " (", round(lda_var[y_idx], 1), "% discriminant variance)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

ggsave(
  file.path(output_dir, paste0("LDA_LD1_LD2_", set_name, ".png")),
  plot = plot_lda("LD1", "LD2"),
  width = 8, height = 6, dpi = 300
)

if (ncol(lda_df) >= 4) {
  ggsave(
    file.path(output_dir, paste0("LDA_LD3_LD4_", set_name, ".png")),
    plot = plot_lda("LD3", "LD4"),
    width = 8, height = 6, dpi = 300
  )
}

# save PCA, CV, and LDA data and results
saveRDS(
  list(
    pca_train = pca_train,
    pca_df = pca_df,
    cv_results = cv_results,
    lda_model = lda_model,
    top_pcs = top_pcs,
    vst_counts = vst_counts,
    col_data = col_data
  ),
  file.path(output_dir, paste0("PCA_CV_LDA_state_", set_name, ".rds"))
)

cat("PCA, CV, and LDA data and results saved at: ", output_dir, "\n")
