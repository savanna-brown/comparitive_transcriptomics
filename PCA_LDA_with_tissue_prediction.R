# ===============================================
# PCA -> LDA : TISSUE PREDICTION
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

species_levels <- c("species1", "species2", "species3", "species4", "species5")
tissue_levels <- c("tissue1", "tissue2", "tissue3", "tissue4",
                   "tissue5", "tissue6", "tissue7", "tissue8")
prediction_tissue <- "tissue1"

max_pcs_to_test <- XX

top_pcs <- XX

# ===============================================
# LOAD AND FILTER DATA
# ===============================================

merged_mat <- read.table(merged_matrix_file, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE)
col_data_all <- read.table(sample_info_file, header = TRUE, row.names = 1,
                           sep = "\t", check.names = FALSE)

col_data <- col_data_all[
  col_data_all$species %in% species_levels &
  col_data_all$tissue %in% tissue_levels, , drop = FALSE
]
valid_samples <- intersect(rownames(col_data), colnames(merged_mat))
col_data <- col_data[valid_samples, , drop = FALSE]
mat <- merged_mat[, valid_samples, drop = FALSE]

col_data$species <- factor(col_data$species, levels = species_levels)
col_data$tissue <- factor(col_data$tissue, levels = tissue_levels)

if (!prediction_tissue %in% tissue_levels) stop("prediction_tissue is not in tissue_levels")
train_samples <- rownames(col_data)[col_data$tissue != prediction_tissue]
prediction_samples <- rownames(col_data)[col_data$tissue == prediction_tissue]
if (!length(prediction_samples)) stop("No samples found for prediction_tissue")
if (length(unique(col_data[train_samples, "tissue"])) < 3) {
  stop("At least three training tissue categories are needed for LD1 vs LD2 plotting")
}

#vst normalization
dds <- DESeqDataSetFromMatrix(round(mat), col_data, design = ~ tissue + species)
dds <- DESeq(dds)
vst_counts <- assay(vst(dds, blind = FALSE))

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

# ===============================================
# PCA ON TRAINING SAMPLES
# ===============================================

train_mat <- vst_counts[, train_samples, drop = FALSE]
pca_train <- prcomp(t(train_mat), scale. = TRUE)


# ===============================================
# PCA-LDA LEAVE-ONE-OUT CROSS-VALIDATION
# ===============================================

expression_data <- t(train_mat)
tissue_labels <- droplevels(col_data[rownames(expression_data), "tissue"])
pc_range <- 2:max_pcs_to_test

predictions <- matrix(NA_character_, nrow = nrow(expression_data),
                      ncol = length(pc_range),
                      dimnames = list(rownames(expression_data), pc_range))

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
  tibble(nPC = n,
         Accuracy = if (anyNA(pred)) NA_real_
         else mean(pred == as.character(tissue_labels)))
})

cv_results <- filter(cv_results, !is.na(Accuracy))
if (!nrow(cv_results)) stop("No PC counts could be evaluated by LOOCV")

write.csv(cv_results,
          file.path(output_dir, paste0("LDA_CV_results_", set_name, ".csv")),
          row.names = FALSE)

p_cv <- ggplot(cv_results, aes(nPC, Accuracy)) +
  geom_line() + geom_point() + theme_minimal(base_size = 12) +
  labs(title = "PCA-LDA LOOCV Accuracy", x = "Number of PCs",
       y = "Classification accuracy")

ggsave(file.path(output_dir, paste0("LDA_CV_accuracy_vs_nPC_", set_name, ".png")),
       p_cv, width = 8, height = 6, dpi = 300)

# ===============================================
# LDA AND HELD-OUT TISSUE PREDICTION
# ===============================================

cat(paste0("Fitting LDA with ", top_pcs, " PCs...\n"))
train_df <- as.data.frame(pca_train$x[, 1:top_pcs, drop = FALSE])
train_df$tissue <- droplevels(col_data[rownames(train_df), "tissue"])
lda_model <- MASS::lda(tissue ~ ., data = train_df)

prediction_pcs <- as.data.frame(predict(
  pca_train, newdata = t(vst_counts[, prediction_samples, drop = FALSE])
)[, 1:top_pcs, drop = FALSE])

lda_pred <- predict(lda_model, newdata = prediction_pcs)
prediction_df <- as.data.frame(lda_pred$posterior) %>%
  rownames_to_column("Sample") %>%
  mutate(PredictedTissue = as.character(lda_pred$class), .after = Sample)

write.csv(prediction_df,
          file.path(output_dir, paste0("LDA_tissue_predictions_", set_name, ".csv")),
          row.names = FALSE)

# ===============================================
# SINGLE-SAMPLE JACKKNIFE PREDICTION STABILITY
# ===============================================

jackknife <- lapply(seq_len(nrow(train_df)), function(i) {
  model <- MASS::lda(tissue ~ ., data = train_df[-i, , drop = FALSE])
  predict(model, newdata = prediction_pcs)
})
names(jackknife) <- rownames(train_df)

jackknife_posteriors <- map_df(names(jackknife), function(i) {
  as.data.frame(jackknife[[i]]$posterior) %>%
    rownames_to_column("Sample") %>%
    mutate(LeftOutSample = i, .before = Sample)
})

write.csv(
  jackknife_posteriors,
  file.path(output_dir, paste0("LOO_single_raw_posteriors_", set_name, ".csv")),
  row.names = FALSE
)

jackknife_predictions <- sapply(jackknife, function(x) as.character(x$class))
rownames(jackknife_predictions) <- prediction_samples

write.csv(jackknife_predictions,
          file.path(output_dir, paste0("LOO_single_pred_matrix_", set_name, ".csv")))

jackknife_summary <- map_df(prediction_samples, function(s) {
  classes <- jackknife_predictions[s, ]
  map_df(levels(train_df$tissue), function(t) {
    probabilities <- sapply(jackknife, function(x) x$posterior[s, t])
    tibble(Sample = s, Tissue = t, Percent = 100 * mean(classes == t),
           MeanPosterior = mean(probabilities),
           Q025 = unname(quantile(probabilities, 0.025)),
           Q975 = unname(quantile(probabilities, 0.975)))
  })
})

write.csv(jackknife_summary,
          file.path(output_dir, paste0("LOO_single_summary_", set_name, ".csv")),
          row.names = FALSE)

# ===============================================
# LDA PROJECTION AND SUMMARY PLOTS
# ===============================================

train_coords <- as.data.frame(predict(lda_model)$x) %>%
  mutate(Sample = rownames(train_df),
         Species = col_data[Sample, "species"],
         Tissue = as.character(col_data[Sample, "tissue"]))

prediction_coords <- as.data.frame(lda_pred$x) %>%
  mutate(Sample = prediction_samples,
         Species = col_data[Sample, "species"],
         Tissue = prediction_tissue)

lda_var <- lda_model$svd^2 / sum(lda_model$svd^2) * 100

p_projection <- ggplot(train_coords, aes(LD1, LD2, color = Tissue, shape = Species)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_point(data = prediction_coords, inherit.aes = FALSE,
             aes(LD1, LD2), shape = 8, color = "black", size = 3) +
  scale_color_manual(values = tissue_colors, drop = FALSE) +
  scale_shape_manual(values = species_shapes, drop = FALSE) +
  labs(title = paste0("LDA projection: ", prediction_tissue),
       x = paste0("LD1 (", round(lda_var[1], 1), "% discriminant variance)"),
       y = paste0("LD2 (", round(lda_var[2], 1), "% discriminant variance)")) +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, paste0("LDA_tissue_projection_", set_name, ".png")),
       p_projection, width = 8, height = 6, dpi = 300)

p_heatmap <- ggplot(jackknife_summary,
                    aes(Tissue, Sample, fill = MeanPosterior)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  labs(title = "Jackknife mean posterior probabilities") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, paste0("LOO_single_summary_heatmap_", set_name, ".png")),
       p_heatmap, width = 8, height = 6, dpi = 300)

saveRDS(list(pca_train = pca_train, cv_results = cv_results,
             lda_model = lda_model, top_pcs = top_pcs,
             prediction_df = prediction_df, jackknife_summary = jackknife_summary,
             vst_counts = vst_counts, col_data = col_data),
        file.path(output_dir, paste0("PCA_CV_LDA_prediction_state_", set_name, ".rds")))

cat("PCA-LDA prediction results saved at: ", output_dir, "\n")
