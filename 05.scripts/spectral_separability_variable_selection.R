#!/usr/bin/env Rscript
# =============================================================================
# Spectral separability and variable selection for the Lake Chilwa 4-class model
# Built 2026-07-05. Mirrors the methodology in Murphy et al. (2026, Forest
# Ecology and Management 618:123985) and the darkwoods_beetles / darkwoods_seedlings
# analyses: separability index + non-parametric tests, PCA, and glmnet importance.
#
# PURPOSE
#   Decide, empirically, which spectral bands and indices best discriminate the
#   four land-cover classes (open water, emergent/flooded vegetation, dry
#   vegetation, bare soil / exposed lakebed), rather than assuming a feature set.
#   Starts with four classes; a fifth (salt-crust vs sandy lakebed) can be added
#   and re-screened.
#
# WORKFLOW (three stages, per the FEM paper)
#   1. Separability index analysis  : M = |mu1 - mu2| / (sd1 + sd2), pairwise.
#   2. Distributional / ANOVA test   : Shapiro-Wilk normality, then
#      Kruskal-Wallis (multi-class) and Wilcoxon rank-sum (pairwise), because
#      spectral variables are typically non-normal.
#   3. Multivariate selection        : PCA (prcomp) for structure and
#      dimensionality, and glmnet multinomial LASSO/ridge for variable importance.
#
# INPUT
#   A data.frame `spec_df` with one row per training pixel: numeric columns for
#   each candidate feature and a factor column `class` with the four class labels.
#   A GEE extraction preamble is provided; or load a saved CSV.
#
# RUN: Rscript 05.scripts/spectral_separability_variable_selection.R
#   The GEE preamble needs an authenticated Earth Engine session; the statistics
#   run locally on the extracted data.frame. Cannot execute here; you run it.
# =============================================================================

suppressMessages({
  library(reticulate); library(dplyr); library(tidyr)
  library(glmnet); library(ggplot2)
})

FEATURES <- c("blue","green","red","nir","swir1","swir2",
              "NDWI","MNDWI","AWEIsh","WRI","NDPI","NDVI")
CLASSES  <- c("water","flooded_veg","dry_veg","bare_soil")

# ------------------- 1. EXTRACT TRAINING SPECTRA (GEE) -----------------------
# Assumes `ref_all` (reference composite with the 12 feature bands) and
# `training` (labelled FeatureCollection with a 'class' property, 4 classes)
# already exist from the manuscript's Methods. Pull the sampled spectra to R.
extract_spec_df <- function(ref_all, training, ee, scale = 30L) {
  samp <- ref_all$select(FEATURES)$sampleRegions(
    collection = training, properties = list("class"), scale = scale,
    geometries = FALSE)
  info <- samp$getInfo()
  df <- do.call(rbind, lapply(info$features, function(f) as.data.frame(f$properties)))
  df$class <- factor(df$class, labels = CLASSES[sort(unique(df$class))])
  df
}
# Example (uncomment in the notebook, after Methods objects exist):
# ee <- reticulate::import("ee"); ee$Initialize(project = Sys.getenv("EE_PROJECT","murphys-deforisk"))
# spec_df <- extract_spec_df(ref_all, training, ee)
# write.csv(spec_df, "03.outputs/training_spectra.csv", row.names = FALSE)
#
# Or load a saved extraction:
# spec_df <- read.csv("03.outputs/training_spectra.csv"); spec_df$class <- factor(spec_df$class)

if (!exists("spec_df")) stop("Provide spec_df (training spectra with a 'class' factor).")

# ------------------- 2. DESCRIPTIVES + NORMALITY -----------------------------
describe_features <- function(df) {
  df %>% pivot_longer(all_of(FEATURES), names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    summarise(n = n(), mean = mean(value), sd = sd(value), median = median(value),
              se = sd(value)/sqrt(n()),
              shapiro_W = tryCatch(shapiro.test(sample(value, min(5000, length(value))))$statistic, error = function(e) NA),
              shapiro_p = tryCatch(shapiro.test(sample(value, min(5000, length(value))))$p.value, error = function(e) NA),
              .groups = "drop")
}
desc <- describe_features(spec_df)
cat("== Feature descriptives and Shapiro-Wilk normality ==\n"); print(desc, n = Inf)
cat("\nNote: low Shapiro p indicates non-normality; justifies the non-parametric\n",
    "Kruskal-Wallis / Wilcoxon tests below (as in Murphy et al. 2026).\n")

# ------------------- 3. SEPARABILITY INDEX M (pairwise) ----------------------
# M = |mu1 - mu2| / (sd1 + sd2). M > 1 indicates good separability; higher is better.
sep_index <- function(df, feature) {
  cl <- levels(df$class); out <- list()
  for (i in 1:(length(cl)-1)) for (j in (i+1):length(cl)) {
    a <- df[[feature]][df$class == cl[i]]; b <- df[[feature]][df$class == cl[j]]
    M <- abs(mean(a) - mean(b)) / (sd(a) + sd(b))
    out[[length(out)+1]] <- data.frame(feature = feature,
      pair = paste(cl[i], "vs", cl[j]), M = round(M, 3))
  }
  do.call(rbind, out)
}
sep_tbl <- do.call(rbind, lapply(FEATURES, function(f) sep_index(spec_df, f)))
sep_wide <- sep_tbl %>% pivot_wider(names_from = pair, values_from = M) %>%
  mutate(M_mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(M_mean))
cat("\n== Separability index M by feature (mean across class pairs, ranked) ==\n")
print(sep_wide, n = Inf)

# ------------------- 4. NON-PARAMETRIC TESTS ---------------------------------
kw_tbl <- do.call(rbind, lapply(FEATURES, function(f) {
  kt <- kruskal.test(spec_df[[f]] ~ spec_df$class)
  data.frame(feature = f, KW_chisq = round(unname(kt$statistic), 2),
             df = unname(kt$parameter), KW_p = signif(kt$p.value, 3))
})) %>% arrange(KW_p)
cat("\n== Kruskal-Wallis (four-class) per feature, ranked by p ==\n"); print(kw_tbl, n = Inf)

# Pairwise Wilcoxon rank-sum (Bonferroni-adjusted) for the top features
pairwise_wilcox <- function(df, feature) {
  suppressWarnings(pairwise.wilcox.test(df[[feature]], df$class, p.adjust.method = "bonferroni"))
}
cat("\n== Pairwise Wilcoxon rank-sum for the two best features ==\n")
for (f in head(kw_tbl$feature, 2)) { cat("\nFeature:", f, "\n"); print(pairwise_wilcox(spec_df, f)$p.value) }

# ------------------- 5. PCA --------------------------------------------------
X <- scale(as.matrix(spec_df[, FEATURES]))
pca <- prcomp(X, center = FALSE, scale. = FALSE)
cat("\n== PCA: variance explained ==\n"); print(summary(pca)$importance[, 1:min(6, ncol(X))])
cat("\n== PCA loadings (PC1-PC3): feature contributions to variance ==\n")
print(round(pca$rotation[, 1:3], 3))
# Feature importance from PCA = sum of squared loadings weighted by variance explained
ve <- (pca$sdev^2) / sum(pca$sdev^2)
pca_importance <- sort(rowSums(sweep(pca$rotation^2, 2, ve, "*")), decreasing = TRUE)
cat("\n== PCA-based feature importance (variance-weighted squared loadings) ==\n")
print(round(pca_importance, 3))

# ------------------- 6. glmnet VARIABLE IMPORTANCE ---------------------------
# Multinomial LASSO (alpha = 1) selects a sparse discriminating feature set;
# ridge (alpha = 0) retains correlated features. 10-fold CV, as in the seedlings analysis.
set.seed(123)
Xg <- as.matrix(spec_df[, FEATURES]); yg <- spec_df$class
lasso_cv <- cv.glmnet(Xg, yg, family = "multinomial", alpha = 1, nfolds = 10)
ridge_cv <- cv.glmnet(Xg, yg, family = "multinomial", alpha = 0, nfolds = 10)

# Importance = mean absolute LASSO coefficient across classes at lambda.1se
lasso_coef <- coef(lasso_cv, s = "lambda.1se")
imp <- sapply(FEATURES, function(f)
  mean(abs(sapply(lasso_coef, function(m) m[f, 1]))))
lasso_importance <- sort(imp[imp > 0], decreasing = TRUE)
cat("\n== glmnet LASSO variable importance (mean |coef| across classes, lambda.1se) ==\n")
print(round(lasso_importance, 4))
cat("\nFeatures with zero LASSO coefficient were dropped by the model as redundant.\n")

# ------------------- 7. CONSENSUS RANKING + RECOMMENDATION -------------------
rank_of <- function(x) setNames(rank(-x, ties.method = "min"), names(x))
consensus <- data.frame(feature = FEATURES) %>%
  left_join(sep_wide %>% transmute(feature, sep_rank = rank(-M_mean)), by = "feature") %>%
  left_join(kw_tbl   %>% transmute(feature, kw_rank  = rank(KW_p)),   by = "feature") %>%
  mutate(pca_rank   = rank_of(pca_importance)[feature],
         lasso_rank = ifelse(feature %in% names(lasso_importance),
                             rank_of(lasso_importance)[feature], NA)) %>%
  rowwise() %>%
  mutate(mean_rank = mean(c(sep_rank, kw_rank, pca_rank,
                            ifelse(is.na(lasso_rank), max(kw_rank), lasso_rank)), na.rm = TRUE)) %>%
  ungroup() %>% arrange(mean_rank)
cat("\n== CONSENSUS feature ranking (lower = better discriminator) ==\n")
print(consensus, n = Inf)
write.csv(consensus, "03.outputs/feature_selection_ranking.csv", row.names = FALSE)

cat("\nRecommendation: retain the features with a mean_rank in the top tier and a\n",
    "separability M > ~1 for the water-vs-vegetation and water-vs-soil pairs, which\n",
    "are the boundaries that define Lake Chilwa's littoral zone. Drop features the\n",
    "LASSO zeroed and that PCA shows loading on the same axis (redundant).\n",
    "Then rebuild the SMA endmembers and the RF feature stack on the retained set.\n")
