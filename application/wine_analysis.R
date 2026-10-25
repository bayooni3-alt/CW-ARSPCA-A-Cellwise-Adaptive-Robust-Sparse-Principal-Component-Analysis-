# ============================================================
# APPLICATION: WINE DATA ANALYSIS
# ============================================================
# Note: Run source("R/cw_arspca_core.R") before running this script.

# 0. LOAD LIBRARIES
if (!require("rattle"))      install.packages("rattle")
if (!require("tidyverse"))   install.packages("tidyverse")
if (!require("rrcov"))       install.packages("rrcov")
if (!require("robustbase"))  install.packages("robustbase")
if (!require("xgboost"))     install.packages("xgboost")
if (!require("pROC"))        install.packages("pROC")
if (!require("rospca"))      install.packages("rospca")
if (!require("cellWise"))    install.packages("cellWise")
if (!require("nsprcomp"))    install.packages("nsprcomp")

library(rattle)
library(dplyr); library(tidyr); library(tibble)
library(tidyverse); library(rrcov); library(robustbase)
library(xgboost); library(pROC)
library(rospca); library(cellWise); library(nsprcomp)

select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; summarise <- dplyr::summarise
arrange <- dplyr::arrange; distinct <- dplyr::distinct

set.seed(42)

# 1. LOAD WINE DATA
data(wine)
X <- as.matrix(wine[, 2:14])
y <- as.numeric(wine$Type == "3")
X_scaled <- scale(X)
X_scaled[!is.finite(X_scaled)] <- 0
storage.mode(X_scaled) <- "double"
ncomp <- min(ncol(X_scaled) - 1, 5)

# 2. HELPER FUNCTION: MANUAL ROSPCA
rospca_manual <- function(X, k = 5, alpha = 0.75) {
  mcd <- rrcov::CovMcd(X, alpha = alpha)
  eig <- eigen(mcd@cov, symmetric = TRUE)
  k_use <- min(k, ncol(X))
  loadings <- eig$vectors[, 1:k_use, drop = FALSE]
  list(loadings = loadings, scores = sweep(X, 2, mcd@center, "-") %*% loadings,
       eigenvalues = eig$values[1:k_use], center = mcd@center, mu = mcd@center, 
       sigma = sqrt(diag(mcd@cov)))
}

# 3. RUN CW-ARSPCA
cat("\n========================================\n")
cat("RUN CW-ARSPCA\n")
cat("========================================\n")
cw_result <- CW_ARSPCA(X_scaled, n_components = ncomp, lambda1 = 0.1, lambda2 = 0.01, n_starts = 5, verbose = TRUE)
cw_var <- cw_result$eigenvalues / sum(cw_result$eigenvalues) * 100

# 4. VARIANCE COMPARISON (ALL 6 METHODS)
fix_len <- function(v, n) {
  if (length(v) > n) v <- v[1:n]
  if (length(v) < n) v <- c(v, rep(NA, n - length(v)))
  v
}
pca_result <- prcomp(X_scaled, center = TRUE, scale. = FALSE)
pca_var <- fix_len(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, ncomp)
rob_result <- PcaHubert(X_scaled, k = ncomp, alpha = 0.75)
rob_var <- fix_len(rob_result@eigenvalues / sum(rob_result@eigenvalues) * 100, ncomp)
spca_result <- nsprcomp(X_scaled, ncomp = ncomp, kappa = 0.5)
spca_var <- fix_len(apply(X_scaled %*% spca_result$rotation, 2, var) / sum(apply(X_scaled, 2, var)) * 100, ncomp)
rospca_result <- rospca_manual(X_scaled, k = ncomp, alpha = 0.75)
rospca_var <- fix_len(rospca_result$eigenvalues / sum(rospca_result$eigenvalues) * 100, ncomp)
macropca_result <- MacroPCA(X_scaled, k = ncomp)
macropca_var <- fix_len(macropca_result$eigenvalues / sum(macropca_result$eigenvalues) * 100, ncomp)
cw_var <- fix_len(cw_var, ncomp)

var_comparison <- data.frame(
  Component = paste0("PC", 1:ncomp),
  CW_ARSPCA = round(cw_var, 2), Classical_PCA = round(pca_var, 2),
  ROBPCA = round(rob_var, 2), SPCA = round(spca_var, 2),
  ROSPCA = round(rospca_var, 2), MacroPCA = round(macropca_var, 2)
)
print(var_comparison)

# 5. DOWNSTREAM ML (5-FOLD CV)
run_ml <- function(X, y, methods, n_folds = 5) {
  set.seed(42)
  folds <- sample(rep(1:n_folds, length.out = nrow(X)))
  per_fold <- lapply(methods, function(m) rep(NA_real_, n_folds))
  names(per_fold) <- methods

  for (k in 1:n_folds) {
    train <- which(folds != k); test <- which(folds == k)
    X_train <- X[train, , drop = FALSE]; X_test <- X[test, , drop = FALSE]
    y_train <- y[train]; y_test <- y[test]

    fit_xgb <- function(tr_x, te_x) {
      dtrain <- xgb.DMatrix(data = tr_x, label = y_train)
      dtest  <- xgb.DMatrix(data = te_x,  label = y_test)
      m <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 6, eta = 0.1), data = dtrain, nrounds = 100, verbose = 0)
      predict(m, dtest)
    }
    safe_auc <- function(pred) { if (is.null(pred)) return(NA_real_); as.numeric(roc(y_test, pred, quiet = TRUE)$auc) }

    per_fold$Raw[k] <- safe_auc(fit_xgb(X_train, X_test))
    
    p <- prcomp(X_train, center = TRUE, scale. = FALSE)
    per_fold$Classical_PCA[k] <- safe_auc(fit_xgb(p$x[, 1:ncomp, drop = FALSE], predict(p, X_test)[, 1:ncomp, drop = FALSE]))
    
    r <- PcaHubert(X_train, k = ncomp, alpha = 0.75)
    te_rob <- scale(X_test, center = r@center, scale = r@scale) %*% r@loadings[, 1:ncomp, drop = FALSE]
    per_fold$ROBPCA[k] <- safe_auc(fit_xgb(r@scores[, 1:ncomp, drop = FALSE], te_rob))
    
    s <- nsprcomp(X_train, ncomp = ncomp, kappa = 0.5)
    per_fold$SPCA[k] <- safe_auc(fit_xgb(X_train %*% s$rotation, X_test %*% s$rotation))
    
    r_m <- rospca_manual(X_train, k = ncomp, alpha = 0.75)
    te_m <- sweep(X_test, 2, r_m$center, "-") %*% r_m$loadings
    per_fold$ROSPCA[k] <- safe_auc(fit_xgb(r_m$scores[, 1:ncomp], te_m[, 1:ncomp]))
    
    mad_vals <- apply(X_train, 2, mad)
    keep <- mad_vals > max(0.05, 0.1 * median(mad_vals))
    if (sum(keep) >= 5) {
      m_macro <- MacroPCA(X_train[, keep], k = min(ncomp, sum(keep) - 1))
      te_macro <- scale(X_test[, keep], center = m_macro$mu, scale = m_macro$scale) %*% m_macro$loadings
      per_fold$MacroPCA[k] <- safe_auc(fit_xgb(m_macro$scores[, 1:min(ncomp, sum(keep)-1)], te_macro))
    }
    
    cw <- CW_ARSPCA(X_train, n_components = ncomp, lambda1 = 0.1, lambda2 = 0.01, n_starts = 5)
    te_cw <- scale(X_test, center = cw$X_center, scale = cw$X_scale) %*% cw$loadings
    per_fold$CW_ARSPCA[k] <- safe_auc(fit_xgb(cw$scores, te_cw))
  }
  return(per_fold)
}

methods <- c("Raw", "Classical_PCA", "ROBPCA", "SPCA", "ROSPCA", "MacroPCA", "CW_ARSPCA")
per_fold <- run_ml(X_scaled, y, methods, n_folds = 5)

# 6. SUMMARY & TABLES
safe_mean <- function(x) mean(x, na.rm = TRUE)
safe_sd   <- function(x) sd(x, na.rm = TRUE)

summary_df <- data.frame(
  Method   = names(per_fold),
  Mean_AUC = round(sapply(per_fold, safe_mean), 4),
  SD_AUC   = round(sapply(per_fold, safe_sd), 4)
)
summary_df <- summary_df[order(-summary_df$Mean_AUC), ]
print(summary_df, row.names = FALSE)

# 7. OUTLIER ANALYSIS
outlier_count <- rowSums(cw_result$outliers)
cat("\nCW-ARSPCA identified", sum(cw_result$outliers), "cellwise outliers across", sum(outlier_count > 0), "rows.\n")
