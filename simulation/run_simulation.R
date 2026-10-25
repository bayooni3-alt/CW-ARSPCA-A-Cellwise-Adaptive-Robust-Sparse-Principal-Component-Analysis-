# ============================================================
# SIMULATION STUDY: CW-ARSPCA
# ============================================================
# Note: Run source("R/cw_arspca_core.R") before running this script.

# 1. DATA GENERATOR
generate_data <- function(n, p, k = 3, epsilon = 0.1,
                          contamination = "cellwise",
                          outlier_magnitude = 5) {
  # True sparse loadings
  W <- matrix(0, p, k)
  for (l in 1:k) {
    idx <- sample(p, 5)
    W[idx, l] <- rnorm(5)
    W[, l] <- W[, l] / sqrt(sum(W[, l]^2))
  }
  
  # Latent scores and clean data
  A <- matrix(rnorm(n * k), n, k)
  X <- A %*% t(W) + matrix(rnorm(n * p, sd = 0.5), n, p)
  
  # Add contamination
  if (contamination %in% c("cellwise", "both")) {
    n_cell <- floor(n * p * epsilon / ifelse(contamination == "both", 2, 1))
    idx <- sample(n * p, n_cell)
    X[idx] <- X[idx] + rnorm(n_cell, sd = outlier_magnitude)
  }
  if (contamination %in% c("rowwise", "both")) {
    n_row <- floor(n * epsilon / ifelse(contamination == "both", 2, 1))
    idx <- sample(n, n_row)
    X[idx, ] <- X[idx, ] + matrix(rnorm(n_row * p, sd = outlier_magnitude), n_row, p)
  }
  
  list(X = X, W = W, A = A)
}

# 2. METRICS AND RUNNER
subspace_error <- function(W_est, W_true) {
  k <- ncol(W_true)
  W_est <- W_est[, 1:k, drop = FALSE]
  P_est <- W_est %*% solve(t(W_est) %*% W_est) %*% t(W_est)
  P_true <- W_true %*% t(W_true)
  norm(P_est - P_true, "F")
}

run_one <- function(n, p, epsilon, contamination, rep_id) {
  set.seed(rep_id)
  d <- generate_data(n, p, k = 3, epsilon, contamination)
  X <- d$X; W_true <- d$W
  
  out <- list()
  
  # Classical PCA
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  out$PCA <- subspace_error(pca$rotation, W_true)
  
  # ROBPCA
  rob <- tryCatch(PcaHubert(X, k = 3, alpha = 0.75), error = function(e) NULL)
  out$ROBPCA <- if (is.null(rob)) NA else subspace_error(rob@loadings, W_true)
  
  # CW-ARSPCA
  cw <- CW_ARSPCA(X, n_components = 3)
  out$CW_ARSPCA <- subspace_error(cw$loadings, W_true)
  
  data.frame(
    n = n, p = p, epsilon = epsilon,
    contamination = contamination, rep = rep_id,
    method = names(out), error = unlist(out)
  )
}

# 3. RUN ALL 108 SCENARIOS (Uncomment to execute)
# scenarios <- expand.grid(
#   n = c(100, 500, 1000), p = c(10, 50, 100),
#   epsilon = c(0, 0.05, 0.1, 0.2), contamination = c("cellwise", "rowwise", "both")
# )
# all_results <- do.call(rbind, lapply(1:nrow(scenarios), function(i) {
#   s <- scenarios[i, ]
#   do.call(rbind, lapply(1:100, function(r) run_one(s$n, s$p, s$epsilon, s$contamination, r)))
# }))
# write.csv(all_results, "simulation_results.csv")

# 4. ABLATION STUDY
ablation_study <- function(n = 500, p = 50, epsilon = 0.1, rep = 50) {
  # Implementation of variants (full, no_cellwise, no_huber, no_net)
}

# 5. WEIGHT FUNCTION COMPARISON
weight_comparison <- function(n = 500, p = 50, epsilon = 0.1, rep = 50) {
  # Implementation comparing three_zone, tukey, huber, hampel
}
