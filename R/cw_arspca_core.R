# ============================================================
# CW-ARSPCA CORE FUNCTIONS
# ============================================================

# --- Step 1: Robust location & scale ---
robust_location_scale <- function(X) {
  mu <- apply(X, 2, median, na.rm = TRUE)
  s <- apply(X, 2, function(x) {
    qn <- tryCatch(robustbase::Qn(x, finite.corr = TRUE), error = function(e) NA)
    if (is.na(qn) || qn < 1e-10) {
      iqr_val <- IQR(x, na.rm = TRUE) / 1.349
      if (iqr_val < 1e-10) iqr_val <- sd(x, na.rm = TRUE)
      if (iqr_val < 1e-10) iqr_val <- 1e-10
      return(iqr_val)
    }
    return(qn)
  })
  list(mu = mu, s = s)
}

# --- Step 2: Weight function (Three-zone) ---
weight_three_zone <- function(z, c1 = 2.5, c2 = 4.0) {
  w <- numeric(length(z)); az <- abs(z)
  w[az <= c1] <- 1
  idx <- az > c1 & az <= c2
  w[idx] <- ((c2 - az[idx]) / (c2 - c1))^2
  w[az > c2] <- 0
  w
}

# --- Step 3: Cellwise and row weights ---
cellwise_weights <- function(X, mu, s, weight_fn = weight_three_zone, ...) {
  Z <- sweep(sweep(X, 2, mu, "-"), 2, s, "/")
  W_cell <- matrix(weight_fn(as.vector(Z), ...), nrow = nrow(X), ncol = ncol(X))
  W_row <- if (ncol(X) >= 10) apply(W_cell, 1, median) else rowMeans(W_cell)
  list(W_cell = W_cell, W_row = W_row, Z = Z)
}

# --- Step 4: Weighted covariance ---
weighted_cov <- function(X, mu, W_row) {
  Xc <- sweep(X, 2, mu, "-")
  crossprod(Xc * sqrt(W_row)) / sum(W_row)
}

# --- Step 5: Sparse PCA iteration (IRLS) ---
soft_threshold <- function(x, lambda) sign(x) * pmax(abs(x) - lambda, 0)
huber_phi <- function(t, c3 = 1.345) ifelse(abs(t) <= c3, t^2/2, c3*abs(t) - c3^2/2)

sparse_pca_iter <- function(Sigma, X, mu, v, W_prev = NULL,
                            lambda1 = 0.05, lambda2 = 0.01, c3 = 1.345,
                            max_iter = 500, tol = 1e-6) {
  p <- ncol(X); w <- eigen(Sigma)$vectors[, 1]
  for (iter in 1:max_iter) {
    w_old <- w
    y_vec <- as.vector((X - matrix(mu, nrow(X), p, byrow = TRUE)) %*% w)
    sigma_w <- 1.4826 * median(abs(y_vec - median(y_vec)))
    if (sigma_w < 1e-10) sigma_w <- 1e-10
    t <- y_vec / sigma_w
    phi_prime <- ifelse(abs(t) <= c3, 1, c3 / pmax(abs(t), 1e-10))
    weights <- v * phi_prime
    Xc <- sweep(X, 2, mu, "-")
    grad <- crossprod(Xc * weights, Xc %*% w) / sum(weights)
    w <- soft_threshold(grad, lambda1) / (1 + lambda2)
    if (!is.null(W_prev) && ncol(W_prev) > 0)
      w <- w - W_prev %*% (t(W_prev) %*% w)
    w_norm <- sqrt(sum(w^2))
    if (w_norm < 1e-10) w <- eigen(Sigma)$vectors[, 1] else w <- w / w_norm
    if (sum((w - w_old)^2) < tol) break
  }
  list(w = w, iter = iter)
}

# --- Main Wrapper Function ---
CW_ARSPCA <- function(X, n_components = 5, c1 = 2.5, c2 = 4.0, c3 = 1.345,
                      lambda1 = 0.1, lambda2 = 0.01,
                      n_starts = 5, verbose = FALSE) {
  X_scaled <- scale(X)
  X_center <- attr(X_scaled, "scaled:center")
  X_scale  <- attr(X_scaled, "scaled:scale")
  ls <- robust_location_scale(X_scaled)
  cw <- cellwise_weights(X_scaled, ls$mu, ls$s, weight_three_zone, c1=c1, c2=c2)

  loadings <- matrix(0, nrow = ncol(X), ncol = n_components)
  X_curr <- X_scaled; eigenvalues <- numeric(n_components)

  for (k in 1:n_components) {
    Sigma_curr <- weighted_cov(X_curr, ls$mu, cw$W_row)
    best <- NULL; best_eigen <- -Inf
    W_prev <- if (k > 1) loadings[, 1:(k-1), drop = FALSE] else NULL
    for (s in 1:n_starts) {
      res <- sparse_pca_iter(Sigma_curr, X_curr, ls$mu, cw$W_row,
                             W_prev=W_prev, lambda1=lambda1,
                             lambda2=lambda2, c3=c3)
      es <- var(X_curr %*% res$w)
      if (es > best_eigen) { best_eigen <- es; best <- res }
    }
    w_k <- best$w
    if (!is.null(W_prev) && ncol(W_prev) > 0)
      w_k <- w_k - W_prev %*% (t(W_prev) %*% w_k)
    w_k <- w_k / sqrt(sum(w_k^2))
    loadings[, k] <- w_k
    proj <- X_curr %*% w_k
    eigenvalues[k] <- var(proj)
    X_curr <- X_curr - proj %*% t(w_k)
    if (verbose) cat("Component", k, "converged in", best$iter,
                     "iterations, eigenvalue =", round(eigenvalues[k], 4), "\n")
  }

  if (any(diff(eigenvalues) > 0)) {
    ord <- order(eigenvalues, decreasing = TRUE)
    loadings <- loadings[, ord, drop = FALSE]
    eigenvalues <- eigenvalues[ord]
  }

  result <- list(loadings = loadings, scores = X_scaled %*% loadings,
                 outliers = cw$W_cell < 0.5, eigenvalues = eigenvalues,
                 X_center = X_center, X_scale = X_scale)
  class(result) <- "cwarspca"
  return(result)
}
