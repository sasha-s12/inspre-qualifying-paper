
# Corrected cross-validation prediction function for inspre
# Uses G_hat for structural-model prediction rather than R_hat.

predict_fixed <- function(res, .X, .beta, .targets) {
  vars <- intersect(colnames(.X), names(.beta))
  .beta <- .beta[vars]
  .X <- .X[, vars, drop = FALSE]
  
  ZB <- outer(.targets, colnames(.X), "==")
  ZB <- t(t(ZB) * .beta)
  
  eps_G <- eps_I <- numeric(length(res$lambda))
  
  for (i in seq_along(res$lambda)) {
    G <- res$G_hat[, , i]
    A <- solve(diag(nrow(G)) - G)
    
    eps_G[i] <- mean((.X - (.X %*% G + ZB))^2)
    eps_I[i] <- mean((.X - ZB %*% A)^2)
  }
  
  list(eps_hat_G = eps_G, eps_hat_I = eps_I)
}

# Temporarily replace the package's prediction function
assignInNamespace("predict_inspre", predict_fixed, ns = "inspre")
