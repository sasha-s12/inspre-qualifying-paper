
library(inspre)
library(dplyr)
library(ggplot2)

set.seed(123)

B <- 100
D <- 10
v <- 0.3

# 1. Save original function and load corrected version
original <- getFromNamespace("predict_inspre", "inspre")

# 2. Store results
results <- matrix(
  NA_real_, nrow = B, ncol = 6,
  dimnames = list(
    NULL,
    c("Raw_MAE_IV", "Raw_MAE_Sparse",
      "MAE_IV", "MAE_Sparse",
      "SHD_IV", "SHD_Sparse")
  )
)

results_nonzero <- matrix(NA_real_, nrow = B, ncol = 4)

colnames(results_nonzero) <- c(
  "Magnitude_IV", "Magnitude_Sparse",
  "MAE_IV", "MAE_Sparse"
)

# 3. Run 100 simulations using corrected CV function
tryCatch({
  
  source("fixed_inspre_cv.R")
  
  for (b in 1:B) {
    
    # Generate data
    d <- generate_dataset(
      D = D,
      N_cont = 1000,
      N_int = 200,
      int_beta = -2,
      graph = "random",
      v = v,
      p = 1/9,
      DAG = TRUE,
      C = 0
    )
    
    # Fit model with five-fold CV
    fit <- fit_inspre_from_X(
      d$Y, d$targets,
      DAG = TRUE,
      cv_folds = 5,
      verbose = 0
    )
    
    # Select lambda
    best <- which.min(fit$eps_hat_G)
    
    # IV estimates
    R_IV <- fit$R_hat
    
    # Sparse reconstruction using U and V
    U_hat <- fit$U[, , best]
    V_hat <- fit$V[, , best]
    
    A <- diag(diag(V_hat)) %*% U_hat
    R_sparse <- sweep(A, 1, diag(A), "/")
    
    # True total causal effects
    R_true <- d$R
    
    diag(R_true) <- 0
    diag(R_IV) <- 0
    diag(R_sparse) <- 0
    
    # Raw MAE (without thresholding)
    offdiag <- row(R_true) != col(R_true)
    
    raw_MAE_IV <- mean(abs(R_IV[offdiag] - R_true[offdiag]))
    raw_MAE_sparse <- mean(abs(R_sparse[offdiag] - R_true[offdiag]))
    
    # Thresholded MAE and SHD
    metrics_IV <- calc_metrics(R_IV, R_true, eps = v/4)
    metrics_sparse <- calc_metrics(R_sparse, R_true, eps = v/4)
    
    # Store results
    results[b, ] <- c(
      raw_MAE_IV,
      raw_MAE_sparse,
      metrics_IV$mae,
      metrics_sparse$mae,
      metrics_IV$shd,
      metrics_sparse$shd
    )
    
    # Evaluate true nonzero effects
    nonzero <- (R_true != 0) &
      (row(R_true) != col(R_true))
    
    results_nonzero[b, ] <- c(
      mean(abs(R_IV[nonzero]) - abs(R_true[nonzero])),
      mean(abs(R_sparse[nonzero]) - abs(R_true[nonzero])),
      mean(abs(R_IV[nonzero] - R_true[nonzero])),
      mean(abs(R_sparse[nonzero] - R_true[nonzero]))
    )
    
    cat("Simulation", b, "of", B, "completed\n")
  }
  
}, finally = {
  # Restore original inspre function
  assignInNamespace("predict_inspre", original, ns = "inspre")
})

# 4. Summary table
summary_table <- data.frame(
  Estimator = c("IV", "Sparse"),
  
  Mean_Raw_MAE = c(
    mean(results[, "Raw_MAE_IV"]),
    mean(results[, "Raw_MAE_Sparse"])
  ),
  
  MCSE_Raw_MAE = c(
    sd(results[, "Raw_MAE_IV"]) / sqrt(B),
    sd(results[, "Raw_MAE_Sparse"]) / sqrt(B)
  ),
  
  Mean_MAE = c(
    mean(results[, "MAE_IV"]),
    mean(results[, "MAE_Sparse"])
  ),
  
  MCSE_MAE = c(
    sd(results[, "MAE_IV"]) / sqrt(B),
    sd(results[, "MAE_Sparse"]) / sqrt(B)
  ),
  
  Mean_SHD = c(
    mean(results[, "SHD_IV"]),
    mean(results[, "SHD_Sparse"])
  ),
  
  MCSE_SHD = c(
    sd(results[, "SHD_IV"]) / sqrt(B),
    sd(results[, "SHD_Sparse"]) / sqrt(B)
  )
)

print(summary_table)

# 5. Paired differences (positive values favour sparse)
gap_raw_MAE <- results[, "Raw_MAE_IV"] -
  results[, "Raw_MAE_Sparse"]

gap_MAE <- results[, "MAE_IV"] -
  results[, "MAE_Sparse"]

gap_SHD <- results[, "SHD_IV"] -
  results[, "SHD_Sparse"]

comparison <- data.frame(
  Metric = c("Raw MAE", "Thresholded MAE", "SHD"),
  
  Mean_Difference = c(
    mean(gap_raw_MAE),
    mean(gap_MAE),
    mean(gap_SHD)
  ),
  
  MCSE = c(
    sd(gap_raw_MAE) / sqrt(B),
    sd(gap_MAE) / sqrt(B),
    sd(gap_SHD) / sqrt(B)
  )
)

print(comparison)

# 6. Heatmap of final simulation replicate
genes <- as.character(1:D)

dimnames(R_true) <- list(genes, genes)
dimnames(R_IV) <- list(genes, genes)
dimnames(R_sparse) <- list(genes, genes)

# Prepare data
plot_data <- rbind(
  data.frame(as.table(R_true), Method = "True R"),
  data.frame(as.table(R_IV), Method = "Direct IV"),
  data.frame(as.table(R_sparse), Method = "Sparse reconstruction")
)

colnames(plot_data) <- c(
  "Gene_from", "Gene_to", "Effect", "Method"
)

plot_data$Gene_from <- factor(
  plot_data$Gene_from, levels = genes
)

plot_data$Gene_to <- factor(
  plot_data$Gene_to, levels = genes
)

plot_data$Method <- factor(
  plot_data$Method,
  levels = c(
    "True R", "Direct IV", "Sparse reconstruction"
  )
)

# Same colour scale for all three heatmaps
limit <- max(abs(plot_data$Effect))

p <- ggplot(
  plot_data,
  aes(x = Gene_to, y = Gene_from, fill = Effect)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~Method, nrow = 1) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-limit, limit)
  ) +
  coord_equal() +
  labs(
    x = "Gene",
    y = "Gene",
    fill = "Total effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

# 7. Nonzero-effect results table
nonzero_table <- data.frame(
  Metric = c("Mean magnitude error", "Raw MAE"),
  
  Direct_IV = c(
    mean(results_nonzero[, "Magnitude_IV"]),
    mean(results_nonzero[, "MAE_IV"])
  ),
  
  Sparse_Reconstruction = c(
    mean(results_nonzero[, "Magnitude_Sparse"]),
    mean(results_nonzero[, "MAE_Sparse"])
  ),
  
  MCSE_IV = c(
    sd(results_nonzero[, "Magnitude_IV"]),
    sd(results_nonzero[, "MAE_IV"])
  ) / sqrt(B),
  
  MCSE_Sparse = c(
    sd(results_nonzero[, "Magnitude_Sparse"]),
    sd(results_nonzero[, "MAE_Sparse"])
  ) / sqrt(B)
)

# Format as mean (MCSE)
formatted_nonzero_table <- data.frame(
  Metric = nonzero_table$Metric,
  
  Direct_IV = sprintf(
    "%.5f (%.5f)",
    nonzero_table$Direct_IV,
    nonzero_table$MCSE_IV
  ),
  
  Sparse_Reconstruction = sprintf(
    "%.5f (%.5f)",
    nonzero_table$Sparse_Reconstruction,
    nonzero_table$MCSE_Sparse
  )
)

print(formatted_nonzero_table, row.names = FALSE)

# 8. Record package versions for reproducibility
sessionInfo()
