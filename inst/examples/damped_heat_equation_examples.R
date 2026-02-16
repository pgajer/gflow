# ============================================================================
# Damped Heat Equation Examples in 1D and 2D
# ============================================================================
#
# This script demonstrates the behavior of the damped heat equation:
#   ∂ρ/∂t = -L₀ρ - β(ρ - u)
#
# through concrete examples that visualize density evolution.
# ============================================================================

library(Matrix)
library(ggplot2)
library(gridExtra)

# ============================================================================
# PART 1: HELPER FUNCTIONS
# ============================================================================

#' Build 1D chain graph Laplacian
#'
#' Creates the graph Laplacian for a 1D chain with n vertices.
#' Vertices are connected sequentially: 1--2--3--...--n
#' This gives a tridiagonal Laplacian matrix.
build_chain_laplacian <- function(n) {
  # Degree matrix: D[i,i] = degree of vertex i
  # Interior vertices have degree 2, endpoints have degree 1
  D <- diag(c(1, rep(2, n-2), 1))

  # Adjacency matrix: A[i,j] = 1 if i and j are neighbors
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    A[i, i+1] <- 1
    A[i+1, i] <- 1
  }

  # Laplacian: L = D - A
  L <- D - A
  return(L)
}

#' Build k-NN graph Laplacian from point cloud
#'
#' Given n points in R^d, connects each point to its k nearest neighbors
#' and builds the combinatorial graph Laplacian.
build_knn_laplacian <- function(X, k) {
  n <- nrow(X)

  # Compute pairwise distances
  dist_matrix <- as.matrix(dist(X, method = "euclidean"))

  # Build adjacency matrix via k-NN
  A <- matrix(0, n, n)
  for (i in 1:n) {
    # Find k nearest neighbors (excluding self)
    neighbors <- order(dist_matrix[i, ])[2:(k+1)]
    A[i, neighbors] <- 1
    A[neighbors, i] <- 1  # Symmetrize
  }

  # Make sure A is symmetric (it should be already, but enforce it)
  A <- (A + t(A)) / 2
  A[A > 0] <- 1

  # Compute Laplacian
  D <- diag(rowSums(A))
  L <- D - A

  return(L)
}

#' Solve one step of damped heat equation using implicit Euler
#'
#' Solves: (I + t(L + βI))ρ_new = ρ_old + tβu
#' where u is the uniform distribution vector (all ones).
damped_heat_step <- function(rho_current, L, t, beta) {
  n <- length(rho_current)

  # System matrix: A = I + t(L + βI)
  I <- diag(n)
  A <- I + t * (L + beta * I)

  # Right-hand side: b = ρ_current + tβu
  u <- rep(1, n)
  b <- rho_current + t * beta * u

  # Solve linear system
  rho_new <- solve(A, b)

  # Normalize to sum to n
  rho_new <- rho_new * n / sum(rho_new)

  # Enforce positivity (clip small negatives)
  rho_new[rho_new < 0] <- 1e-15

  return(rho_new)
}

#' Evolve density through multiple time steps
evolve_density <- function(rho_init, L, t_step, beta, n_steps) {
  n <- length(rho_init)
  history <- matrix(0, nrow = n_steps + 1, ncol = n)
  history[1, ] <- rho_init

  rho <- rho_init
  for (step in 1:n_steps) {
    rho <- damped_heat_step(rho, L, t_step, beta)
    history[step + 1, ] <- rho
  }

  return(history)
}

# ============================================================================
# PART 2: 1D EXAMPLES
# ============================================================================

cat("=" , rep("=", 70), "\n", sep="")
cat("PART 2: 1D Chain Graph Examples\n")
cat("=" , rep("=", 70), "\n", sep="")

# Setup: 1D chain with 50 vertices
n_1d <- 50
L_1d <- build_chain_laplacian(n_1d)

# Compute spectral gap for automatic parameter selection
eig <- eigen(L_1d, symmetric = TRUE)
eigenvalues <- sort(eig$values)
lambda_2 <- eigenvalues[2]  # Second smallest eigenvalue (first is ~0)

cat("\nSpectral gap λ₂ =", lambda_2, "\n")
cat("Suggested t range: [", 0.1/lambda_2, ",", 1.0/lambda_2, "]\n")

# Example 1: Localized initial density (spike at center)
cat("\n--- Example 1A: Localized density, NO damping (β = 0) ---\n")
rho_init_spike <- rep(0.1, n_1d)
rho_init_spike[25] <- 10  # Spike at center
rho_init_spike <- rho_init_spike * n_1d / sum(rho_init_spike)

t_step <- 0.5 / lambda_2
beta_none <- 0.0
n_steps <- 30

history_undamped <- evolve_density(rho_init_spike, L_1d, t_step, beta_none, n_steps)

cat("Initial entropy:", -sum((rho_init_spike/n_1d) * log(rho_init_spike/n_1d)), "\n")
cat("Final entropy:", -sum((history_undamped[n_steps+1, ]/n_1d) *
                             log(history_undamped[n_steps+1, ]/n_1d)), "\n")
cat("Final density range: [", min(history_undamped[n_steps+1, ]), ",",
    max(history_undamped[n_steps+1, ]), "]\n")

# Example 1B: Same initial density WITH damping
cat("\n--- Example 1B: Localized density, WITH damping (β = 0.1/t) ---\n")
beta_damped <- 0.1 / t_step

history_damped <- evolve_density(rho_init_spike, L_1d, t_step, beta_damped, n_steps)

cat("Initial entropy:", -sum((rho_init_spike/n_1d) * log(rho_init_spike/n_1d)), "\n")
cat("Final entropy:", -sum((history_damped[n_steps+1, ]/n_1d) *
                             log(history_damped[n_steps+1, ]/n_1d)), "\n")
cat("Final density range: [", min(history_damped[n_steps+1, ]), ",",
    max(history_damped[n_steps+1, ]), "]\n")

# Create visualization comparing damped vs undamped
plot_data_1d <- data.frame(
  vertex = rep(1:n_1d, 4),
  density = c(rho_init_spike,
              history_undamped[n_steps+1, ],
              rho_init_spike,
              history_damped[n_steps+1, ]),
  type = rep(c("Undamped", "Undamped", "Damped", "Damped"), each = n_1d),
  time = rep(c("Initial", "Final", "Initial", "Final"), each = n_1d)
)

plot_data_1d$time <- factor(plot_data_1d$time, levels = c("Initial", "Final"))

p1 <- ggplot(plot_data_1d, aes(x = vertex, y = density, color = time)) +
  geom_line(linewidth = 1) +
  facet_wrap(~type, ncol = 1) +
  labs(title = "1D Chain: Damped vs Undamped Heat Equation",
       subtitle = sprintf("Localized initial spike (t = %.3f, β = %.3f for damped)",
                         t_step, beta_damped),
       x = "Vertex Position", y = "Density ρ(x)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(p1)

# Example 2: Two-cluster initial density
cat("\n--- Example 2: Two-cluster initial density ---\n")
rho_init_clusters <- rep(0.5, n_1d)
rho_init_clusters[10:15] <- 3
rho_init_clusters[35:40] <- 3
rho_init_clusters <- rho_init_clusters * n_1d / sum(rho_init_clusters)

history_clusters_undamped <- evolve_density(rho_init_clusters, L_1d, t_step, 0.0, n_steps)
history_clusters_damped <- evolve_density(rho_init_clusters, L_1d, t_step, beta_damped, n_steps)

plot_data_clusters <- data.frame(
  vertex = rep(1:n_1d, 4),
  density = c(rho_init_clusters,
              history_clusters_undamped[n_steps+1, ],
              rho_init_clusters,
              history_clusters_damped[n_steps+1, ]),
  type = rep(c("Undamped", "Undamped", "Damped", "Damped"), each = n_1d),
  time = rep(c("Initial", "Final", "Initial", "Final"), each = n_1d)
)

plot_data_clusters$time <- factor(plot_data_clusters$time, levels = c("Initial", "Final"))

p2 <- ggplot(plot_data_clusters, aes(x = vertex, y = density, color = time)) +
  geom_line(linewidth = 1) +
  facet_wrap(~type, ncol = 1) +
  labs(title = "1D Chain: Two-Cluster Initial Density",
       subtitle = "Undamped preserves bimodality; Damped smooths toward uniform",
       x = "Vertex Position", y = "Density ρ(x)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(p2)

# ============================================================================
# PART 3: 2D EXAMPLES
# ============================================================================

cat("\n", rep("=", 72), "\n", sep="")
cat("PART 3: 2D Point Cloud Examples\n")
cat(rep("=", 72), "\n", sep="")

# Example 3A: Uniform random points in [-1, 1]^2
cat("\n--- Example 3A: Uniform random sample in [-1,1]^2 ---\n")
set.seed(42)
n_2d <- 100
X_uniform <- matrix(runif(n_2d * 2, -1, 1), ncol = 2)

k_nn <- 8
L_uniform <- build_knn_laplacian(X_uniform, k_nn)

# Spectral analysis
eig_uniform <- eigen(L_uniform, symmetric = TRUE)
lambda_2_uniform <- sort(eig_uniform$values)[2]
cat("Spectral gap λ₂ =", lambda_2_uniform, "\n")

# Localized initial density (spike at center)
center_idx <- which.min(rowSums((X_uniform - c(0, 0))^2))
rho_init_2d <- rep(0.5, n_2d)
rho_init_2d[center_idx] <- 20
rho_init_2d <- rho_init_2d * n_2d / sum(rho_init_2d)

t_step_2d <- 0.5 / lambda_2_uniform
beta_2d <- 0.1 / t_step_2d
n_steps_2d <- 20

history_2d_undamped <- evolve_density(rho_init_2d, L_uniform, t_step_2d, 0.0, n_steps_2d)
history_2d_damped <- evolve_density(rho_init_2d, L_uniform, t_step_2d, beta_2d, n_steps_2d)

# Visualization
plot_data_2d_uniform <- rbind(
  data.frame(X_uniform,
             density = rho_init_2d,
             type = "Initial",
             damping = "Both"),
  data.frame(X_uniform,
             density = history_2d_undamped[n_steps_2d+1, ],
             type = "Final",
             damping = "Undamped (β=0)"),
  data.frame(X_uniform,
             density = history_2d_damped[n_steps_2d+1, ],
             type = "Final",
             damping = sprintf("Damped (β=%.3f)", beta_2d))
)

colnames(plot_data_2d_uniform)[1:2] <- c("x", "y")

p3 <- ggplot(plot_data_2d_uniform, aes(x = x, y = y, color = density, size = density)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~type + damping, ncol = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(1, 8)) +
  labs(title = "2D Uniform Random Points: k-NN Graph (k=8)",
       subtitle = "Central spike diffuses and spreads") +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_fixed()

print(p3)

# Example 3B: Structured/sparse point cloud (two well-separated clusters)
cat("\n--- Example 3B: Two-cluster structure in 2D ---\n")
set.seed(123)
n_cluster <- 50

# Cluster 1: centered at (-0.5, -0.5)
cluster1 <- matrix(rnorm(n_cluster * 2, mean = 0, sd = 0.15), ncol = 2) +
  matrix(c(-0.5, -0.5), nrow = n_cluster, ncol = 2, byrow = TRUE)

# Cluster 2: centered at (0.5, 0.5)
cluster2 <- matrix(rnorm(n_cluster * 2, mean = 0, sd = 0.15), ncol = 2) +
  matrix(c(0.5, 0.5), nrow = n_cluster, ncol = 2, byrow = TRUE)

X_clusters <- rbind(cluster1, cluster2)
n_2d_clusters <- nrow(X_clusters)

L_clusters <- build_knn_laplacian(X_clusters, k_nn)

eig_clusters <- eigen(L_clusters, symmetric = TRUE)
lambda_2_clusters <- sort(eig_clusters$values)[2]
cat("Spectral gap λ₂ =", lambda_2_clusters, "\n")
cat("Note: Smaller spectral gap indicates bottleneck between clusters\n")

# Initial density: concentrated in cluster 1
rho_init_clusters_2d <- rep(0.2, n_2d_clusters)
rho_init_clusters_2d[1:n_cluster] <- 3  # High density in cluster 1
rho_init_clusters_2d <- rho_init_clusters_2d * n_2d_clusters / sum(rho_init_clusters_2d)

t_step_clusters <- 1.0 / lambda_2_clusters  # Longer time due to bottleneck
beta_clusters <- 0.1 / t_step_clusters

history_clusters_2d_undamped <- evolve_density(rho_init_clusters_2d, L_clusters,
                                               t_step_clusters, 0.0, n_steps_2d)
history_clusters_2d_damped <- evolve_density(rho_init_clusters_2d, L_clusters,
                                             t_step_clusters, beta_clusters, n_steps_2d)

plot_data_2d_clusters <- rbind(
  data.frame(X_clusters,
             density = rho_init_clusters_2d,
             type = "Initial",
             damping = "Both"),
  data.frame(X_clusters,
             density = history_clusters_2d_undamped[n_steps_2d+1, ],
             type = "Final",
             damping = "Undamped (β=0)"),
  data.frame(X_clusters,
             density = history_clusters_2d_damped[n_steps_2d+1, ],
             type = "Final",
             damping = sprintf("Damped (β=%.3f)", beta_clusters))
)

colnames(plot_data_2d_clusters)[1:2] <- c("x", "y")

p4 <- ggplot(plot_data_2d_clusters, aes(x = x, y = y, color = density, size = density)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~type + damping, ncol = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(1, 8)) +
  labs(title = "2D Two-Cluster Structure: k-NN Graph (k=8)",
       subtitle = "Density flows between clusters; damping prevents over-concentration") +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_fixed()

print(p4)

# ============================================================================
# PART 4: QUANTITATIVE COMPARISON
# ============================================================================

cat("\n", rep("=", 72), "\n", sep="")
cat("PART 4: Quantitative Comparison of Damping Effects\n")
cat(rep("=", 72), "\n", sep="")

compute_metrics <- function(rho, n) {
  # Entropy
  probs <- rho / n
  probs[probs < 1e-15] <- 1e-15  # Avoid log(0)
  entropy <- -sum(probs * log(probs))

  # Concentration: ratio of max to mean
  concentration <- max(rho) / mean(rho)

  # Uniformity: L2 distance to uniform
  uniformity <- sqrt(mean((rho - 1)^2))

  return(list(entropy = entropy, concentration = concentration, uniformity = uniformity))
}

cat("\n1D Chain Examples:\n")
cat("  Spike initial (undamped final):\n")
metrics <- compute_metrics(history_undamped[n_steps+1, ], n_1d)
cat(sprintf("    Entropy: %.4f, Concentration: %.2fx, Uniformity: %.4f\n",
           metrics$entropy, metrics$concentration, metrics$uniformity))

cat("  Spike initial (damped final):\n")
metrics <- compute_metrics(history_damped[n_steps+1, ], n_1d)
cat(sprintf("    Entropy: %.4f, Concentration: %.2fx, Uniformity: %.4f\n",
           metrics$entropy, metrics$concentration, metrics$uniformity))

cat("\n2D Uniform Examples:\n")
cat("  Central spike (undamped final):\n")
metrics <- compute_metrics(history_2d_undamped[n_steps_2d+1, ], n_2d)
cat(sprintf("    Entropy: %.4f, Concentration: %.2fx, Uniformity: %.4f\n",
           metrics$entropy, metrics$concentration, metrics$uniformity))

cat("  Central spike (damped final):\n")
metrics <- compute_metrics(history_2d_damped[n_steps_2d+1, ], n_2d)
cat(sprintf("    Entropy: %.4f, Concentration: %.2fx, Uniformity: %.4f\n",
           metrics$entropy, metrics$concentration, metrics$uniformity))

cat("\n2D Cluster Examples:\n")
cat("  Cluster imbalance (undamped final):\n")
metrics <- compute_metrics(history_clusters_2d_undamped[n_steps_2d+1, ], n_2d_clusters)
cat(sprintf("    Entropy: %.4f, Concentration: %.2fx, Uniformity: %.4f\n",
           metrics$entropy, metrics$concentration, metrics$uniformity))

cat("  Cluster imbalance (damped final):\n")
metrics <- compute_metrics(history_clusters_2d_damped[n_steps_2d+1, ], n_2d_clusters)
cat(sprintf("    Entropy: %.4f, Concentration: %.2fx, Uniformity: %.4f\n",
           metrics$entropy, metrics$concentration, metrics$uniformity))

cat("\n", rep("=", 72), "\n", sep="")
cat("KEY INSIGHTS:\n")
cat(rep("=", 72), "\n", sep="")
cat("1. Undamped diffusion (β=0) smooths density but converges very slowly\n")
cat("   to uniform distribution. Can maintain bimodality/structure.\n\n")
cat("2. Damped diffusion (β>0) adds explicit pull toward uniformity,\n")
cat("   achieving faster convergence while still respecting geometry.\n\n")
cat("3. In clustered data, damping prevents over-concentration within\n")
cat("   clusters while allowing some geometric differentiation.\n\n")
cat("4. The spectral gap λ₂ indicates graph connectivity: smaller λ₂\n")
cat("   (bottlenecks) requires longer diffusion time t.\n\n")
cat("5. The choice β ≈ 0.1/t balances geometric sensitivity with stability.\n")
cat(rep("=", 72), "\n", sep="")
