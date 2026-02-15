# Test Stage 0: Initialization only
# Tests: initialize_from_knn()
#   - kNN graph construction
#   - Geometric pruning
#   - Initial density computation
#   - Initial metric construction
#   - Initial Laplacian assembly


source("/Users/pgajer/current_projects/gflow/tests/manual/00-setup.R")

cat("=" , rep("=", 70), "\n", sep = "")
cat("TEST STAGE 0: INITIALIZATION\n")
cat("=" , rep("=", 70), "\n", sep = "")

# Generate test data
test.data <- generate.test.data.1d(n.pts = 50, noise_sd = 1)
cat("\nTest data: n =", nrow(test.data$X), "\n")

# Plot test data
plot.test.data(test.data)

timestamp <- generate.timestamp()
file <- "~/current_projects/gflow/tests/manual/pics/"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=12, height=6)
op <- par(mar=c(3.5,3.5,1.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot.test.data(test.data)
par(op)
dev.off()
system(paste0("open ",file))
## "~/current_projects/gflow/tests/manual/pics/2025_10_07_175445.pdf"

# Run initialization only (test_stage = 0)
cat("\nRunning fit_riem_dcx_regression with test_stage = 0...\n")

rdcx.res <- fit.knn.riem.graph.regression(
    test.data$X,
    test.data$y,
    k = 3,
    n.eigenpairs = 10,
    filter.type = "heat_kernel", # "tikhonov", "cubic_spline"),
    test.stage = -1
)

cat("\n--- Structure of result ---\n")
str(rdcx.res, max.level = 1)

plot(test.data$X[,1], rdcx.res$fitted.values, type = "l")

plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, rdcx.res$fitted.values)))
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], rdcx.res$fitted.values, col = "red")

timestamp <- generate.timestamp()
file <- "~/current_projects/gflow/tests/manual/pics/"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=6, height=6)
op <- par(mar=c(2.25,2.25,0.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, rdcx.res$fitted.values)), las = 1)
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], rdcx.res$fitted.values, col = "red")
par(op)
dev.off()
system(paste0("open ",file))
## "~/current_projects/gflow/tests/manual/pics/2025_10_07_180302.pdf"


dcx.res <- fit.knn.riem.graph.regression(
    test.data$X,
    test.data$y,
    k = 3,
    n.eigenpairs = 10,
    filter.type = "heat_kernel", # "tikhonov", "cubic_spline"),
    use.counting.measure = TRUE,
    density.normalization = 0,
    t.diffusion = 0,
    density.uniform.pull = 0,
    response.penalty.exp = 1.0,
    epsilon.y = 1e-4,
    epsilon.rho = 1e-4,
    max.iterations = 10,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.5,
    test.stage = -1
)

plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, rdcx.res$fitted.values)))
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], dcx.res$fitted.values, col = "red")

timestamp <- generate.timestamp()
file <- "~/current_projects/gflow/tests/manual/pics/10itrs_"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=6, height=6)
op <- par(mar=c(2.25,2.25,0.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, rdcx.res$fitted.values)), las = 1)
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], dcx.res$fitted.values, col = "red")
par(op)
dev.off()
system(paste0("open ",file))
## "~/current_projects/gflow/tests/manual/pics/10itrs_2025_10_07_181659.pdf"

dcx.res <- fit.knn.riem.graph.regression(
    test.data$X,
    test.data$y,
    k = 3,
    n.eigenpairs = 10,
    max.iterations = 50,
    filter.type = "heat_kernel", # "tikhonov", "cubic_spline"),
    use.counting.measure = TRUE,
    density.normalization = 0,
    t.diffusion = 0,
    density.uniform.pull = 0,
    response.penalty.exp = 1.0,
    epsilon.y = 1e-4,
    epsilon.rho = 1e-4,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.5,
    test.stage = -1
)

plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, rdcx.res$fitted.values)))
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], dcx.res$fitted.values, col = "red")

timestamp <- generate.timestamp()
file <- "~/current_projects/gflow/tests/manual/pics/50itrs_"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=6, height=6)
op <- par(mar=c(2.25,2.25,0.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, rdcx.res$fitted.values)), las = 1)
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], dcx.res$fitted.values, col = "red")
par(op)
dev.off()
system(paste0("open ",file))
## "~/current_projects/gflow/tests/manual/pics/50itrs_2025_10_07_181824.pdf"



# Compare with standalone ikNN construction
cat("\n--- Comparing with standalone ikNN graph ---\n")
iknn.res <- create.single.iknn.graph(test.data$X, k = 3)
cat("\nikNN summary:\n")
summary(iknn.res)

# Validation checks
cat("\n--- Validation Checks ---\n")

# Check 1: Number of vertices
n_vertices_rdcx <- length(rdcx.res$vertices)  # Adjust based on actual structure
n_vertices_iknn <- length(iknn.res$vertices)
cat("Vertices: rdcx =", n_vertices_rdcx, ", iknn =", n_vertices_iknn, "\n")

# Check 2: Number of edges (after pruning)
n_edges_rdcx <- nrow(rdcx.res$edges)  # Adjust based on actual structure
n_edges_iknn <- nrow(iknn.res$edges)
cat("Edges (pruned): rdcx =", n_edges_rdcx, ", iknn =", n_edges_iknn, "\n")

# Check 3: Edge pruning ratio
initial_edges <- choose(3, 2) * nrow(test.data$X)  # Approximate
pruning_ratio <- 1 - (n_edges_rdcx / initial_edges)
cat("Pruning ratio: ", round(pruning_ratio * 100, 1), "%\n", sep = "")

# Check 4: Basic structure integrity
cat("\nStructure checks:\n")
cat("  - Vertices exist:", !is.null(rdcx.res$vertices), "\n")
cat("  - Edges exist:", !is.null(rdcx.res$edges), "\n")
cat("  - Densities exist:", !is.null(rdcx.res$densities), "\n")  # Adjust field name

# Visual comparison of graphs
if (require(igraph, quietly = TRUE)) {
  par(mfrow = c(1, 2))

  # Plot rdcx graph (you'll need to extract edge list)
  # plot_graph(rdcx.res, main = "RDCX Initialization")

  # Plot iknn graph
  # plot_graph(iknn.res, main = "iKNN Graph")

  par(mfrow = c(1, 1))
}

cat("\n=== STAGE 0 TEST COMPLETE ===\n\n")

# Save results for next test
save(test.data, rdcx.res, iknn.res,
     file = "stage0_results.RData")
cat("Results saved to stage0_results.RData\n")
