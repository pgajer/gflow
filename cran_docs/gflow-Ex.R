### Name: adaptive_mean_shift_gfa
### Title: Fully Adaptive Mean Shift with Gradient Field Averaging
### Aliases: adaptive_mean_shift_gfa

### ** Examples

set.seed(1)
X <- matrix(rnorm(200), ncol = 2)
out <- adaptive_mean_shift_gfa(
  X, k = 10, density_k = 10, n_steps = 5, initial_step_size = 0.1
)
length(out$X_traj)



### Name: adaptive.uggmalo
### Title: Adaptive Uniform Grid Graph Model-Averaged Local Linear
###   Regression
### Aliases: adaptive.uggmalo

### ** Examples

## Not run: 
##D # Create a simple graph with 3 vertices
##D adj.list <- list(c(2), c(1, 3), c(2))
##D weight.list <- list(c(1), c(1, 1), c(1))
##D y <- c(1, 2, 3)
##D 
##D # Run basic analysis
##D result <- adaptive.uggmalo(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   min.path.size = 2,
##D   n.grid.vertices = 5,
##D   n.bws = 10
##D )
##D 
##D # Run with bootstrap confidence intervals
##D result.boot <- adaptive.uggmalo(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   min.path.size = 2,
##D   n.grid.vertices = 5,
##D   n.bws = 10,
##D   n.bb = 100
##D )
## End(Not run)




### Name: add.grad.ED.arrows
### Title: Add Gradient Arrows to 3D Plot
### Aliases: add.grad.ED.arrows

### ** Examples

## Not run: 
##D S <- matrix(rnorm(30), ncol = 3)
##D rownames(S) <- paste0("Point", 1:10)
##D grad.ED <- matrix(rnorm(30) * 0.1, ncol = 3)
##D rownames(grad.ED) <- rownames(S)
##D 
##D plot3D.plain(S)
##D add.grad.ED.arrows(c("Point1", "Point5"), S, grad.ED, C = 2.5)
## End(Not run)




### Name: agemalo
### Title: Adaptive Graph Geodesic Model-Averaged Local Linear Regression
### Aliases: agemalo

### ** Examples

## Not run: 
##D # Create a simple graph with 3 vertices
##D adj.list <- list(c(2), c(1, 3), c(2))
##D weight.list <- list(c(1), c(1, 1), c(1))
##D y <- c(1, 2, 3)
##D 
##D # Run basic analysis
##D result <- agemalo(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   min.path.size = 2,
##D   n.packing.vertices = 5,
##D   n.bws = 10
##D )
##D 
##D # Run with bootstrap confidence intervals
##D result.boot <- agemalo(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   min.path.size = 2,
##D   n.packing.vertices = 5,
##D   n.bws = 10,
##D   n.bb = 100
##D )
## End(Not run)




### Name: amagelo
### Title: Adaptive Model Averaged GEodesic LOcal linear smoothing
### Aliases: amagelo

### ** Examples

## Not run: 
##D # Simulate data with smooth trend and noise
##D x <- seq(0, 10, length.out = 200)
##D y <- sin(x) + 0.2 * rnorm(length(x))
##D 
##D # Apply AMAGELO smoothing
##D result <- amagelo(x, y, grid_size = 100)
##D 
##D # Plot results
##D plot(x, y, pch = 16, col = "gray")
##D lines(x[result$order], result$predictions, col = "red", lwd = 2)
##D 
##D # Examine local extrema
##D extrema <- result$local_extrema
##D points(extrema[,"x"], extrema[,"y"],
##D        pch = ifelse(extrema[,"is_max"] == 1, 24, 25),
##D        bg = ifelse(extrema[,"is_max"] == 1, "red", "blue"),
##D        cex = 2 * extrema[,"range_rel_depth"])
## End(Not run)




### Name: analyze.categorical.proportions
### Title: Analyze Categorical Proportions with Mixed Effects Support
### Aliases: analyze.categorical.proportions

### ** Examples

## Not run: 
##D # Simple analysis without repeated measures
##D set.seed(123)
##D cst <- factor(sample(c("I", "II", "III", "IV"), 200, replace = TRUE))
##D outcome <- rbinom(200, 1, c(0.1, 0.3, 0.2, 0.4)[as.numeric(cst)])
##D 
##D result <- analyze.categorical.proportions(
##D   x = cst,
##D   y = outcome,
##D   pos.label = "Case",
##D   neg.label = "Control"
##D )
##D 
##D # Analysis with repeated measures
##D subjects <- factor(rep(1:50, each = 4))
##D cst_repeated <- factor(sample(c("I", "II", "III", "IV"), 200, replace = TRUE))
##D outcome_repeated <- rbinom(200, 1, 0.3)
##D 
##D result_mixed <- analyze.categorical.proportions(
##D   x = cst_repeated,
##D   y = outcome_repeated,
##D   subj.ids = subjects
##D )
## End(Not run)




### Name: analyze.harmonic.extensions
### Title: Analyze Harmonic Extensions for Local Extrema
### Aliases: analyze.harmonic.extensions

### ** Examples

## Not run: 
##D result <- analyze.harmonic.extensions(
##D   graph$adj_list,
##D   graph$weight_list,
##D   y,
##D   gsf_res$predictions,
##D   gsf_b_cx$extrema_df,
##D   hop_offset = 2,
##D   min_hop_idx = 1,
##D   max_hop_idx = 5
##D )
## End(Not run)




### Name: analyze.hierarchical.differences
### Title: Perform Hierarchical Bayesian Analysis of Method Differences
### Aliases: analyze.hierarchical.differences

### ** Examples

## Not run: 
##D # Example: Compare three methods with hierarchical model
##D set.seed(123)
##D bb.integrals <- list(
##D   method_A = rnorm(100, mean = 0.5, sd = 0.1),
##D   method_B = rnorm(100, mean = 0.52, sd = 0.12),
##D   method_C = rnorm(100, mean = 0.48, sd = 0.08)
##D )
##D 
##D # Fit hierarchical model
##D fit <- analyze.hierarchical.differences(bb.integrals)
##D 
##D # Check convergence
##D rstan::check_hmc_diagnostics(fit)
##D 
##D # Extract posterior summaries
##D print(fit, pars = c("mu", "sigma", "tau"))
##D 
##D # Get pairwise probabilities
##D prob_matrix <- rstan::extract(fit, "prob_diff")$prob_diff
##D apply(prob_matrix, c(2,3), mean)  # Mean probabilities
## End(Not run)



### Name: analyze.weighted.cst
### Title: Time-Weighted CST Analysis for Binary Outcomes
### Aliases: analyze.weighted.cst

### ** Examples

## Not run: 
##D # Example with gestational age time points
##D set.seed(123)
##D 
##D # Generate data for 20 subjects
##D subjects <- rep(1:20, each = 5)
##D weeks <- rep(c(12, 16, 20, 28, 32), 20)  # Gestational weeks
##D 
##D # Simulate CSTs with some subjects more stable than others
##D csts <- character(100)
##D for (i in 1:20) {
##D   if (i <= 10) {  # Stable subjects
##D     csts[(i-1)*5 + 1:5] <- sample(c("I", "III"), 5, replace = TRUE, prob = c(0.8, 0.2))
##D   } else {  # Variable subjects
##D     csts[(i-1)*5 + 1:5] <- sample(c("I", "III", "IV"), 5, replace = TRUE)
##D   }
##D }
##D 
##D # Binary outcome
##D outcomes <- rep(c(rep(0, 10), rep(1, 10)), each = 5)
##D 
##D # Analyze with time weighting
##D result <- analyze.weighted.cst(
##D   x = csts,
##D   y = outcomes,
##D   subj.ids = subjects,
##D   time.points = weeks
##D )
##D 
##D # View results
##D print(result$summary)
##D 
##D # Visualize distributions
##D library(ggplot2)
##D prop_df <- as.data.frame(result$weighted_proportions)
##D prop_df$outcome <- result$subject_data$outcome
##D prop_df$subject <- rownames(prop_df)
##D 
##D # Plot CST I proportions by outcome
##D ggplot(prop_df, aes(x = outcome, y = I)) +
##D   geom_boxplot() +
##D   geom_point(position = position_jitter(width = 0.1)) +
##D   labs(y = "Proportion of time in CST I",
##D        title = "Time-weighted CST proportions by outcome")
## End(Not run)




### Name: angle.3D
### Title: Computes the angle between two 3D vectors in radians
### Aliases: angle.3D

### ** Examples

## Not run: 
##D v <- c(0,1,3)
##D w <- c(2,3,4)
##D angle <- angle.3D(v, w)
## End(Not run)



### Name: angle.alignment
### Title: Align Estimated Angles with True Angles
### Aliases: angle.alignment

### ** Examples

# Generate some true angles around a circle
true.angles <- seq(0, 2*pi - 0.1, length.out = 20)

# Create estimated angles with an arbitrary shift and some noise
est.angles <- (true.angles + 1.5) %% (2*pi) + rnorm(20, 0, 0.1)

# Align the estimated angles
aligned <- angle.alignment(true.angles, est.angles)

# Plot results
plot(true.angles, ylim=c(-0.5, 2*pi+0.5), pch=19, col="blue",
     main="Angle Alignment")
points(est.angles, pch=19, col="red")
points(aligned$angles, pch=19, col="green")
legend("topright", legend=c("True", "Estimated", "Aligned"),
       col=c("blue", "red", "green"), pch=19)




### Name: apply.harmonic.extension
### Title: Apply Harmonic Extension to Local Subgraph
### Aliases: apply.harmonic.extension

### ** Examples

## Not run: 
##D # Assumes subgraph was created by extract.extrema.subgraphs
##D result <- apply.harmonic.extension(subgraph, method = "harmonic_eigen")
## End(Not run)




### Name: basins.merge
### Title: Merge Two Basins in a Gradient Flow Complex
### Aliases: basins.merge

### ** Examples

## Not run: 
##D # Merge descending basins M2 and M3, with M2 absorbing M3
##D merged_flow <- basins.merge(flow, "M2", "M3")
##D 
##D # Check the updated structure
##D summary(merged_flow)
## End(Not run)




### Name: bayes.factors.symbols
### Title: Format Matrix of Bayes Factors into Simplified Notation
### Aliases: bayes.factors.symbols

### ** Examples

# Create example Bayes factor matrix
bf.matrix <- matrix(c(1, 15.2, 0.05, 1200,
                      0.066, 1, 0.002, 8.5,
                      20.1, 500, 1, 2.1,
                      0.0008, 0.118, 0.476, 1),
                    nrow = 4, byrow = TRUE)
rownames(bf.matrix) <- colnames(bf.matrix) <- c("A", "B", "C", "D")

# Format with default thresholds
bayes.factors.symbols(bf.matrix)

# Format with custom thresholds
bayes.factors.symbols(bf.matrix,
                     thresholds = c(strong = 100, moderate = 10, weak = 3),
                     digits = 1)

# Use in knitr/RMarkdown for LaTeX output
# knitr::kable(bayes.factors.symbols(bf.matrix), escape = FALSE)



### Name: bbmwd.over.hHN.graphs
### Title: Calculate Bayesian Bootstrap Mean Wasserstein Distance over
###   k-Hop Neighbor Graphs
### Aliases: bbmwd.over.hHN.graphs

### ** Examples

## Not run: 
##D # Assuming IkNN.graphs is pre-computed
##D set.seed(123)
##D n <- 100
##D y <- sample(0:1, n, replace = TRUE)
##D 
##D # Run BBMWD analysis over a subset of k and hop values
##D results <- bbmwd.over.hHN.graphs(
##D   y = y,
##D   IkNN.graphs = IkNN.graphs,
##D   k.values = 10:15,
##D   hop.values = 2:5,
##D   n.BB = 50
##D )
##D 
##D # Access specific result
##D bbmwd_k10_h3 <- results[[10]][[3]]$bbmwd
## End(Not run)




### Name: bbox.hcube.tiling
### Title: Creates a hypercube tiling of a bounding box
### Aliases: bbox.hcube.tiling

### ** Examples

## Not run: 
##D w <- 1
##D L <- c(0, 0)
##D R <- c(2, 3)
##D sub_boxes <- bbox.hcube.tiling(w, L, R)
## End(Not run)




### Name: bin.segments3d
### Title: Add 3D Line Segments for Binary Variable
### Aliases: bin.segments3d

### ** Examples

## Not run: 
##D X <- matrix(rnorm(300), ncol = 3)
##D y <- sample(0:1, 100, replace = TRUE)
##D names(y) <- rownames(X) <- paste0("Sample", 1:100)
##D 
##D plot3D.plain(X)
##D bin.segments3d(X, y, offset = c(0, 0, 0.1))
## End(Not run)




### Name: boa.overlap
### Title: Calculate Overlap Coefficients Between Basins of Attraction
### Aliases: boa.overlap

### ** Examples

ids <- paste0("id", 1:500)
x <- c(rep(1, 200), rep(2, 300))
names(x) <- ids
y <- c(rep(1, 150), rep(2, 350))
names(y) <- ids
x.labels <- c("1" = "Peak A", "2" = "Peak B")
y.labels <- c("1" = "Max 1", "2" = "Max 2")
boa.overlap(x, y, min.BoA.size = 100, x.labels = x.labels, y.labels = y.labels)




### Name: box.tiling.of.box
### Title: Creates box tiling of a box
### Aliases: box.tiling.of.box

### ** Examples

## Not run: 
##D w <- 1
##D box <- list()
##D box$L <- c(0, 0)
##D box$R <- c(2, 3)
##D box$i <- 1
##D sub_boxes <- box.tiling.of.box(w, box)
## End(Not run)




### Name: box.tiling
### Title: Creates a box tiling of the bounding box of some set
### Aliases: box.tiling

### ** Examples

## Not run: 
##D   X <- matrix(runif(100), ncol = 2)
##D   boxes <- box.tiling(X, n = 50, eps = 0.05, p.exp = 0.05)
## End(Not run)



### Name: boxcox.mle
### Title: Box-Cox MLE for the power parameter lambda
### Aliases: boxcox.mle

### ** Examples

set.seed(1)
n <- 60
x <- runif(n)
y <- exp(1 + 2 * x + rnorm(n, sd = 0.2))  # positive response
d <- data.frame(y = y, x = x)

fit <- boxcox.mle(y ~ x, data = d)
fit$lambda
head(fit$loglik)
fit$ci95

# Transform y with the estimated lambda
y.bc <- boxcox.transform(d$y, fit$lambda)




### Name: boxcox.transform
### Title: Box-Cox power transform
### Aliases: boxcox.transform

### ** Examples

y <- rexp(50, rate = 2)          # positive data
boxcox.transform(y, lambda = 0)  # log-transform
boxcox.transform(y, lambda = 0.5)




### Name: bump.fn
### Title: Bump Function
### Aliases: bump.fn

### ** Examples

## Not run: 
##D x <- seq(-2, 2, length.out = 100)
##D plot(x, bump.fn(x), type = "l")
## End(Not run)



### Name: calculate.edit.distances
### Title: Calculate Edit Distances Between Sequential Graphs
### Aliases: calculate.edit.distances

### ** Examples

## Not run: 
##D # Assuming graph.data has been loaded
##D k.vals <- c(5, 10, 15, 20)
##D distances <- calculate.edit.distances(k.vals, graph.data, offset = 1)
##D plot(distances$k.values, distances$distances, type = "b",
##D      xlab = "k", ylab = "Edit Distance")
## End(Not run)




### Name: centroid.shift
### Title: Subtracts the centroid from each point of the given dataset
### Aliases: centroid.shift

### ** Examples

# Example 1: Auto-compute centroid
X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
centroid.shift(X)

# Example 2: Provide custom centroid
X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
custom_centroid <- c(2, 5)
centroid.shift(X, centroid = custom_centroid)



### Name: circle.plot
### Title: Create Circle Plot of Matrix Data
### Aliases: circle.plot

### ** Examples

## Not run: 
##D X <- matrix(rnorm(50), ncol = 5)
##D colnames(X) <- paste0("Var", 1:5)
##D circle.plot(X)
## End(Not run)




### Name: circular.synthetic.mixture.of.gaussians
### Title: Create a Synthetic Mixture of Circular Gaussians
### Aliases: circular.synthetic.mixture.of.gaussians

### ** Examples

# Generate a simple circular Gaussian mixture
x <- seq(0, 2*pi, length.out = 100)
y <- circular.synthetic.mixture.of.gaussians(
  x,
  x.knot = c(0, pi, 1.5*pi),
  y.knot = c(5, 8, 2.5),
  sd.knot = 0.3
)

# Plot the result
plot(x, y, type = "l", xlab = "Angle (radians)", ylab = "Value")




### Name: classifier.based.divergence
### Title: Estimate Divergence Using a Classifier-based Approach
### Aliases: classifier.based.divergence

### ** Examples

X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- classifier.based.divergence(X, Y)
print(result)




### Name: clusters.reorder
### Title: Reorder Cluster IDs by Size
### Aliases: clusters.reorder

### ** Examples

# Numeric clusters
cltr <- c(1, 2, 2, 3, 3, 3)
clusters.reorder(cltr)
# Returns: [1] 3 2 2 1 1 1

# Character clusters
cltr_char <- c("A", "B", "B", "C", "C", "C")
clusters.reorder(cltr_char)
# Returns: [1] "C" "B" "B" "A" "A" "A"

# With names preserved
named_cltr <- c(a = 1, b = 2, c = 2, d = 3, e = 3, f = 3)
clusters.reorder(named_cltr)




### Name: compare.adj.lists
### Title: Compare Two Adjacency Lists
### Aliases: compare.adj.lists

### ** Examples

adj.list1 <- list(c(2, 3), c(1, 3), c(1, 2))
adj.list2 <- list(c(3, 2), c(3, 1), c(2, 1))
compare.adj.lists(adj.list1, adj.list2) # returns TRUE

adj.list1 <- list(c(2, 3, 4), c(1, 3), c(1, 2), c(1))
adj.list2 <- list(c(3, 2), c(3, 1), c(2, 1), c(1))
compare.adj.lists(adj.list1, adj.list2) # returns FALSE




### Name: compare.graph.vs.nerve.cx.filtering
### Title: Compare Graph vs. Nerve Complex Spectral Filtering
### Aliases: compare.graph.vs.nerve.cx.filtering

### ** Examples

## Not run: 
##D # Generate test data
##D set.seed(123)
##D coords <- matrix(runif(200), ncol = 2)
##D f_true <- function(x) sin(2*pi*x[1]) * cos(2*pi*x[2])
##D y_true <- apply(coords, 1, f_true)
##D y_noisy <- y_true + rnorm(length(y_true), 0, 0.2)
##D 
##D # Create nerve complex
##D complex <- create.nerve.complex(coords, k = 8, max.dim = 2)
##D complex <- set.complex.function.values(complex, y_noisy)
##D 
##D # Compare filtering approaches
##D comparison <- compare.graph.vs.nerve.cx.filtering(
##D   complex, y_noisy, y_true,
##D   dim.weights.complex = c(1.0, 0.5, 0.25),
##D   verbose = TRUE
##D )
##D 
##D # Print improvement
##D cat("GCV improvement:", comparison$gcv_improvement_pct, "%\n")
##D cat("MSE improvement:", comparison$mse_improvement_pct, "%\n")
## End(Not run)




### Name: compare.harmonic.methods
### Title: Compare Multiple Harmonic Extension Methods
### Aliases: compare.harmonic.methods

### ** Examples

## Not run: 
##D # Assumes subgraph was created by extract.extrema.subgraphs
##D comparison <- compare.harmonic.methods(
##D   subgraph,
##D   methods = c("harmonic_iterative", "harmonic_eigen")
##D )
## End(Not run)




### Name: compare.paths
### Title: Compare Paths Across Different Hop Limits
### Aliases: compare.paths

### ** Examples

## Not run: 
##D pgs <- create.path.graph.series(graph, edge.lengths, h.values = c(1, 2, 3))
##D compare.paths(pgs, from = 1, to = 3)
## End(Not run)




### Name: compute.bayesian.effects
### Title: Compute Standardized Effect Sizes with Bayesian Inference
### Aliases: compute.bayesian.effects

### ** Examples

# Create example data
bb.integrals <- list(
  method1 = rnorm(100, mean = 0, sd = 1),
  method2 = rnorm(100, mean = 0.5, sd = 1),
  method3 = rnorm(100, mean = -0.3, sd = 1.2)
)

results <- compute.bayesian.effects(bb.integrals)

# View results for specific method
results["method2", ]

# Compare probabilities
results[, c("prob_smaller", "prob_larger", "prob_practical_equiv")]




### Name: compute.diff.profiles
### Title: Create Higher-Order Difference Profiles
### Aliases: compute.diff.profiles

### ** Examples

## Not run: 
##D # Generate example curve
##D x <- seq(0, 1, length.out = 50)
##D Eyg <- sin(2*pi*x)
##D 
##D # Get profiles up to 3rd order
##D profiles <- compute.diff.profiles(Eyg, k = 3)
##D 
##D # The result contains:
##D # - Original values (50 elements)
##D # - First differences (49 elements)
##D # - Second differences (48 elements)
##D # - Third differences (47 elements)
##D # Total: 194 elements
## End(Not run)




### Name: compute.gradient.trajectory
### Title: Compute Gradient Trajectory
### Aliases: compute.gradient.trajectory

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) -((x-0.5)^2 + (y-0.5)^2)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D traj <- compute.gradient.trajectory(15, 15, f.grid)
##D str(traj)
## End(Not run)



### Name: compute.graph.diameter
### Title: Compute the Diameter of a Weighted Undirected Graph
### Aliases: compute.graph.diameter

### ** Examples

# Example with a simple graph
adj.list <- list(c(2,3), c(1,3), c(1,2))
weight.list <- list(c(1,2), c(1,3), c(2,3))
result <- compute.graph.diameter(adj.list, weight.list)
print(result$message)




### Name: compute.graph.distance
### Title: Calculate Shortest Path Distance Between Two Vertices in a Graph
### Aliases: compute.graph.distance

### ** Examples

## Not run: 
##D # Using star_obj
##D star.obj <- generate.star.dataset(n.points = 10, n.arms = 3)
##D dist <- compute.graph.distance(star.obj = star.obj, i = 1, j = 5)
##D 
##D # Using adjacency list and edge lengths directly
##D # Assuming adj_list and edge_lengths are defined
##D dist <- compute.graph.distance(i = 1, j = 5,
##D                               adj.list = adj_list,
##D                               edge.lengths = edge_lengths)
## End(Not run)



### Name: compute.graph.gradient.flow
### Title: Compute Gradient Flow Information for Vertices in a Graph
### Aliases: compute.graph.gradient.flow

### ** Examples

## Not run: 
##D # Create example graph
##D adj.list <- list(c(2), c(1,3), c(2,4,5), c(3), c(3))
##D y <- c(0.2, 0.5, 0.8, 1.0, 0.1)
##D 
##D # Define extrema
##D lmax.list <- list(list(lmax=4, vertices=c(4), label="peak"))
##D lmin.list <- list(list(lmin=5, vertices=c(5), label="valley"))
##D 
##D # Compute gradient flow
##D flow_info <- compute.graph.gradient.flow(adj.list, NULL, y,
##D                                          lmax.list, lmin.list)
##D 
##D # Examine results
##D print(flow_info$basin_stats)
## End(Not run)




### Name: compute.local.distance.fidelity
### Title: Assess Fidelity of Graph-Based Geodesic Distances to Euclidean
###   Geometry
### Aliases: compute.local.distance.fidelity

### ** Examples

## Not run: 
##D set.seed(1)
##D X <- matrix(rnorm(300 * 2), ncol = 2)
##D g <- build_graph(X)  # your function to generate adj.list and weight.list
##D taus <- seq(0.05, 0.2, by = 0.01)
##D fidelity <- compute.local.distance.fidelity(X, g$adj.list, g$weight.list, taus)
## End(Not run)




### Name: compute.morse.smale.cells
### Title: Compute Morse-Smale Cells
### Aliases: compute.morse.smale.cells

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D cells <- compute.morse.smale.cells(f.grid)
##D image(cells)
## End(Not run)



### Name: compute.pairwise.bayes.factors
### Title: Compute Pairwise Bayes Factors for Method Comparisons
### Aliases: compute.pairwise.bayes.factors

### ** Examples

# Example: Compare three methods
set.seed(123)
bb.integrals <- list(
  method_A = rnorm(1000, mean = 0.5, sd = 0.1),
  method_B = rnorm(1000, mean = 0.52, sd = 0.1),
  method_C = rnorm(1000, mean = 0.48, sd = 0.1)
)

bf_matrix <- compute.pairwise.bayes.factors(bb.integrals)
print(bf_matrix)

# Interpret results
# bf_matrix[1,2] > 1 suggests evidence for method_A over method_B
# bf_matrix[2,3] > 1 suggests evidence for method_B over method_C



### Name: compute.persistence
### Title: Compute Persistent Homology for Graphs
### Aliases: compute.persistence

### ** Examples

# Simple chain graph
adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
y <- c(10, 7, 9, 2, 5)

# Compute persistence
result <- compute.persistence(adj.list, y)
print(result$persistence)




### Name: construct.graph.gradient.flow
### Title: Graph Gradient Flow Analysis
### Aliases: construct.graph.gradient.flow

### ** Examples

## Not run: 
##D # Create a simple grid graph
##D n <- 10
##D adj.list <- list()
##D weight.list <- list()
##D 
##D # Build grid adjacency (simplified example)
##D for(i in 1:(n*n)) {
##D   neighbors <- integer()
##D   weights <- numeric()
##D 
##D   # Add horizontal neighbors
##D   if(i %% n != 1) {
##D     neighbors <- c(neighbors, i-1)
##D     weights <- c(weights, 1.0)
##D   }
##D   if(i %% n != 0) {
##D     neighbors <- c(neighbors, i+1)
##D     weights <- c(weights, 1.0)
##D   }
##D 
##D   adj.list[[i]] <- neighbors
##D   weight.list[[i]] <- weights
##D }
##D 
##D # Define a function on the graph (e.g., a Gaussian peak)
##D coords <- expand.grid(x = 1:n/n, y = 1:n/n)
##D z <- exp(-10 * ((coords$x - 0.5)^2 + (coords$y - 0.5)^2))
##D 
##D # Compute gradient flow
##D flow <- construct.graph.gradient.flow(
##D   adj.list,
##D   weight.list,
##D   z,
##D   scale = 2.0,
##D   with.trajectories = TRUE
##D )
##D 
##D # Examine results
##D print(flow$local_extrema)
##D summary(flow)
## End(Not run)




### Name: cont.segments3d
### Title: Add 3D Line Segments Color-Coded by Continuous Variable
### Aliases: cont.segments3d

### ** Examples

## Not run: 
##D X <- matrix(rnorm(300), ncol = 3)
##D y <- runif(100)
##D names(y) <- rownames(X) <- paste0("Sample", 1:100)
##D 
##D plot3D.plain(X)
##D cont.segments3d(X, y, offset = c(0, 0, 0.1))
## End(Not run)




### Name: convert.adjacency.list.to.adjacency.matrix
### Title: Converts a graph adjacency list to an adjacency matrix
### Aliases: convert.adjacency.list.to.adjacency.matrix

### ** Examples

adj.list <- list(c(2, 4), c(1, 3), c(2, 4), c(1, 3))
weights.list <- list(c(2, 4), c(2, 3), c(3, 1), c(4, 1))
convert.adjacency.list.to.adjacency.matrix(adj.list)
convert.adjacency.list.to.adjacency.matrix(adj.list, weights.list)



### Name: convert.adjacency.to.edge.matrix
### Title: Convert Adjacency List to Edge Matrix
### Aliases: convert.adjacency.to.edge.matrix

### ** Examples

adj.list <- list(c(2,3), c(1,3), c(1,2))
result <- convert.adjacency.to.edge.matrix(adj.list)

# With weights
weights.list <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
result_weighted <- convert.adjacency.to.edge.matrix(adj.list, weights.list)




### Name: convert.to.undirected
### Title: Converts a directed graph to an undirected graph
### Aliases: convert.to.undirected

### ** Examples

directed_graph <- list(
  "A" = c("B", "C"),
  "B" = c("C"),
  "C" = c("D"),
  "D" = c("A")
)

undirected_graph <- convert.to.undirected(directed_graph)
print(undirected_graph)




### Name: count.edges
### Title: Count Edges in an Undirected Graph
### Aliases: count.edges

### ** Examples

# Create an adjacency list for a simple undirected graph
# with 3 vertices and 2 edges: (1-2) and (2-3)
adj <- list(
  c(2),    # Vertex 1 is connected to vertex 2
  c(1, 3), # Vertex 2 is connected to vertices 1 and 3
  c(2)     # Vertex 3 is connected to vertex 2
)
count.edges(adj)  # Should return 2




### Name: create.adaptive.tiled.X.grid.xD
### Title: Create an Adaptive Tiled Grid Representation of X
### Aliases: create.adaptive.tiled.X.grid.xD

### ** Examples

## Not run: 
##D X <- matrix(rnorm(1000), ncol = 2)
##D grid_obj <- create.adaptive.tiled.X.grid.xD(X, gSf = 2, gRf = 1.5, verbose = TRUE)
##D plot(grid_obj$X.grid, col = "red", pch = 20)
##D points(X, col = "blue", pch = 1)
## End(Not run)




### Name: create.basin.cx
### Title: Create a Gradient Flow Basin Complex on a Weighted Graph
### Aliases: create.basin.cx

### ** Examples

## Not run: 
##D # Create a graph with adjacency list and weights
##D adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
##D weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
##D 
##D # Define a function on vertices
##D y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
##D 
##D # Create basin complex
##D basin_cx <- create.basin.cx(adj_list, weight_list, y, basin.merge.overlap.thld = 0.15)
##D 
##D # Examine results
##D summary(basin_cx)
##D 
##D # Visualize original vs. simplified function
##D plot(basin_cx, type = "comparison")
##D 
##D # Extract vertices from a specific basin
##D m1_vertices <- get_basin_vertices(basin_cx, "m1")
## End(Not run)




### Name: create.bi.kNN.chain.graph
### Title: Creates a bi-k-NN chain graph
### Aliases: create.bi.kNN.chain.graph

### ** Examples

# Create a 5-chain graph with 10 vertices
graph <- create.bi.kNN.chain.graph(10, 5)

# Create a graph based on x-coordinates
x <- runif(10)
y <- sin(x)
graph.with.coords <- create.bi.kNN.chain.graph(10, 2, x, y)




### Name: create.bipartite.graph
### Title: Create a Bipartite Graph
### Aliases: create.bipartite.graph

### ** Examples

bipartite_graph <- create.bipartite.graph(3, 4)




### Name: create.chain.graph
### Title: Create a Chain Graph Structure
### Aliases: create.chain.graph

### ** Examples

# Create a simple 4-vertex chain
chain1 <- create.chain.graph(4)



### Name: create.cmst.graph
### Title: Construct a Minimal Spanning Tree (MST) Completion Graph
### Aliases: create.cmst.graph

### ** Examples

## Not run: 
##D # Generate sample data
##D set.seed(123)
##D X <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
##D 
##D # Create MST completion graph with default parameters
##D graph <- create.cmst.graph(X)
##D 
##D # Create graph with PCA dimensionality reduction
##D X_high <- matrix(rnorm(100 * 200), nrow = 100, ncol = 200)
##D graph_pca <- create.cmst.graph(X_high, pca.dim = 50, variance.explained = 0.95)
##D 
##D # Print summary
##D print(graph_pca)
##D summary(graph_pca)
## End(Not run)



### Name: create.complete.graph
### Title: Create a Complete Graph
### Aliases: create.complete.graph

### ** Examples

complete_graph <- create.complete.graph(5)




### Name: create.delta1.Delta1.df
### Title: Create Summary Data Frame from First-Order Association Test
###   Results
### Aliases: create.delta1.Delta1.df

### ** Examples

## Not run: 
##D # Run multiple tests
##D results <- list()
##D for (i in 1:10) {
##D   x <- runif(100)
##D   y <- sin(2*pi*x) + rnorm(100, sd = 0.3)
##D   results\code{[[paste0("var", i)]}] <- fassoc.test(x, y, order = 1)
##D }
##D 
##D # Create summary
##D summary_df <- create.delta1.Delta1.df(results)
##D print(summary_df$sign.d1D1.df)
## End(Not run)




### Name: create.distance.plot
### Title: Create Distance Plot
### Aliases: create.distance.plot

### ** Examples

## Not run: 
##D # Create sample data
##D k.vals <- seq(5, 50, by = 5)
##D distances <- 100 / k.vals + rnorm(length(k.vals), sd = 2)
##D 
##D # Basic plot
##D create.distance.plot(k.vals, distances)
##D 
##D # Plot with marked minimum
##D min.k <- k.vals[which.min(distances)]
##D create.distance.plot(k.vals, distances, mark.x = min.k,
##D                     main = "Edit Distance vs k")
## End(Not run)



### Name: create.ED.grid.xD
### Title: Create an Equidistant Grid in x-Dimensional Bounding Box (C
###   Interface)
### Aliases: create.ED.grid.xD

### ** Examples

## Not run: 
##D edge.length <- 1.0
##D lower.bounds <- c(0, 0, 0)
##D upper.bounds <- c(2, 3, 4)
##D grid <- create.ED.grid.xD(edge.length, lower.bounds, upper.bounds, round.up = TRUE)
##D print(grid)
## End(Not run)



### Name: create.empty.graph
### Title: Create an Empty Graph
### Aliases: create.empty.graph

### ** Examples

empty_graph <- create.empty.graph(5)




### Name: create.final.grid
### Title: Creates a Final Grid within Selected Boxes
### Aliases: create.final.grid

### ** Examples

## Not run: 
##D boxes <- create.ED.boxes(0.5, c(0,0), c(1,1))
##D X <- matrix(runif(20), ncol = 2)
##D w <- 0.1
##D epsilon <- 0.5
##D final_grid <- create.final.grid(boxes, X, w, epsilon)
## End(Not run)



### Name: create.gflow.cx
### Title: Create Graph Flow Complex with Harmonic Extension
### Aliases: create.gflow.cx

### ** Examples

## Not run: 
##D # Create a simple triangle graph
##D adj.list <- list(c(2,3), c(1,3), c(1,2))
##D weight.list <- list(c(1,1), c(1,1), c(1,1))
##D 
##D # Function with a spurious maximum at vertex 2
##D y <- c(0.5, 1.0, 0.7)
##D 
##D # Smooth out spurious extrema
##D result <- create.gflow.cx(
##D   adj.list, weight.list, y,
##D   hop.idx.thld = 1,
##D   smoother.type = 0,  # Weighted Mean
##D   verbose = FALSE
##D )
##D 
##D # Check smoothed values
##D print(result$harmonic_predictions)
##D 
##D # More complex example with detailed recording
##D # Create a larger graph (grid)
##D n <- 5
##D adj.list <- vector("list", n*n)
##D weight.list <- vector("list", n*n)
##D 
##D for(i in 1:n) {
##D   for(j in 1:n) {
##D     idx <- (i-1)*n + j
##D     neighbors <- c()
##D     weights <- c()
##D 
##D     # Add edges to grid neighbors
##D     if(i > 1) { neighbors <- c(neighbors, (i-2)*n + j); weights <- c(weights, 1) }
##D     if(i < n) { neighbors <- c(neighbors, i*n + j); weights <- c(weights, 1) }
##D     if(j > 1) { neighbors <- c(neighbors, (i-1)*n + j-1); weights <- c(weights, 1) }
##D     if(j < n) { neighbors <- c(neighbors, (i-1)*n + j+1); weights <- c(weights, 1) }
##D 
##D     adj.list[[idx]] <- neighbors
##D     weight.list[[idx]] <- weights
##D   }
##D }
##D 
##D # Create function with multiple extrema
##D y <- sin(seq(0, 2*pi, length.out = n*n)) + 0.1*rnorm(n*n)
##D 
##D # Apply smoothing with detailed recording
##D result <- create.gflow.cx(
##D   adj.list, weight.list, y,
##D   hop.idx.thld = 2,
##D   smoother.type = 2,  # Harmonic Eigen
##D   detailed.recording = TRUE,
##D   verbose = TRUE
##D )
##D 
##D # Examine extrema
##D print(result$extrema_df)
## End(Not run)




### Name: create.grid.graph
### Title: Create a Refined Graph with Uniformly Spaced Grid Vertices
### Aliases: create.grid.graph

### ** Examples

# Create a simple path graph with 3 vertices
adj <- list(c(2L), c(1L, 3L), c(2L))
weights <- list(c(1.0), c(1.0, 2.0), c(2.0))

# Add approximately 5 grid vertices
result <- create.grid.graph(adj, weights, grid.size = 5)

# Examine the results
cat("Original graph had", length(adj), "vertices\\n")
cat("Refined graph has", length(result$adj.list), "vertices\\n")
cat("Grid vertices added:", length(result$grid.vertices), "\\n")
cat("Indices of grid vertices:", result$grid.vertices, "\\n")

## Not run: 
##D # Create a more complex graph (a cycle with varying edge weights)
##D n <- 6
##D adj <- lapply(1:n, function(i) c(ifelse(i == 1, n, i - 1),
##D                                  ifelse(i == n, 1, i + 1)))
##D weights <- lapply(1:n, function(i) runif(2, 0.5, 2.0))
##D 
##D # Refine with more grid vertices
##D result <- create.grid.graph(adj, weights, grid.size = 20)
##D 
##D g <- igraph::graph_from_adj_list(result$adj.list, mode = "undirected")
##D plot(g, vertex.color = ifelse(1:length(result$adj.list) %in%
##D                                 result$grid.vertices, "red", "blue"))
## End(Not run)




### Name: create.grid
### Title: Create Uniform 2D Grid
### Aliases: create.grid

### ** Examples

grid <- create.grid(50)
str(grid)



### Name: create.hHN.graph
### Title: Create a k-hop Neighborhood (hHN) Graph
### Aliases: create.hHN.graph

### ** Examples

# Create a simple 4-vertex graph
graph <- list(
  c(2, 3),      # Vertex 1 connects to vertices 2 and 3
  c(1, 3, 4),   # Vertex 2 connects to vertices 1, 3, and 4
  c(1, 2, 4),   # Vertex 3 connects to vertices 1, 2, and 4
  c(2, 3)       # Vertex 4 connects to vertices 2 and 3
)

edge.lengths <- list(
  c(1, 2),      # Edge weights from vertex 1
  c(1, 1, 3),   # Edge weights from vertex 2
  c(2, 1, 1),   # Edge weights from vertex 3
  c(3, 1)       # Edge weights from vertex 4
)

# Generate the 2-hop neighborhood graph
h <- 2
khn.graph <- create.hHN.graph(graph, edge.lengths, h)

# Access the results
khn.graph$adj_list   # Adjacency list of the hHN graph
khn.graph$dist_list  # Shortest path distances in the hHN graph




### Name: create.hop.nbhd.extrema.df
### Title: Create Data Frame of Extrema Information from Hop Neighborhoods
### Aliases: create.hop.nbhd.extrema.df

### ** Examples

## Not run: 
##D # Create example graph and compute hop neighborhoods
##D adj.list <- list(c(2,3), c(1,3,4), c(1,2,4), c(2,3))
##D weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1,1))
##D y <- c(0, 1, 0.5, 2)  # Minima at vertices 1,3 and maxima at vertices 2,4
##D 
##D # Using with create.gflow.cx
##D result <- create.gflow.cx(adj.list, weight.list, y, verbose = FALSE)
##D extrema_df <- create.hop.nbhd.extrema.df(result)
##D print(extrema_df)
##D 
##D # Exclude spurious extrema
##D significant_extrema <- create.hop.nbhd.extrema.df(
##D   result,
##D   include_spurious = FALSE,
##D   threshold = 2
##D )
##D print(significant_extrema)
## End(Not run)



### Name: create.iknn.graphs
### Title: Create Intersection k-Nearest Neighbor Graphs with Dual Pruning
###   Methods
### Aliases: create.iknn.graphs

### ** Examples

## Not run: 
##D # Generate sample data
##D X <- matrix(rnorm(100 * 5), 100, 5)
##D 
##D # Basic usage
##D result <- create.iknn.graphs(
##D   X, kmin = 3, kmax = 10,
##D   compute.full = FALSE
##D )
##D 
##D # With custom pruning parameters
##D result <- create.iknn.graphs(
##D   X, kmin = 3, kmax = 10,
##D   max.path.edge.ratio.deviation.thld = 0.1,
##D   path.edge.ratio.percentile = 0.5,
##D   compute.full = TRUE,
##D   verbose = TRUE
##D )
##D 
##D # View statistics for each k
##D print(result$k_statistics)
## End(Not run)




### Name: create.intersection.graph
### Title: Create an Intersection Graph
### Aliases: create.intersection.graph

### ** Examples

adj.list <- list(c(1,2), c(0,2,3), c(0,1), c(1))
result <- create.intersection.graph(adj.list, 0.5)
print(result)




### Name: create.latex.table
### Title: Create LaTeX Table from Matrix or Data Frame
### Aliases: create.latex.table

### ** Examples

## Not run: 
##D # Create example data
##D mat <- matrix(1:12, nrow = 3)
##D colnames(mat) <- paste0("Group", 1:4)
##D rownames(mat) <- paste0("Category", 1:3)
##D 
##D # Generate LaTeX table
##D create.latex.table(
##D   data = mat,
##D   file = "output.tex",
##D   label = "results",
##D   caption = "Example results table"
##D )
## End(Not run)




### Name: create.lmin.lmax.contingency.table
### Title: Create Contingency Table of Local Maxima and Minima
###   Relationships
### Aliases: create.lmin.lmax.contingency.table

### ** Examples

## Not run: 
##D result <- create.lmin.lmax.contingency.table(MS.res, lmin.labels, lmax.labels)
##D print(result$contingency_table)
##D cat(result$caption)
## End(Not run)



### Name: create.lmin.lmax.label.indicators
### Title: Create Label Tables and Indicator Vectors for Morse-Smale
###   Complex Critical Points
### Aliases: create.lmin.lmax.label.indicators

### ** Examples

## Not run: 
##D # Given lmin.labels, lmax.labels and state space matrix state.space
##D result <- create.lmin.lmax.label.indicators(lmin.labels, lmax.labels, state.space)
##D print(head(result$lmin.lab.tbl))
##D print(sum(result$lmin.ind)) # Number of local minima
## End(Not run)



### Name: create.lmin.lmax.labels
### Title: Create Labels for Local Maxima and Minima in Morse-Smale Complex
### Aliases: create.lmin.lmax.labels

### ** Examples

## Not run: 
##D labels <- create.lmin.lmax.labels(MS.res, state.space, taxonomy,
##D                                   freq.thld = 100, min.relAb.thld = 0.05)
##D print(labels$lmax.labels)
##D print(labels$lmin.labels)
## End(Not run)



### Name: create.maximal.packing
### Title: Create a Maximal Packing of Vertices in a Graph
### Aliases: create.maximal.packing

### ** Examples

## Not run: 
##D # Create a simple triangle graph
##D adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
##D weight.list <- list(c(1, 1), c(1, 1), c(1, 1))
##D 
##D # Create grid graph with grid size 2
##D result <- create.maximal.packing(adj.list, weight.list, grid.size = 2)
##D 
##D # View the vertices in the maximal packing
##D print(result$grid_vertices)
##D 
##D # Access the graph diameter and packing radius
##D cat("Graph diameter:", result$graph_diameter, "\n")
##D cat("Packing radius:", result$max_packing_radius, "\n")
## End(Not run)




### Name: create.mknn.graph
### Title: Compute a Mutual k-Nearest Neighbor Graph with Weights
### Aliases: create.mknn.graph

### ** Examples

## Not run: 
##D # Generate sample 2D data
##D set.seed(123)
##D X <- matrix(rnorm(100 * 2), ncol = 2)
##D 
##D # Create mutual 5-NN graph
##D graph <- create.mknn.graph(X, k = 5)
##D 
##D # Print basic statistics
##D cat("Number of vertices:", graph$n_vertices, "\n")
##D cat("Number of edges:", graph$n_edges, "\n")
##D 
##D # Examine connections for first vertex
##D cat("Vertex 1 is connected to:", graph$adj_list[[1]], "\n")
##D cat("With distances:", round(graph$weight_list[[1]], 3), "\n")
## End(Not run)




### Name: create.mknn.graphs
### Title: Create Multiple Mutual kNN Graphs with Geometric Pruning
### Aliases: create.mknn.graphs

### ** Examples

## Not run: 
##D # Generate sample data with clusters
##D set.seed(123)
##D n <- 150
##D X <- rbind(
##D   matrix(rnorm(n * 2, mean = 0), ncol = 2),
##D   matrix(rnorm(n * 2, mean = 5), ncol = 2),
##D   matrix(rnorm(n * 2, mean = c(2.5, 5)), ncol = 2)
##D )
##D 
##D # Create MkNN graphs for k from 5 to 15 with pruning
##D result <- create.mknn.graphs(
##D   X,
##D   kmin = 5,
##D   kmax = 15,
##D   max.path.edge.ratio.thld = 1.2,
##D   path.edge.ratio.percentile = 0.5,
##D   compute.full = TRUE,
##D   verbose = TRUE
##D )
##D 
##D # Examine the summary statistics
##D print(result$k_statistics)
##D 
##D # Get the graph for k=10
##D k10_graph <- result$pruned_graphs[["10"]]
##D cat("Graph with k=10 has", k10_graph$n_edges, "edges\n")
## End(Not run)

# High-dimensional example with PCA
## Not run: 
##D # Generate high-dimensional data
##D set.seed(456)
##D X_highdim <- matrix(rnorm(200 * 1000), nrow = 200, ncol = 1000)
##D 
##D # Apply PCA before graph construction
##D result_pca <- create.mknn.graphs(
##D   X_highdim,
##D   kmin = 10,
##D   kmax = 20,
##D   pca.dim = 50,
##D   variance.explained = 0.95,
##D   verbose = TRUE
##D )
##D 
##D # Check PCA information
##D pca_info <- attr(result_pca, "pca")
##D cat("Used", pca_info$n_components, "components explaining",
##D     round(pca_info$variance_explained * 100, 2), "% of variance\n")
## End(Not run)




### Name: create.morse.smale.complex
### Title: Create Morse-Smale Complex Object
### Aliases: create.morse.smale.complex

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D complex <- create.morse.smale.complex(f.grid)
##D str(complex)
## End(Not run)



### Name: create.path.graph
### Title: Create a Path Graph with Limited Hop Distance
### Aliases: create.path.graph

### ** Examples

## Not run: 
##D # Create a simple graph with 3 vertices
##D graph <- list(
##D   c(2),    # Vertex 1 connected to 2
##D   c(1, 3), # Vertex 2 connected to 1 and 3
##D   c(2)     # Vertex 3 connected to 2
##D )
##D 
##D # Define edge lengths
##D edge.lengths <- list(
##D   c(1.0),      # Length of edge 1->2
##D   c(1.0, 2.0), # Lengths of edges 2->1 and 2->3
##D   c(2.0)       # Length of edge 3->2
##D )
##D 
##D # Create path graph with maximum 2 hops
##D pg <- create.path.graph(graph, edge.lengths, h = 2)
##D 
##D # Print the path graph
##D print(pg)
##D 
##D # Get shortest path between vertices 1 and 3
##D path <- get.shortest.path(pg, 1, 3)
## End(Not run)




### Name: create.path.graph.series
### Title: Create a Series of Path Graphs with Different Hop Limits
### Aliases: create.path.graph.series

### ** Examples

## Not run: 
##D # Create a simple graph
##D graph <- list(
##D   c(2),    # Vertex 1 connected to 2
##D   c(1, 3), # Vertex 2 connected to 1 and 3
##D   c(2)     # Vertex 3 connected to 2
##D )
##D edge.lengths <- list(c(1.0), c(1.0, 2.0), c(2.0))
##D 
##D # Create path graphs for hop limits 1, 2, and 3
##D pgs <- create.path.graph.series(graph, edge.lengths, h.values = c(1, 2, 3))
##D 
##D # Access individual path graphs
##D pg.h2 <- pgs[[2]]  # Path graph with h=2
## End(Not run)




### Name: create.plm.graph
### Title: Create a Path Length Matrix Graph Structure
### Aliases: create.plm.graph

### ** Examples

## Not run: 
##D # Create a simple path graph with 3 vertices: 1 -- 2 -- 3
##D graph <- list(c(2), c(1, 3), c(2))
##D edge.lengths <- list(c(1), c(1, 1), c(1))
##D 
##D # Create PLM graph with paths up to length 3
##D plm_graph <- create.plm.graph(graph, edge.lengths, h = 3)
##D 
##D print(plm_graph)
## End(Not run)




### Name: create.pruned.isize.list
### Title: Create Pruned Intersection Size List
### Aliases: create.pruned.isize.list

### ** Examples

## Not run: 
##D # Example usage:
##D adj.list <- list(c(2,3,4), c(1,3), c(1,2,4), c(1,3))
##D isize.list <- list(c(2,1,3), c(2,1), c(1,1,2), c(3,2))
##D pruned.adj.list <- list(c(2,4), c(1), c(4), c(1,3))
##D pruned.isize.list <- create.pruned.isize.list(adj.list, isize.list, pruned.adj.list)
##D print(pruned.isize.list)
## End(Not run)




### Name: create.random.graph
### Title: Create a Random Graph
### Aliases: create.random.graph

### ** Examples

graph <- create.random.graph(100, 4)




### Name: create.single.iknn.graph
### Title: Create a Single Intersection-weighted k-Nearest Neighbors Graph
### Aliases: create.single.iknn.graph

### ** Examples

# Create sample data
set.seed(123)
X <- matrix(rnorm(100), ncol = 2)
result <- create.single.iknn.graph(X, k = 3)
summary(result)



### Name: create.star.graph
### Title: Create a star graph by joining chain graphs
### Aliases: create.star.graph

### ** Examples

star_graph <- create.star.graph(c(3, 4, 2))
print(star_graph)




### Name: create.subgraph
### Title: Create a subgraph from a given graph
### Aliases: create.subgraph

### ** Examples

X <- runif.sphere(20, 2)
graph <- create.single.iknn.graph(X, k = 3, compute.full = TRUE, verbose = FALSE)
graph$dist_list <- graph$weight_list
graph$weight_list <- NULL
subgraph <- create.subgraph(graph, id.indices = c(1:10))
subgraph_sequential <- create.subgraph(graph, id.indices = c(1:10, 15:16),
                                       use.sequential.indices = TRUE)



### Name: create.threshold.distance.graph
### Title: Create Threshold Distance Graph from Distance Matrix
### Aliases: create.threshold.distance.graph

### ** Examples

# Example distance matrix
dist <- matrix(c(
  0.0000000, 0.02834008, 0.05050505, 0.12500000, 0.1086957,
  0.02834008, 0.00000000, 0.88888889, 0.54166667, 1.0000000,
  0.05050505, 0.88888889, 0.00000000, 0.04166667, 0.1086957,
  0.12500000, 0.54166667, 0.04166667, 0.00000000, 1.0000000,
  0.10869565, 1.00000000, 0.10869565, 1.00000000, 0.0000000
), nrow=5, byrow=TRUE)
rownames(dist) <- colnames(dist) <- c("M1", "M2", "M3", "M4", "M5")

# Force symmetry by averaging with transpose
dist <- (dist + t(dist)) / 2

# Create graph with threshold 0.15
graph <- create.threshold.distance.graph(dist, 0.15)




### Name: create.vertex.label.indicators
### Title: Create Label Tables and Indicator Vectors for Vertices
### Aliases: create.vertex.label.indicators

### ** Examples

## Not run: 
##D # Given vertex.labels and state space matrix state.space
##D result <- create.vertex.label.indicators(vertex.labels, state.space)
##D print(head(result$lab.tbl))
##D print(sum(result$ind)) # Number of vertices
## End(Not run)



### Name: create.vertex.labels
### Title: Create Labels for Vertices in State Space
### Aliases: create.vertex.labels

### ** Examples

## Not run: 
##D # Keep only top 5 species in profiles
##D result <- create.vertex.labels(c(1,3,5), state.space, taxonomy,
##D                               min.relAb.thld = 0.05, profile.length = 5)
## End(Not run)



### Name: create.X.grid
### Title: Creates a grid around a state space
### Aliases: create.X.grid

### ** Examples

## Not run: 
##D # Let X be a low-dimensional model of a state space.
##D res <- create.X.grid(X, gSf=5, gRf=5, min.K=10, med.dK.divf=5, max.dx.C=1)
##D str(res)
## End(Not run)



### Name: create.X.grid.xD
### Title: Creates a Uniform Grid in a Tubular Neighborhood of a State
###   Space
### Aliases: create.X.grid.xD

### ** Examples

## Not run: 
##D # Let X be a low-dimensional model of a state space.
##D res <- create.X.grid.xD(X, gSf=5, gRf=5, p.exp=0.05, wC=0.1)
##D str(res)
## End(Not run)



### Name: critical.points.plot
### Title: Plot Critical Points
### Aliases: critical.points.plot

### ** Examples

# Create a function with known critical points
f <- function(x, y) x^2 - y^2  # Saddle at origin
gradient <- function(x, y) c(2*x, -2*y)
grid <- create.grid(30)
critical <- find.critical.points.continuous(gradient, grid)
critical.points.plot(critical)



### Name: cross.prod
### Title: Computes a cross product between two 3D vectors
### Aliases: cross.prod

### ** Examples

## Not run: 
##D x <- c(0,1,3)
##D y <- c(2,3,4)
##D z <- cross.prod(x, y)
## End(Not run)



### Name: cv.imputation
### Title: Cross-Validation Imputation on Graphs
### Aliases: cv.imputation

### ** Examples

## Not run: 
##D # Example with a small graph
##D graph <- list(c(2,3), c(1,3), c(1,2,4), c(3))
##D y <- c(1, 0, 1, 0)
##D test.set <- c(2, 4)
##D result <- cv.imputation(test.set, graph, y = y, y.binary = TRUE,
##D                         imputation.method = "LOCAL_MEAN_THRESHOLD")
##D print(result)
## End(Not run)




### Name: deg0.lowess.graph.smoothing
### Title: Iterative Degree 0 LOWESS Graph Smoothing
### Aliases: deg0.lowess.graph.smoothing

### ** Examples

## Not run: 
##D # Create a graph and data matrix
##D graph <- create.iknn.graph(X, k = 10, pruning.thld = 0.1)
##D 
##D # Apply degree 0 LOWESS graph smoothing - traditional approach
##D result1 <- deg0.lowess.graph.smoothing(
##D   adj.list = graph$adj_list,
##D   weight.list = graph$weight_list,
##D   X = X,
##D   max.iterations = 10,
##D   convergence.threshold = 1e-4,
##D   convergence.type = 1,  # MAX.ABSOLUTE.DIFF
##D   k = 10,
##D   pruning.thld = 0.1,
##D   n.bws = 10,
##D   n.folds = 5,
##D   verbose = TRUE
##D )
##D 
##D # Apply degree 0 LOWESS graph smoothing with boosting
##D result2 <- deg0.lowess.graph.smoothing(
##D   adj.list = graph$adj_list,
##D   weight.list = graph$weight_list,
##D   X = X,
##D   max.iterations = 10,
##D   convergence.threshold = 1e-4,
##D   convergence.type = 1,  # MAX.ABSOLUTE.DIFF
##D   k = 10,
##D   pruning.thld = 0.1,
##D   n.bws = 10,
##D   n.folds = 5,
##D   switch.to.residuals.after = 2,  # Switch to boosting after 2 iterations
##D   verbose = TRUE
##D )
##D 
##D # Access final smoothed data matrix
##D X.smoothed <- result2$smoothed.X[[length(result2$smoothed.X)]]
##D 
##D # Plot convergence metrics for both approaches
##D plot(result1$convergence.metrics, type = "b", col = "blue",
##D      xlab = "Iteration", ylab = "Convergence Metric")
##D lines(result2$convergence.metrics, type = "b", col = "red")
##D legend("topright", legend = c("Traditional", "Boosting"),
##D        col = c("blue", "red"), lty = 1)
## End(Not run)




### Name: delta.indices
### Title: Calculate Higher-Order Functional Association Indices
### Aliases: delta.indices

### ** Examples

## Not run: 
##D # Generate example conditional mean curve
##D x <- seq(0, 1, length.out = 100)
##D Eyg <- sin(4*pi*x) + 0.5*x
##D 
##D # Calculate indices up to order 5
##D indices <- delta.indices(Eyg, k = 5)
##D print(indices$Delta)
##D print(indices$delta)
##D 
##D # Plot delta values by order
##D plot(1:5, indices$delta, type = "b",
##D      xlab = "Order", ylab = "delta",
##D      main = "Functional Association by Order")
## End(Not run)




### Name: derivative.second.order.method
### Title: Second-order accurate method of derivative of a function
###   estimate
### Aliases: derivative.second.order.method

### ** Examples

## Not run: 
##D # Example with a quadratic function
##D x <- seq(0, 2*pi, length.out = 100)
##D y <- sin(x)
##D dx <- x[2] - x[1]
##D dy <- derivative.second.order.method(y, dx)
##D 
##D # Compare with analytical derivative
##D dy_true <- cos(x)
##D plot(x, dy, type = "l", col = "blue", main = "Numerical vs Analytical Derivative")
##D lines(x, dy_true, col = "red", lty = 2)
##D legend("topright", c("Numerical", "Analytical"), col = c("blue", "red"), lty = c(1, 2))
## End(Not run)



### Name: detect.local.extrema
### Title: Detect Local Extrema in a Graph
### Aliases: detect.local.extrema

### ** Examples

# Create a simple chain graph
adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
y <- c(1, 3, 2, 5, 1)  # Function values with peaks at vertices 2 and 4

# Detect maxima
maxima <- detect.local.extrema(adj.list, weight.list, y,
                               max.radius = 2,
                               min.neighborhood.size = 2)
print(maxima$vertices)  # Should identify vertices 2 and 4

# Detect minima
minima <- detect.local.extrema(adj.list, weight.list, y,
                               max.radius = 2,
                               min.neighborhood.size = 2,
                               detect.maxima = FALSE)
print(minima$vertices)  # Should identify vertices 1, 3, and 5




### Name: disk.to.sphere
### Title: Maps points from the interior of a unit disk onto the unit
###   sphere of the same dimension as the disk
### Aliases: disk.to.sphere

### ** Examples

# Example 1: Map 2D disk points to 3D sphere
X <- matrix(c(0.5, 0.3, -0.2, 0.4, 0, 0), nrow = 3, ncol = 2, byrow = TRUE)
sphere_points <- disk.to.sphere(X)
# Verify all points are on unit sphere
apply(sphere_points, 1, function(x) sqrt(sum(x^2)))

# Example 2: Center point maps to north pole
origin <- matrix(c(0, 0), nrow = 1)
disk.to.sphere(origin)  # Returns (0, 0, 1)



### Name: dist.to.knn
### Title: Produces k-NN distance and index matrices associated with a
###   distance matrix
### Aliases: dist.to.knn

### ** Examples

# Create a simple distance matrix
d <- as.matrix(dist(matrix(rnorm(20), ncol=2)))
knn_result <- dist.to.knn(d, k=3)
# knn_result$nn.i contains indices of 3 nearest neighbors for each point
# knn_result$nn.d contains distances to those neighbors




### Name: dlaplace
### Title: Laplace Distribution Density Function
### Aliases: dlaplace

### ** Examples

x <- seq(-5, 5, by = 0.1)
y <- dlaplace(x, location = 0, scale = 1)
plot(x, y, type = "l", main = "Laplace Density")




### Name: draw.3d.line
### Title: Draw 3D Line Segment
### Aliases: draw.3d.line

### ** Examples

## Not run: 
##D plot3D.plain(matrix(0, ncol = 3))
##D draw.3d.line(c(1, 1, 1), length = 2, col = "red")
##D draw.3d.line(c(1, 0, 0), length = 1.5, col = "blue")
## End(Not run)




### Name: draw.axes
### Title: Draw 3D Axes
### Aliases: draw.axes

### ** Examples

## Not run: 
##D plot3D.plain(matrix(rnorm(30), ncol = 3))
##D draw.axes(delta = 0.5, axes.color = "black", half.axes = TRUE)
## End(Not run)




### Name: draw.dashed.line3d
### Title: Draw Dashed Line in 3D Space
### Aliases: draw.dashed.line3d

### ** Examples

## Not run: 
##D # Draw a dashed line from origin to (1, 1, 1)
##D plot3D.plain(matrix(c(0,0,0,1,1,1), ncol = 3, byrow = TRUE))
##D draw.dashed.line3d(0, 0, 0, 1, 1, 1, col = "red", lwd = 3)
## End(Not run)




### Name: E.geodesic.X
### Title: Estimate Mean Abundance Along a Geodesic
### Aliases: E.geodesic.X

### ** Examples

## Not run: 
##D # Create example data
##D set.seed(123)
##D X <- matrix(rnorm(1000), nrow=100, ncol=10)
##D rownames(X) <- paste0("sample", 1:100)
##D colnames(X) <- paste0("var", 1:10)
##D 
##D # Define strand
##D strand.ids <- paste0("sample", 1:80)
##D dist.along.strand.geodesic <- seq(0, 1, length.out=80)
##D names(dist.along.strand.geodesic) <- strand.ids
##D 
##D # Run analysis
##D result <- E.geodesic.X(X = X,
##D                        strand.ids = strand.ids,
##D                        dist.along.strand.geodesic = dist.along.strand.geodesic,
##D                        data.dir = tempdir(),
##D                        prev.thld = 0.5)
##D 
##D # Check results
##D str(result)
## End(Not run)




### Name: ecdf.cpp
### Title: Empirical Cumulative Distribution Function (ECDF)
### Aliases: ecdf.cpp

### ** Examples

x <- c(1.0, 2.0, 3.0, 4.0, 5.0)
result <- ecdf.cpp(x)
print(result)




### Name: edge.diff
### Title: Compute Edge Difference Between Two Graphs
### Aliases: edge.diff

### ** Examples

graph1 <- list(c(2,3), c(1,3), c(1,2), c(5), c(4))
graph2 <- list(c(2), c(1), c(), c(5), c(4))
diff <- edge.diff(graph1, graph2)
print(diff)




### Name: eigen.ulogit
### Title: Fit Univariate Logistic Regression Using Eigen
### Aliases: eigen.ulogit

### ** Examples

# Generate example data with non-linear relationship
set.seed(456)
x <- runif(150, -3, 3)
true_prob <- 1/(1 + exp(-(1 + 2*x - 0.5*x^2)))
y <- rbinom(150, 1, prob = true_prob)

# Fit linear model
fit_linear <- eigen.ulogit(x, y)

# Fit quadratic model
fit_quad <- eigen.ulogit(x, y, fit.quadratic = TRUE)

# Compare models using AIC
cat("Linear model AIC:", fit_linear$aic, "\n")
cat("Quadratic model AIC:", fit_quad$aic, "\n")

# Plot both fits
ord <- order(x)
plot(x, y, pch = 16, col = ifelse(y == 1, "blue", "red"),
     main = "Linear vs Quadratic Logistic Regression")
lines(x[ord], fit_linear$predictions[ord], lwd = 2, col = "green")
lines(x[ord], fit_quad$predictions[ord], lwd = 2, col = "purple")
legend("topleft", c("y = 1", "y = 0", "Linear", "Quadratic"),
       col = c("blue", "red", "green", "purple"),
       pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1))

# Example with weights and without errors for efficiency
w <- runif(150, 0.5, 2)
fit_fast <- eigen.ulogit(x, y, w = w, with.errors = FALSE)

# Access coefficients
cat("Quadratic model coefficients:\n")
cat("  Intercept:", fit_quad$beta[1], "\n")
cat("  Linear:", fit_quad$beta[2], "\n")
cat("  Quadratic:", fit_quad$beta[3], "\n")




### Name: elapsed.time
### Title: Format and Print Elapsed Time
### Aliases: elapsed.time

### ** Examples

## Not run: 
##D start <- proc.time()
##D Sys.sleep(2)  # Do some work
##D elapsed.time(start, "Processing complete")
## End(Not run)



### Name: energy.distance
### Title: Compute Energy Distance Between Two Datasets
### Aliases: energy.distance

### ** Examples

X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- energy.distance(X, Y)
print(result)




### Name: entropy.difference
### Title: Estimate Relative Entropy Using Entropy Difference
### Aliases: entropy.difference

### ** Examples

X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- entropy.difference(X, Y)
print(result)




### Name: estimate_density
### Title: Define a function to estimate the density rho_S(X)
### Aliases: estimate_density

### ** Examples

# Current placeholder usage
X <- matrix(rnorm(100), ncol=2)
point <- c(0, 0)
density <- estimate_density(point, X)  # Returns 1




### Name: estimate.geodesic.distances
### Title: Estimate Geodesic Distances Between Points
### Aliases: estimate.geodesic.distances

### ** Examples

# Generate sample data on a Swiss roll manifold
n <- 100
t <- seq(0, 4*pi, length.out = n)
points <- cbind(
  x = t * cos(t),
  y = 10 * runif(n),
  z = t * sin(t)
)

# Estimate geodesic distances using k-NN graph
geo_dist_knn <- estimate.geodesic.distances(points, k = 5)

# Estimate geodesic distances using MST
geo_dist_mst <- estimate.geodesic.distances(points, method = "mst")

# Compare with Euclidean distances
euc_dist <- as.matrix(dist(points))

# Geodesic distances are typically larger than Euclidean for manifold data
mean(geo_dist_knn > euc_dist, na.rm = TRUE)




### Name: estimate.optimal.bandwidth.from.extrema.elbow
### Title: Estimate Optimal Bandwidth Using Local Extrema Count Elbow
###   Method
### Aliases: estimate.optimal.bandwidth.from.extrema.elbow

### ** Examples

# Simulated extrema counts for increasing bandwidths
extrema.counts <- c(150, 142, 128, 95, 72, 58, 45, 38, 35, 33, 32, 31)

# Estimate optimal bandwidth
opt.idx <- estimate.optimal.bandwidth.from.extrema.elbow(
  extrema.counts,
  sd.multiplier = 1.5,
  plot.results = TRUE
)

print(paste("Optimal bandwidth index:", opt.idx))




### Name: euclidean.distance
### Title: Function to calculate the Euclidean distance between two points
### Aliases: euclidean.distance

### ** Examples

# Distance between two 2D points
point1 <- c(0, 0)
point2 <- c(3, 4)
euclidean.distance(point1, point2)  # Returns 5

# Distance between two 3D points
point1 <- c(1, 2, 3)
point2 <- c(4, 6, 8)
euclidean.distance(point1, point2)  # Returns sqrt(50)




### Name: evaluate.function.on.grid.as.vector
### Title: Evaluate Function on Grid as Vector
### Aliases: evaluate.function.on.grid.as.vector

### ** Examples

grid <- create.grid(50)
f <- function(x, y) x^2 + y^2
values <- evaluate.function.on.grid.as.vector(f, grid)



### Name: evaluate.function.on.grid
### Title: Evaluate Function on Grid
### Aliases: evaluate.function.on.grid

### ** Examples

## Not run: 
##D grid <- create.grid(50)
##D f <- function(x, y) sin(x) * cos(y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D image(grid$x, grid$y, f.grid)
## End(Not run)



### Name: ext.graph.diffusion.smoother
### Title: Extended Graph Diffusion Smoother
### Aliases: ext.graph.diffusion.smoother

### ** Examples

## Not run: 
##D # Create a simple chain graph
##D n <- 20
##D graph <- vector("list", n)
##D edge.lengths <- vector("list", n)
##D for(i in 1:(n-1)) {
##D   graph[[i]] <- c(graph[[i]], i+1)
##D   graph[[i+1]] <- c(graph[[i+1]], i)
##D   edge.lengths[[i]] <- c(edge.lengths[[i]], 1)
##D   edge.lengths[[i+1]] <- c(edge.lengths[[i+1]], 1)
##D }
##D 
##D # Create noisy signal
##D true.signal <- sin(seq(0, 2*pi, length.out = n))
##D y <- true.signal + rnorm(n, 0, 0.2)
##D 
##D # Apply diffusion smoothing
##D result <- ext.graph.diffusion.smoother(
##D   graph = graph,
##D   edge.lengths = edge.lengths,
##D   y = y,
##D   n.time.steps = 50,
##D   step.factor = 0.1,
##D   n.CVs = 5,
##D   verbose = TRUE
##D )
##D 
##D # Plot results
##D plot(y, type = "b", col = "red", main = "Graph Diffusion Smoothing")
##D lines(result$y.optimal, type = "b", col = "blue")
##D legend("topright", c("Original", "Smoothed"), col = c("red", "blue"), lty = 1)
## End(Not run)




### Name: extract.derivatives
### Title: Extract Derivative Information from assoc1 Object
### Aliases: extract.derivatives

### ** Examples

## Not run: 
##D result <- fassoc1.test(x, y)
##D deriv.est <- extract.derivatives(result)
##D deriv.ci <- extract.derivatives(result, type = "credible.interval")
## End(Not run)



### Name: extract.extrema.subgraphs
### Title: Extract Local Subgraphs Around Extrema
### Aliases: extract.extrema.subgraphs

### ** Examples

## Not run: 
##D # Create a simple graph
##D adj_list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
##D weight_list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
##D y <- c(1, 2, 3, 1)
##D predictions <- c(1.1, 2.1, 2.9, 1.1)
##D 
##D # Create extrema data frame
##D extrema_df <- data.frame(
##D   vertex = 3,
##D   hop_idx = 1,
##D   is_max = TRUE,
##D   label = "Max1",
##D   fn_value = 3,
##D   stringsAsFactors = FALSE
##D )
##D 
##D # Extract subgraph
##D result <- extract.extrema.subgraphs(
##D   adj_list, weight_list, y, predictions, extrema_df, hop.offset = 1
##D )
## End(Not run)



### Name: extract.xy
### Title: Extracts Numerical x and y Values from a Character Vector
### Aliases: extract.xy

### ** Examples

## Not run: 
##D s <- c("1,2", "3,4", "5,6")
##D result <- extract_xy(s)
##D x <- result$x
##D y <- result$y
## End(Not run)



### Name: fassoc.test
### Title: Functional Association Test
### Aliases: fassoc.test functional.association.test

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 200
##D x <- runif(n)
##D y <- sin(2*pi*x) + rnorm(n, sd = 0.3)
##D 
##D # Test for zero-order functional association
##D result0 <- fassoc.test(x, y, order = 0, n.cores = 2)
##D 
##D # Test for first-order functional association
##D result1 <- fassoc.test(x, y, order = 1, n.cores = 2)
##D 
##D # Compare p-values
##D print(result0$p.value)
##D print(result1$p.value)
## End(Not run)




### Name: fassoc0.test
### Title: Tests for Existence of Non-trivial Functional Association
###   Between Two Variables
### Aliases: fassoc0.test zofam.test

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 200
##D x <- runif(n)
##D y <- sin(2*pi*x) + rnorm(n, sd = 0.3)
##D 
##D # Test for functional association
##D result <- fassoc0.test(x, y, n.cores = 2, n.perms = 1000)
##D 
##D # Custom plot
##D plot(result, plot = "Exy")
## End(Not run)




### Name: fassoc1.test
### Title: Tests for First-Order Functional Association Between Two
###   Variables
### Aliases: fassoc1.test fofam.test

### ** Examples

## Not run: 
##D # Generate example data with nonlinear relationship
##D set.seed(123)
##D n <- 200
##D x <- runif(n)
##D y <- sin(4*pi*x) + rnorm(n, sd = 0.2)
##D 
##D # Test for first-order functional association
##D result <- fassoc1.test(x, y, n.cores = 2, n.perms = 500)
##D 
##D # Plot derivative of conditional mean
##D plot(result, plot = "dExy")
## End(Not run)




### Name: find.box.containing.x
### Title: Finds the box within a box list covering some data set that
###   contains a given point
### Aliases: find.box.containing.x

### ** Examples

## Not run: 
##D boxes <- create.ED.boxes(0.5, c(0,0), c(1,1))
##D x <- c(0.2, 0.3)
##D containing_box <- find.box.containing.x(boxes, x)
## End(Not run)



### Name: find.boxes.containing.x
### Title: Finds boxes within a box list covering some data set that (after
###   expansion by a factor p.exp) contain a given point
### Aliases: find.boxes.containing.x

### ** Examples

## Not run: 
##D boxes <- create.ED.boxes(0.5, c(0,0), c(1,1))
##D x <- c(0.2, 0.3)
##D containing_box <- find.box.containing.x(boxes, x)
## End(Not run)



### Name: find.critical.points
### Title: Find Critical Points in Function Grid
### Aliases: find.critical.points

### ** Examples

## Not run: 
##D grid <- create.grid(50)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D critical <- find.critical.points(f.grid)
##D print(critical)
## End(Not run)



### Name: find.gflow.basins
### Title: Find Gradient-Flow Basins on a Weighted Graph
### Aliases: find.gflow.basins

### ** Examples

# Create a simple graph
adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
weight.list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
y <- c(1, 0, 2, 1.5)

## No test: 
# Find basins (requires compiled C++ backend)
result <- find.gflow.basins(adj.list, weight.list, y,
                            min.basin.size = 1,
                            min.path.size = 2,
                            q.edge.thld = 0.5)
summary(result)
## End(No test)




### Name: find.inflection.pts
### Title: Finds inflection points of a vector
### Aliases: find.inflection.pts

### ** Examples

## Not run: 
##D # Example 1: Find inflection points in a sigmoid-like curve
##D x <- seq(-5, 5, length.out = 100)
##D y <- 1 / (1 + exp(-x)) + 0.1 * sin(2*x)  # Sigmoid with oscillation
##D 
##D result <- find.inflection.pts(y, x, method = "magelo")
##D 
##D # Plot the results
##D plot(x, y, type = "l", main = "Inflection Points")
##D points(x[result$infl.pts], result$infl.pts.vals, col = "red", pch = 19)
##D 
##D # Example 2: Using indices as x-coordinates
##D y <- c(1, 2, 5, 10, 15, 18, 19, 19.5, 19.8, 20)  # Growth curve
##D result <- find.inflection.pts(y, method = "mabilo")
##D print(result$infl.pts)
## End(Not run)



### Name: find.local.extrema
### Title: Detect Local Extrema on a Graph with Labeled Basins
### Aliases: find.local.extrema

### ** Examples

# Create a simple graph
adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
weight.list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
y <- c(1, 0, 2, 1.5)

## No test: 
# Find local extrema (requires compiled C++ backend)
result <- find.local.extrema(adj.list, weight.list, y, min.basin.size = 1)
print(result$local_extrema)
## End(No test)




### Name: find.local.minima
### Title: Find Local Minima in Distance Sequence
### Aliases: find.local.minima

### ** Examples

# Create sample data with clear minima
k <- 1:20
dist <- sin(k / 3) + 0.1 * rnorm(20)

# Find local minima
minima <- find.local.minima(dist, k)
print(minima)

# Plot with minima marked
plot(k, dist, type = "b")
points(minima$k, minima$distance, col = "red", pch = 19, cex = 1.5)

# Example with NA handling
dist_na <- dist
dist_na[c(5, 15)] <- NA
minima_na <- find.local.minima(dist_na, k, na.rm = TRUE)




### Name: find.longest.true.stretch
### Title: Finds the Longest Stretch of Consecutive TRUE Values in a
###   Logical Vector
### Aliases: find.longest.true.stretch

### ** Examples

l <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
result <- find.longest.true.stretch(l)
print(result)
# $start_index
# [1] 4
# $length
# [1] 3

## Two longest TRUE stretches
l <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
find.longest.true.stretch(l)
# $start.index
# [1] 1
# $length
# [1] 3

## No TRUE values
l <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
find.longest.true.stretch(l)
# $start.index
# [1] NA
# $length
# [1] 0




### Name: find.optimal.k
### Title: Find Optimal k Parameter Using Edge Persistence Analysis
### Aliases: find.optimal.k

### ** Examples

## Not run: 
##D # Assuming we have birth-death data from create.iknn.graphs()
##D result <- create.iknn.graphs(X, kmin = 3, kmax = 10)
##D 
##D # Find optimal k
##D stability <- find.optimal.k(result$birth_death_matrix, kmin = 3, kmax = 10)
##D 
##D # Plot stability scores
##D plot(stability$k.values, stability$stability.scores, type = "l",
##D      xlab = "k", ylab = "Stability Score")
##D abline(v = stability$opt.k, col = "red", lty = 2)
## End(Not run)




### Name: find.points.within.box
### Title: Find the points of a set X within a box
### Aliases: find.points.within.box

### ** Examples

## Not run: 
##D X <- matrix(runif(100), ncol = 2)
##D box <- list()
##D box$L <- c(0.2, 0.2)
##D box$R <- c(0.5, 0.5)
##D points_within_box <- find.points.within.box(X, box, eps = 0)
## End(Not run)




### Name: find.shortest.paths.within.radius
### Title: Find All Shortest Paths Within a Radius
### Aliases: find.shortest.paths.within.radius

### ** Examples

# Create a simple graph with 5 vertices
adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 5), c(2, 5), c(3, 4))
weight.list <- list(c(1, 2), c(1, 1, 3), c(2, 1, 2), c(3, 1), c(2, 1))

# Find all shortest paths within radius 3 from vertex 1
(res <- find.shortest.paths.within.radius(adj.list, weight.list, 1, 3))




### Name: fit.pwlm
### Title: Fit a Piecewise Linear Model
### Aliases: fit.pwlm

### ** Examples

# Generate sample data
x <- 1:20
y <- x + 2*x^2 + rnorm(20, 0, 10)

# Single breakpoint
fit1 <- fit.pwlm(x, y)

# Two breakpoints
fit2 <- fit.pwlm(x, y, n_breakpoints = 2)



### Name: function.contours.plot
### Title: Plot Function Contours
### Aliases: function.contours.plot

### ** Examples

## Not run: 
##D grid <- create.grid(50)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D function.contours.plot(grid, f.grid)
## End(Not run)



### Name: gaussian.second.derivative
### Title: Calculate Second Derivative of Gaussian Function
### Aliases: gaussian.second.derivative

### ** Examples

# Example 1: Basic usage
x <- seq(-5, 5, length.out = 100)
result <- gaussian.second.derivative(x, mu = 0, sd = 1)

# Example 2: Verify inflection points at mu +/- sd
mu <- 2
sd <- 1.5
x_inflection <- c(mu - sd, mu + sd)
gaussian.second.derivative(x_inflection, mu, sd)  # Should be ~0

# Example 3: Visualize the Gaussian and its second derivative
x <- seq(-4, 4, length.out = 200)
f <- exp(-(x^2)/2)  # Original function (mu=0, sd=1)
f_double_prime <- gaussian.second.derivative(x, mu = 0, sd = 1)

par(mfrow = c(2, 1))
plot(x, f, type = "l", main = "Gaussian Function", ylab = "f(x)")
abline(v = c(-1, 1), col = "red", lty = 2)  # Inflection points
plot(x, f_double_prime, type = "l", main = "Second Derivative", ylab = "f''(x)")
abline(h = 0, col = "gray")
abline(v = c(-1, 1), col = "red", lty = 2)  # Zero crossings



### Name: gdensity
### Title: Kernel Density Estimation on a Uniform Grid
### Aliases: gdensity

### ** Examples

data <- rnorm(100)
result <- gdensity(data, grid.size = 200, poffset = 0.1, bw = 0, kernel = 5)
plot(result$x, result$y, type = "l")




### Name: generate.1d.circular.gaussian.mixture
### Title: Generate a 1D Circular Gaussian Mixture
### Aliases: generate.1d.circular.gaussian.mixture

### ** Examples

# Generate a simple circular Gaussian mixture
result <- generate.1d.circular.gaussian.mixture(
  n.points = 200,
  x.knots = c(0, pi, 1.5*pi),
  y.knots = c(5, 8, 2.5),
  sd.knot = 0.3
)

# Plot the result
plot(result$x, result$y, type = "l", xlab = "Angle (radians)", ylab = "Value")




### Name: generate.1d.gaussian.mixture
### Title: Generate a 1D Gaussian Mixture with Direct Knot Specification
### Aliases: generate.1d.gaussian.mixture

### ** Examples

# Generate mixture with default parameters
result1 <- generate.1d.gaussian.mixture()
plot(result1$x, result1$y, type = 'l', main = "Default 1D Gaussian Mixture")

# Generate mixture with custom parameters
result2 <- generate.1d.gaussian.mixture(
  n.points = 200,
  x.knots = c(-10, 0, 5, 15),
  y.knots = c(2, 5, 8, 3),
  sd.knot = 0.8,
  x.offset = 5
)

# Generate mixture with noise
result3 <- generate.1d.gaussian.mixture(
  x.knots = c(-2, 2),
  y.knots = c(1, -1),
  add.noise = TRUE,
  noise.fraction = 0.15
)




### Name: generate.binary.sample.from.gaussian.mixture
### Title: Generate Binary Sample from Gaussian Mixture
### Aliases: generate.binary.sample.from.gaussian.mixture

### ** Examples

# First generate a Gaussian mixture
gm <- generate.1d.gaussian.mixture(
  x.knots = c(-2, 0, 2),
  y.knots = c(-1, 2, -1.5)
)

# Generate binary sample with default settings
binary.data <- generate.binary.sample.from.gaussian.mixture(gm)

# Plot the results
plot(gm$x, gm$y, type = 'l', col = 'blue',
     main = "Gaussian Mixture and Binary Sample")
points(binary.data$x, binary.data$y.binary, pch = 19,
       col = ifelse(binary.data$y.binary == 1, "green", "red"))
lines(binary.data$x, binary.data$y.prob, col = 'orange', lwd = 2)

# Generate regular grid sample with wider logit range
binary.regular <- generate.binary.sample.from.gaussian.mixture(
  gm,
  n.sample.points = 100,
  sample.method = "regular",
  logit.range = c(-5, 5)  # Steeper probability transitions
)




### Name: generate.circle.data
### Title: Generate Points on a Circle
### Aliases: generate.circle.data

### ** Examples

# Generate 100 equally spaced points on a circle with radius 2
df <- generate.circle.data(100, radius = 2)

# Generate 50 random points with Laplace noise
df_noisy <- generate.circle.data(50, noise = 0.1, type = "random", noise.type = "laplace")




### Name: generate.circle.graph
### Title: Generate a Circle Graph
### Aliases: generate.circle.graph

### ** Examples

# Generate a circle graph with 5 vertices and uniform angles
graph <- generate.circle.graph(5, type = "uniform")

# Generate a circle graph with 10 vertices and random angles
graph <- generate.circle.graph(10, type = "random", seed = 123)




### Name: generate.graph.gaussian.mixture
### Title: Generate a Gaussian Mixture Function on a Graph
### Aliases: generate.graph.gaussian.mixture

### ** Examples

## Not run: 
##D # Create a grid graph
##D grid <- create.grid.graph(10, 10)
##D 
##D # Generate a mixture of two Gaussians
##D centers <- c(1, 100)  # Corner and center vertices
##D amplitudes <- c(1.0, 0.7)
##D sigmas <- c(2.0, 3.0)
##D 
##D y <- generate.graph.gaussian.mixture(
##D   grid$adj.list,
##D   grid$weight.list,
##D   centers,
##D   amplitudes,
##D   sigmas
##D )
## End(Not run)




### Name: generate.noisy.circle.embedding
### Title: Generate Noisy Circle Data in Higher Dimensions
### Aliases: generate.noisy.circle.embedding

### ** Examples

# Generate 100 points on a noisy circle in 3D space
data <- generate.noisy.circle.embedding(100, dim = 3, radius = 1, noise = 0.1)




### Name: generate.partition
### Title: Generate Random Partition of an Integer
### Aliases: generate.partition

### ** Examples

generate.partition(20, 4, 3)  # Generates a partition of 20 into 4 parts, each >= 3
generate.partition(10, 3, 2)  # Generates a partition of 10 into 3 parts, each >= 2




### Name: generate.star.dataset
### Title: Generate Synthetic Dataset on a Star-Shaped Geometric Space
### Aliases: generate.star.dataset

### ** Examples

# Generate a 2D star graph with 5 arms and exponential function centered at origin
result <- generate.star.dataset(
  n.points = 20,
  n.arms = 5,
  min.n.pts.within.arm = 3,
  min.arm.length = 1,
  max.arm.length = 3,
  dim = 2,
  fn = "exp",
  fn.center = c(0, 0),
  fn.scale = 1,
  noise = "norm",
  noise.sd = 0.1
)




### Name: generate.synthetic.function.higher.dim
### Title: Generate a Synthetic Smooth Function in Higher Dimensions
### Aliases: generate.synthetic.function.higher.dim

### ** Examples

## Not run: 
##D # Example usage
##D X <- matrix(runif(20), ncol = 2) # Random 2D points
##D synthetic_fn <- generate_synthetic_function_higher_dim(X)
##D point <- c(0.5, 0.5) # A point in 2D space
##D value <- synthetic_fn(point) # Evaluate the function at the given point
## End(Not run)



### Name: geodesic.knn
### Title: Estimates geodesic (shortest path) nearest neighbors
### Aliases: geodesic.knn

### ** Examples

## Not run: 
##D # Generate sample data on a spiral
##D t <- seq(0, 4*pi, length.out = 100)
##D X <- cbind(t*cos(t), t*sin(t))
##D 
##D # Find 5 geodesic nearest neighbors
##D gnn <- geodesic.knn(X, k = 5, K = 10)
##D 
##D # Compare with Euclidean nearest neighbors
##D enn <- get.knn(X, k = 5)
## End(Not run)




### Name: geodesic.knnx
### Title: Estimates geodesic (shortest path) nearest neighbors of X.grid
###   points in X That is for each point of X.grid the k-NN's within X are
###   returned
### Aliases: geodesic.knnx

### ** Examples

## Not run: 
##D # Create sample data
##D X <- matrix(rnorm(200), ncol=2)
##D 
##D # Create a regular grid
##D x_seq <- seq(min(X[,1]), max(X[,1]), length.out=10)
##D y_seq <- seq(min(X[,2]), max(X[,2]), length.out=10)
##D X.grid <- as.matrix(expand.grid(x_seq, y_seq))
##D 
##D # Find 3 nearest neighbors in X for each grid point
##D result <- geodesic.knnx(X, X.grid, k=3)
##D 
##D # Access the nearest neighbor indices and distances
##D nn_indices <- result$nn.index
##D nn_distances <- result$nn.dist
## End(Not run)




### Name: get.1d.binary.models.MABs
### Title: Estimate Mean Absolute Bias for Multiple Non-linear Regression
###   Models with Binary Outcome
### Aliases: get.1d.binary.models.MABs

### ** Examples

## Not run: 
##D # Generate synthetic data
##D set.seed(123)
##D n <- 200
##D xt <- seq(0, 10, length.out = n)
##D df <- data.frame(
##D   linear = xt,
##D   quadratic = xt^2 / 10,
##D   sine = sin(xt) * 5,
##D   exponential = exp(xt/5) - 1
##D )
##D 
##D # Compare models for binary outcomes
##D results <- get.1d.binary.models.MABs(xt, df, n.subsamples = 100, n.cores = 1)
##D print(results$MAB.df)
## End(Not run)




### Name: get.1d.models.MABs
### Title: Estimate Mean Absolute Bias for Multiple Non-linear Regression
###   Models
### Aliases: get.1d.models.MABs

### ** Examples

## Not run: 
##D # Generate synthetic data
##D set.seed(123)
##D n <- 200
##D xt <- seq(0, 10, length.out = n)
##D df <- data.frame(
##D   linear = xt,
##D   quadratic = xt^2 / 10,
##D   sine = sin(xt),
##D   exponential = exp(xt/5) - 1
##D )
##D 
##D # Compare models
##D results <- get.1d.models.MABs(xt, df, sd = 0.5, n.subsamples = 100, n.cores = 1)
##D print(results$MAB.df)
## End(Not run)




### Name: get.edge.weights
### Title: Get Unique Edge Weights from a Weighted Graph in Parallel
### Aliases: get.edge.weights

### ** Examples

## Not run: 
##D # Create a simple undirected weighted graph
##D adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
##D weight.list <- list(c(5, 10), c(5, 7), c(10, 7))
##D 
##D # Get weights of all edges using 2 cores
##D edge_weights <- get.edge.weights(adj.list, weight.list, n.cores = safe_cores(1))
## End(Not run)




### Name: get.gam.spline.MAB
### Title: Estimate Mean Absolute Bias (MAB) of GAM Spline Model
### Aliases: get.gam.spline.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- runif(100, 0, 10)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- seq(0, 10, length.out = 50)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Run the function
##D result <- get.gam.spline.MAB(x, y, xt, yt, folds)
##D print(result$MAB)
##D 
##D # Example with binary data
##D y.binary <- rbinom(100, 1, plogis(sin(x)))
##D yt.binary <- plogis(sin(xt))
##D result.binary <- get.gam.spline.MAB(x, y.binary, xt, yt.binary, folds, y.binary = TRUE)
##D print(result.binary$MAB)
## End(Not run)




### Name: get.gaussian.mixture.2d
### Title: Generate 2D Gaussian Mixture
### Aliases: get.gaussian.mixture.2d

### ** Examples

## Not run: 
##D # Generate random 2D mixture
##D gm1 <- get.gaussian.mixture.2d(n.components = 3)
##D plot(gm1)
##D 
##D # Generate grid arrangement with custom amplitudes
##D gm2 <- get.gaussian.mixture.2d(
##D   n.components = 4,
##D   centers.strategy = "grid",
##D   amplitudes.strategy = "custom",
##D   custom.amplitudes = c(1, -0.5, 0.8, -0.3)
##D )
##D 
##D # Custom centers with anisotropic Gaussians
##D gm3 <- get.gaussian.mixture.2d(
##D   n.components = 2,
##D   centers.strategy = "custom",
##D   custom.centers = matrix(c(0.3, 0.3, 0.7, 0.7), ncol = 2, byrow = TRUE),
##D   correlation = 0.5
##D )
## End(Not run)



### Name: get.gaussian.mixture
### Title: Generate Synthetic Data from a Mixture of Gaussians
### Aliases: get.gaussian.mixture

### ** Examples

# Basic usage with default parameters
data <- get.gaussian.mixture()

# Custom mixture with specific parameters
data <- get.gaussian.mixture(
  n.points = 200,
  n.components = 3,
  x.knots.strategy = "jittered",
  sd.strategy = "stoch.adaptive",
  y.strategy = "mixed",
  noise.fraction = 0.1
)

# Using sigma-scaled offset positioning
data <- get.gaussian.mixture(
  n.components = 4,
  x.knots.strategy = "offset",
  sigma = 2,
  xmin.factor = 2,
  xmax.factor = 5,
  sd.strategy = "adaptive"
)

# Example with fixed standard deviations and custom amplitudes
data <- get.gaussian.mixture(
  n.components = 2,
  x.knots.strategy = "uniform",
  sd.strategy = "custom",
  custom.sd.knots = c(0.05, 0.15),
  y.strategy = "custom",
  custom.y.knots = c(1, -0.5)
)




### Name: get.gausspr.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Gaussian Process Regression
###   Model
### Aliases: get.gausspr.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D library(kernlab)
##D set.seed(123)
##D x <- matrix(rnorm(100), ncol = 1)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- matrix(rnorm(50), ncol = 1)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Run the function
##D result <- get.gausspr.MAB(x, y, xt, yt, folds)
##D print(result$MAB)
## End(Not run)




### Name: get.glmnet.poly.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Regularized Polynomial
###   Regression Model
### Aliases: get.glmnet.poly.MAB

### ** Examples

## Not run: 
##D # Example with continuous outcome
##D set.seed(123)
##D x <- seq(-2, 2, length.out = 100)
##D y <- x^2 + rnorm(100, sd = 0.5)
##D xt <- seq(-2, 2, length.out = 50)
##D yt <- xt^2
##D 
##D # Create cross-validation folds (not used internally)
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Fit model
##D result <- get.glmnet.poly.MAB(x, y, xt, yt, folds, alpha = 1, y.binary = FALSE)
##D print(result$MAB)
##D 
##D # Example with binary outcome
##D y.binary <- rbinom(100, 1, plogis(x^2))
##D yt.binary <- plogis(xt^2)
##D result.binary <- get.glmnet.poly.MAB(x, y.binary, xt, yt.binary, folds, y.binary = TRUE)
##D print(result.binary$MAB)
## End(Not run)




### Name: get.gpredictions.CrI
### Title: Generates Bayesian bootstrap credible intervals of gpredictions
### Aliases: get.gpredictions.CrI

### ** Examples

# TBD



### Name: get.gpredictionss
### Title: Estimates gpredictions Means Matrix Over a Uniform Grid for
###   Different bandwidths
### Aliases: get.gpredictionss

### ** Examples

## Not run: 
##D res <- get.gpredictionss(bws, nn.i, nn.d, nn.x, nn.y, y.binary, degree, min.K, xgrid)
##D str(res)
## End(Not run)



### Name: get.krr.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Kernel Ridge Regression
###   Model
### Aliases: get.krr.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- matrix(rnorm(100), ncol = 1)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- matrix(rnorm(50), ncol = 1)
##D yt <- sin(xt)
##D 
##D # Run the function
##D result <- get.krr.MAB(x, y, xt, yt)
##D print(result$MAB)
## End(Not run)




### Name: get.krr.MAB.xD
### Title: Estimate Mean Absolute Bias (MAB) of Kernel Ridge Regression for
###   Multi-dimensional Data
### Aliases: get.krr.MAB.xD

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 100
##D X <- matrix(rnorm(n * 3), ncol = 3)
##D y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
##D yt <- X[,1] + X[,2]^2  # True values
##D 
##D # Fit KRR model and compute MAB
##D result <- get.krr.MAB.xD(X, y, yt)
##D print(result$MAB)
##D 
##D # With custom parameters
##D result2 <- get.krr.MAB.xD(X, y, yt, kernel = 'polydot',
##D                          lambdas = 10^seq(-6, -2, length = 5))
##D print(result2$MAB)
## End(Not run)




### Name: get.locpoly.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Local Polynomial Regression
### Aliases: get.locpoly.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- runif(100, 0, 10)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- seq(0, 10, length.out = 50)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Define bandwidth values
##D bandwidths <- seq(0.1, 1, by = 0.1)
##D 
##D # Run the function
##D result <- get.locpoly.MAB(x, y, xt, yt, folds, bandwidths)
##D print(result$MAB)
## End(Not run)




### Name: get.loess.MAB
### Title: Estimate Mean Absolute Bias (MAB) of LOESS Model
### Aliases: get.loess.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- runif(100, 0, 10)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- seq(0, 10, length.out = 50)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Fit LOESS model and compute MAB
##D result <- get.loess.MAB(x, y, xt, yt, folds, deg = 2)
##D print(result$MAB)
## End(Not run)




### Name: get.magelo.MAB
### Title: Estimate Mean Absolute Bias (MAB) of MAGELO Model
### Aliases: get.magelo.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- runif(100, 0, 10)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- seq(0, 10, length.out = 50)
##D yt <- sin(xt)
##D 
##D # Fit magelo model and compute MAB
##D result <- get.magelo.MAB(x, y, xt, yt, deg = 2)
##D print(result$MAB)
## End(Not run)




### Name: get.np.MAB.xD
### Title: Estimate Mean Absolute Bias (MAB) of Nonparametric Kernel
###   Regression Model
### Aliases: get.np.MAB.xD

### ** Examples

## Not run: 
##D # Generate example data
##D library(np)
##D set.seed(123)
##D n <- 100
##D X <- matrix(rnorm(n * 2), ncol = 2)
##D colnames(X) <- c("x1", "x2")
##D y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
##D yt <- X[,1] + X[,2]^2  # True values
##D 
##D # Fit nonparametric model and compute MAB
##D result <- get.np.MAB.xD(X, y, yt)
##D print(result$MAB)
##D 
##D # Example with binary data
##D y.binary <- rbinom(n, 1, plogis(y))
##D yt.binary <- plogis(yt)
##D result.binary <- get.np.MAB.xD(X, y.binary, yt.binary, y.binary = TRUE)
##D print(result.binary$MAB)
## End(Not run)




### Name: get.path.data
### Title: Compute Graph Path Data with Kernel Weights
### Aliases: get.path.data

### ** Examples

## Not run: 
##D # Create a simple graph with 5 vertices
##D adj.list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2), c(3))
##D weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1), c(1))
##D y <- c(1.5, 2.0, 0.5, 1.0, 1.5)
##D 
##D # Find paths centered around vertex 2
##D paths <- get.path.data(adj.list, weight.list, y,
##D                       ref.vertex = 2, bandwidth = 2)
## End(Not run)



### Name: get.random.forest.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Random Forest Model
### Aliases: get.random.forest.MAB

### ** Examples

## Not run: 
##D # Regression example
##D set.seed(123)
##D x <- matrix(rnorm(200), ncol = 2)
##D y <- x[,1] + x[,2]^2 + rnorm(100, sd = 0.5)
##D xt <- matrix(rnorm(100), ncol = 2)
##D yt <- xt[,1] + xt[,2]^2
##D 
##D # Basic usage
##D result <- get.random.forest.MAB(x, y, xt, yt, ntree = 500)
##D print(result$MAB)
##D 
##D # With optimization
##D result.opt <- get.random.forest.MAB(x, y, xt, yt, optimize.ntree = TRUE)
##D print(result.opt$parameters$optimal.ntree)
##D 
##D # Classification example
##D y.class <- factor(ifelse(y > median(y), "high", "low"))
##D yt.prob <- pnorm(yt, mean = mean(yt), sd = sd(yt))
##D result.class <- get.random.forest.MAB(x, y.class, xt, yt.prob)
##D print(result.class$MAB)
## End(Not run)




### Name: get.random.forest.MAB.xD
### Title: Estimate Mean Absolute Bias (MAB) of Random Forest Model for
###   Multi-dimensional Data
### Aliases: get.random.forest.MAB.xD

### ** Examples

## Not run: 
##D # Regression example
##D set.seed(123)
##D n <- 100
##D X <- matrix(rnorm(n * 3), ncol = 3)
##D y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
##D yt <- X[,1] + X[,2]^2  # True values
##D 
##D result <- get.random.forest.MAB.xD(X, y, yt, ntree = 500)
##D print(result$MAB)
##D 
##D # Classification example
##D y.class <- factor(ifelse(y > median(y), "high", "low"))
##D yt.prob <- plogis(yt - median(yt))  # True probabilities
##D result.class <- get.random.forest.MAB.xD(X, y.class, yt.prob, ntree = 500)
##D print(result.class$MAB)
## End(Not run)




### Name: get.region.boundary
### Title: Get Boundary Vertices of a Region in a Graph
### Aliases: get.region.boundary

### ** Examples

## Not run: 
##D # Create a simple grid graph adjacency list (5x5 grid)
##D create_grid_adj_list <- function(n_rows, n_cols) {
##D   n <- n_rows * n_cols
##D   adj_list <- vector("list", n)
##D   for (i in 1:n_rows) {
##D     for (j in 1:n_cols) {
##D       v <- (i-1) * n_cols + j
##D       neighbors <- numeric(0)
##D 
##D       # Add neighbors (up, down, left, right)
##D       if (i > 1) neighbors <- c(neighbors, (i-2) * n_cols + j)  # up
##D       if (i < n_rows) neighbors <- c(neighbors, i * n_cols + j)  # down
##D       if (j > 1) neighbors <- c(neighbors, (i-1) * n_cols + (j-1))  # left
##D       if (j < n_cols) neighbors <- c(neighbors, (i-1) * n_cols + (j+1))  # right
##D 
##D       adj_list[[v]] <- neighbors
##D     }
##D   }
##D   return(adj_list)
##D }
##D 
##D # Create a 5x5 grid graph
##D grid_adj_list <- create_grid_adj_list(5, 5)
##D 
##D # Define a region (central 3x3 subgrid)
##D central_region <- c(7:9, 12:14, 17:19)
##D 
##D # Find boundary vertices of the central region
##D boundary <- get.region.boundary(grid_adj_list, central_region)
##D print(boundary)
##D # Expected output: c(7, 8, 9, 12, 14, 17, 18, 19)
##D # (all except the center vertex 13)
## End(Not run)




### Name: get.shortest.path
### Title: Get Shortest Path Between Two Vertices
### Aliases: get.shortest.path

### ** Examples

## Not run: 
##D # Create a path graph
##D pg <- create.path.graph(graph, edge.lengths, h = 3)
##D 
##D # Get shortest path from vertex 1 to vertex 3
##D path_info <- get.shortest.path(pg, 1, 3)
##D if (!is.null(path_info)) {
##D   cat("Path:", path_info$path, "\n")
##D   cat("Length:", path_info$length, "\n")
##D }
## End(Not run)




### Name: get.smooth.spline.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Smooth Spline Model
### Aliases: get.smooth.spline.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- runif(100, 0, 10)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- seq(0, 10, length.out = 50)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Run the function
##D result <- get.smooth.spline.MAB(x, y, xt, yt, folds)
##D print(result$MAB)
## End(Not run)




### Name: get.sphere.degree.props
### Title: Calculate Degree Distribution Properties for Random Points on a
###   Sphere
### Aliases: get.sphere.degree.props

### ** Examples

## Not run: 
##D # Calculate degree distribution properties for 1000 points on a circle
##D circle_props <- get.sphere.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 2)
##D 
##D # Calculate for points on a sphere
##D sphere_props <- get.sphere.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 3)
## End(Not run)



### Name: get.spline.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Spline Model
### Aliases: get.spline.MAB

### ** Examples

## Not run: 
##D # Example with continuous outcome
##D set.seed(123)
##D x <- seq(0, 2*pi, length.out = 100)
##D y <- sin(x) + rnorm(100, sd = 0.2)
##D xt <- seq(0, 2*pi, length.out = 50)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 5, list = TRUE, returnTrain = TRUE)
##D 
##D # Fit spline model
##D result <- get.spline.MAB(x, y, xt, yt, folds, y.binary = FALSE)
##D print(result$MAB)
##D 
##D # Example with binary outcome
##D y.binary <- rbinom(100, 1, plogis(sin(x)))
##D yt.binary <- plogis(sin(xt))
##D result.binary <- get.spline.MAB(x, y.binary, xt, yt.binary, folds, y.binary = TRUE)
##D print(result.binary$MAB)
## End(Not run)




### Name: get.spline.MAB.xD
### Title: Estimate Mean Absolute Bias (MAB) of GAM Spline Model for
###   Multi-dimensional Data
### Aliases: get.spline.MAB.xD

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 100
##D X <- matrix(rnorm(n * 2), ncol = 2)
##D colnames(X) <- c("x1", "x2")
##D y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
##D yt <- X[,1] + X[,2]^2  # True values
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Fit GAM model and compute MAB
##D result <- get.spline.MAB.xD(X, y, yt, folds)
##D print(result$MAB)
##D 
##D # Example with binary data
##D y.binary <- rbinom(n, 1, plogis(y))
##D yt.binary <- plogis(yt)
##D result.binary <- get.spline.MAB.xD(X, y.binary, yt.binary, folds, y.binary = TRUE)
##D print(result.binary$MAB)
## End(Not run)




### Name: get.subj.ids
### Title: Get subject-specific visit information and indices
### Aliases: get.subj.ids

### ** Examples

# Example usage:
# S.3d <- matrix(1:12, nrow=4, dimnames=list(c("s1_v1", "s1_v2", "s2_v1", "s2_v2")))
# subjID <- c("s1", "s1", "s2", "s2")
# visit <- c(1, 2, 1, 2)
# get.subj.ids("s1", S.3d, subjID, visit)



### Name: get.svm.MAB
### Title: Estimate Mean Absolute Bias (MAB) of Support Vector Machine
###   Model
### Aliases: get.svm.MAB

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- matrix(rnorm(100), ncol = 1)
##D y <- sin(x) + rnorm(100, sd = 0.1)
##D xt <- matrix(rnorm(50), ncol = 1)
##D yt <- sin(xt)
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Define hyperparameter grid
##D cost <- 10^seq(-1, 1, length.out = 3)
##D gamma <- 10^seq(-1, 1, length.out = 3)
##D 
##D # Run the function
##D result <- get.svm.MAB(x, y, xt, yt, folds, cost, gamma)
##D print(result$MAB)
## End(Not run)




### Name: get.svm.MAB.xD
### Title: Estimate Mean Absolute Bias (MAB) of SVM Model for
###   Multi-dimensional Data
### Aliases: get.svm.MAB.xD

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 100
##D X <- matrix(rnorm(n * 3), ncol = 3)
##D y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
##D yt <- X[,1] + X[,2]^2  # True values
##D 
##D # Create cross-validation folds
##D folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
##D 
##D # Define hyperparameter grid
##D cost <- 10^seq(-1, 1, length.out = 3)
##D gamma <- 10^seq(-1, 1, length.out = 3)
##D 
##D # Fit SVM model and compute MAB
##D result <- get.svm.MAB.xD(X, y, yt, folds, cost, gamma)
##D print(result$MAB)
## End(Not run)




### Name: get.TN.S1.degree.props
### Title: Generate Degree Distribution Properties in Tubular Neighborhood
###   of Unit Circle
### Aliases: get.TN.S1.degree.props

### ** Examples

res <- get.TN.S1.degree.props(
  n.pts = 100,
  n.sims = 10,
  k = 5,
  noise = 0.05,
  noise.type = "normal"
)
str(res)



### Name: get.torus.degree.props
### Title: Calculate Degree Distribution Properties for Random Points on a
###   Torus
### Aliases: get.torus.degree.props

### ** Examples

## Not run: 
##D # Calculate degree distribution properties for 1000 points on a circle
##D circle_props <- get.torus.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 1)
##D 
##D # Calculate for points on a 2-torus
##D torus_props <- get.torus.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 2)
## End(Not run)



### Name: ggaussian
### Title: Generalized Gaussian Function
### Aliases: ggaussian

### ** Examples

## Not run: 
##D x <- seq(-5, 5, length.out = 100)
##D plot(x, ggaussian(x, p = 2), type = "l")  # Standard Gaussian
##D plot(x, ggaussian(x, p = 1), type = "l")  # Laplace distribution (heavier tails)
##D plot(x, ggaussian(x, p = 3), type = "l")  # Lighter tails than Gaussian
## End(Not run)



### Name: gradient.trajectories.plot
### Title: Plot Gradient Trajectories
### Aliases: gradient.trajectories.plot

### ** Examples

## Not run: 
##D grid <- create.grid(50)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D gradient.trajectories.plot(grid, f.grid, spacing = 15)
## End(Not run)



### Name: gradient.trajectory.plot
### Title: Plot Gradient Trajectory
### Aliases: gradient.trajectory.plot

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) -((x-0.5)^2 + (y-0.5)^2)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D traj <- compute.gradient.trajectory(15, 15, f.grid)
##D gradient.trajectory.plot(grid, traj)
## End(Not run)



### Name: graph.adj.mat
### Title: Creates weights adjacency matrix given a matrix of edges and
###   position matrix of the corresponding vertices
### Aliases: graph.adj.mat

### ** Examples

# Create a simple triangle graph
X <- matrix(c(0,0, 1,0, 0,1), ncol=2, byrow=TRUE)
E <- matrix(c(1,2, 2,3, 3,1), ncol=2, byrow=TRUE)
A <- graph.adj.mat(X, E)
print(A)
# Should show distances: A[1,2]=A[2,1]=1, A[2,3]=A[3,2]=sqrt(2), A[3,1]=A[1,3]=1




### Name: graph.connected.components
### Title: Assign Vertices to Connected Components
### Aliases: graph.connected.components

### ** Examples

adj.list <- list(c(2,3), c(1), c(1), c(5), c(4))
components <- graph.connected.components(adj.list)
print(components)




### Name: graph.constrained.gradient.flow.trajectories
### Title: Compute Constrained Gradient Flow Trajectories
### Aliases: graph.constrained.gradient.flow.trajectories

### ** Examples

## Not run: 
##D trajectories <- graph.constrained.gradient.flow.trajectories(graph, core.graph, Ey)
##D local.extrema <- trajectories$lext
##D paths <- trajectories$trajectories
## End(Not run)




### Name: graph.deg0.lowess.cv.mat
### Title: Matrix Version of Graph-Based LOWESS with Cross-Validation
### Aliases: graph.deg0.lowess.cv.mat

### ** Examples

## Not run: 
##D # Create a simple graph with 3 response variables
##D adj.list <- list(c(2,3), c(1,3), c(1,2))
##D weight.list <- list(c(1,1), c(1,1), c(1,1))
##D Y <- list(c(1,2,3), c(4,5,6), c(7,8,9))
##D 
##D # Run the algorithm
##D result <- graph.deg0.lowess.cv.mat(
##D   adj.list, weight.list, Y,
##D   min.bw.factor = 0.1, max.bw.factor = 0.5,
##D   n.bws = 10, log.grid = TRUE, kernel.type = 2,
##D   dist.normalization.factor = 1.1, n.folds = 3,
##D   with.bw.predictions = FALSE, precision = 1e-6, verbose = TRUE
##D )
##D 
##D # Access results
##D result$predictions  # Smoothed values for each response variable
##D result$opt.bws      # Optimal bandwidths for each response variable
## End(Not run)




### Name: graph.deg0.lowess
### Title: Degree 0 LOWESS (Locally Weighted Average) on Graphs
### Aliases: graph.deg0.lowess

### ** Examples

## Not run: 
##D # Create a simple ring graph
##D n <- 100
##D set.seed(123)
##D 
##D # Create a ring graph
##D adj.list <- vector("list", n)
##D weight.list <- vector("list", n)
##D 
##D for (i in 1:n) {
##D   neighbors <- c(i-1, i+1)
##D   # Handle wrap-around for ring structure
##D   neighbors[neighbors == 0] <- n
##D   neighbors[neighbors == n+1] <- 1
##D 
##D   adj.list[[i]] <- neighbors
##D   weight.list[[i]] <- rep(1, length(neighbors))
##D }
##D 
##D # Generate response values with spatial pattern plus noise
##D y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
##D 
##D # Apply degree 0 LOWESS with different bandwidths
##D bw1 <- 5
##D bw2 <- 10
##D bw3 <- 20
##D 
##D result1 <- graph.deg0.lowess(adj.list, weight.list, y, bw1)
##D result2 <- graph.deg0.lowess(adj.list, weight.list, y, bw2)
##D result3 <- graph.deg0.lowess(adj.list, weight.list, y, bw3)
##D 
##D # Plot results
##D plot(y, type="p", col="gray", main="Graph Degree 0 LOWESS")
##D lines(result1, col="blue", lwd=2)
##D lines(result2, col="red", lwd=2)
##D lines(result3, col="green", lwd=2)
##D legend("topright", legend=c(paste0("BW=", bw1), 
##D                             paste0("BW=", bw2),
##D                             paste0("BW=", bw3)),
##D        col=c("blue", "red", "green"), lwd=2)
## End(Not run)




### Name: graph.diffusion.smoother
### Title: Graph Diffusion Smoother
### Aliases: graph.diffusion.smoother

### ** Examples

## Not run: 
##D # Create simple chain graph
##D adj.list <- list(2, c(1,3), c(2,4), c(3,5), 4)
##D weight.list <- lapply(adj.list, function(x) rep(1, length(x)))
##D y <- c(1, 0, 1, 0, 1) + rnorm(5, 0, 0.1)
##D 
##D result <- graph.diffusion.smoother(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   n.time.steps = 20,
##D   n.CVs = 3,
##D   verbose = TRUE
##D )
## End(Not run)




### Name: graph.edit.distance
### Title: Calculate Graph Edit Distance (Pure R Implementation)
### Aliases: graph.edit.distance

### ** Examples

# Create two simple graphs
graph1.adj <- list(c(2, 3), c(1, 3), c(1, 2))
graph1.wts <- list(c(1, 2), c(1, 3), c(2, 3))

graph2.adj <- list(c(2, 3), c(1, 3), c(1, 2))
graph2.wts <- list(c(1.5, 2), c(1.5, 2.5), c(2, 2.5))

# Calculate distance
dist <- graph.edit.distance(graph1.adj, graph1.wts,
                           graph2.adj, graph2.wts)
print(dist)  # Weight differences only




### Name: graph.embedding
### Title: Create a 2D or 3D Embedding of a Graph
### Aliases: graph.embedding

### ** Examples

# Unweighted graph
adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
embedding <- graph.embedding(adj.list, dim = 2, method = "fr")

# Weighted graph
weights.list <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
embedding_weighted <- graph.embedding(adj.list, weights.list, dim = 2, method = "kk")




### Name: graph.kernel.smoother
### Title: Graph Kernel Smoother
### Aliases: graph.kernel.smoother

### ** Examples

## Not run: 
##D # Tiny chain graph with unit weights
##D n <- 5L
##D adj <- vector("list", n)
##D for (i in seq_len(n)) {
##D   adj[[i]] <- unique(c(if (i > 1) i - 1L else integer(), if (i < n) i + 1L else integer()))
##D }
##D w <- lapply(adj, function(nei) rep(1, length(nei)))
##D y <- as.numeric(1:n)
##D 
##D fit <- graph.kernel.smoother(adj, w, y, bandwidth = 2L, with_details = FALSE)
##D str(fit)
## End(Not run)




### Name: graph.kmean.cv
### Title: Cross-validation for Kernel-weighted Graph Neighbor Mean
### Aliases: graph.kmean.cv

### ** Examples

# Triangle graph example
graph <- list(c(2, 3), c(1, 3), c(1, 2))
edge.lengths <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
y <- c(1.0, 2.0, 3.0)

# Perform cross-validation
cv.errors <- graph.kmean.cv(graph, edge.lengths, y,
                            kernel = "epanechnikov",
                            n.CVs = 5, n.CV.folds = 3)
print(cv.errors)




### Name: graph.kmean
### Title: Compute the Kernel-weighted Graph Neighbor Mean
### Aliases: graph.kmean

### ** Examples

## Not run: 
##D # Simple triangle graph
##D graph <- list(c(2, 3), c(1, 3), c(1, 2))
##D edge.lengths <- list(c(0.1, 0.2), c(0.1, 0.3), c(0.2, 0.3))
##D y <- c(1.0, 2.0, 3.0)
##D 
##D # Compute kernel-weighted means with Epanechnikov kernel
##D y.kwmean <- graph.kmean(graph, edge.lengths, y, kernel = 1)
##D print(y.kwmean)
## End(Not run)




### Name: graph.low.pass.filter
### Title: Compute Low-Pass Filter of a Function Over a Graph
### Aliases: graph.low.pass.filter

### ** Examples

# Create example eigenvectors (2 vertices, 2 eigenvectors)
evectors <- matrix(c(1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2)),
                  nrow = 2, ncol = 2)

# Example GFT coefficients
y.gft <- matrix(c(3, 1), nrow = 2, ncol = 1)

# Apply low-pass filter starting from the first eigenvector
filtered <- graph.low.pass.filter(1, evectors, y.gft)

# Apply low-pass filter using only the second eigenvector
filtered_high <- graph.low.pass.filter(2, evectors, y.gft)




### Name: graph.mad
### Title: Calculate Median Absolute Deviation for Graph Vertices
### Aliases: graph.mad

### ** Examples

# Create a simple chain graph with 3 vertices
g <- list(c(2), c(1,3), c(2))
values <- c(1, 10, 2)
graph.mad(g, values)




### Name: graph.spectral.embedding
### Title: Generate Spectral Embedding of a Graph
### Aliases: graph.spectral.embedding

### ** Examples

# Create example eigenvectors for 5 vertices
set.seed(123)
evectors <- matrix(rnorm(5 * 5), nrow = 5, ncol = 5)

# Apply Gram-Schmidt to ensure orthogonality
evectors <- qr.Q(qr(evectors))

# Generate 2D embedding without eigenvalue scaling
embedding_2d <- graph.spectral.embedding(evectors, dim = 2)

# Generate 3D embedding with eigenvalue scaling
evalues <- c(5, 3, 2, 0.5, 0)  # Example eigenvalues in descending order
embedding_3d_scaled <- graph.spectral.embedding(evectors, dim = 3, evalues)




### Name: graph.spectral.filter
### Title: Spectral Filtering for Graph Signals
### Aliases: graph.spectral.filter

### ** Examples

## Not run: 
##D # Create a simple chain graph with 30 vertices (for CRAN check timing)
##D n <- 30
##D adj_list <- vector("list", n)
##D weight_list <- vector("list", n)
##D 
##D # Build adjacency structure for chain
##D for (i in seq_len(n)) {
##D   adj_list[[i]] <- integer(0)
##D   weight_list[[i]] <- numeric(0)
##D 
##D   if (i > 1) {
##D     adj_list[[i]] <- c(adj_list[[i]], i - 1L)
##D     weight_list[[i]] <- c(weight_list[[i]], 1.0)
##D   }
##D   if (i < n) {
##D     adj_list[[i]] <- c(adj_list[[i]], i + 1L)
##D     weight_list[[i]] <- c(weight_list[[i]], 1.0)
##D   }
##D }
##D 
##D # Create a noisy signal
##D set.seed(123)
##D x <- seq(0, 1, length.out = n)
##D y <- sin(2 * pi * x) + rnorm(n, 0, 0.2)
##D 
##D # Standard Laplacian with heat kernel filter
##D result1 <- graph.spectral.filter(
##D   adj.list = adj_list,
##D   weight.list = weight_list,
##D   y = y,
##D   laplacian.type = 0L,    # STANDARD
##D   filter.type = 0L,       # HEAT
##D   laplacian.power = 1L,
##D   n.candidates = 50       # Reduced for faster execution
##D )
##D 
##D # Cubic spline-like smoothing
##D result2 <- graph.spectral.filter(
##D   adj.list = adj_list,
##D   weight.list = weight_list,
##D   y = y,
##D   laplacian.type = 0L,    # STANDARD
##D   filter.type = 3L,       # CUBIC_SPLINE
##D   laplacian.power = 2L,
##D   n.candidates = 50
##D )
##D 
##D # Compare results
##D plot(x, y, pch = 16, col = "gray", main = "Graph Spectral Filtering",
##D      xlab = "Position", ylab = "Value")
##D lines(x, result1$predictions, col = "blue", lwd = 2)
##D lines(x, result2$predictions, col = "red", lwd = 2)
##D legend("topright",
##D        legend = c("Noisy data", "Heat kernel", "Cubic spline"),
##D        pch = c(16, NA, NA),
##D        lty = c(NA, 1, 1),
##D        col = c("gray", "blue", "red"),
##D        lwd = c(NA, 2, 2))
##D 
##D # Print summary
##D summary(result1)
## End(Not run)




### Name: graph.spectral.lowess.mat
### Title: Matrix version of Spectral Graph Local Polynomial Regression
### Aliases: graph.spectral.lowess.mat

### ** Examples

## Not run: 
##D # Create sample graph
##D n <- 100
##D g <- make.example.graph(n)
##D 
##D # Create multiple response variables (e.g., 3 columns)
##D Y <- matrix(rnorm(n*3), ncol=3)
##D 
##D # Apply spectral lowess to all columns at once
##D result <- graph.spectral.lowess.mat(g$adj.list, g$weight.list, Y)
##D 
##D # Plot results for first response variable
##D plot(Y[,1], result$predictions[,1], main="Original vs. Smoothed (Var 1)")
##D abline(0, 1, col="red", lty=2)
## End(Not run)




### Name: graph.spectral.lowess
### Title: Local Regression on Graphs Using Spectral Embedding
### Aliases: graph.spectral.lowess

### ** Examples

## Not run: 
##D # Create a simple graph with 100 vertices
##D n <- 100
##D set.seed(123)
##D 
##D # Create a ring graph
##D adj.list <- vector("list", n)
##D weight.list <- vector("list", n)
##D 
##D for (i in 1:n) {
##D   neighbors <- c(i-1, i+1)
##D   # Handle wrap-around for ring structure
##D   neighbors[neighbors == 0] <- n
##D   neighbors[neighbors == n+1] <- 1
##D 
##D   adj.list[[i]] <- neighbors
##D   weight.list[[i]] <- rep(1, length(neighbors))
##D }
##D 
##D # Generate response values with spatial pattern
##D y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
##D 
##D # Apply spectral LOWESS
##D result <- graph.spectral.lowess(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   n.evectors = 5,
##D   verbose = TRUE
##D )
##D 
##D # Plot results
##D plot(y, type="l", col="gray", main="Graph Spectral LOWESS")
##D lines(result$predictions, col="red", lwd=2)
## End(Not run)




### Name: graph.spectral.ma.lowess
### Title: Local Regression on Graphs Using Spectral Embedding
### Aliases: graph.spectral.ma.lowess

### ** Examples

## Not run: 
##D # Create a simple graph with 100 vertices
##D n <- 100
##D set.seed(123)
##D 
##D # Create a ring graph
##D adj.list <- vector("list", n)
##D weight.list <- vector("list", n)
##D 
##D for (i in 1:n) {
##D   neighbors <- c(i-1, i+1)
##D   # Handle wrap-around for ring structure
##D   neighbors\[neighbors == 0\] <- n
##D   neighbors\[neighbors == n+1\] <- 1
##D 
##D   adj.list\[\[i\]\] <- neighbors
##D   weight.list\[\[i\]\] <- rep(1, length(neighbors))
##D }
##D 
##D # Generate response values with spatial pattern
##D y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
##D 
##D # Apply spectral LOWESS
##D result <- graph.spectral.lowess(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   n.evectors = 5,
##D   verbose = TRUE
##D )
##D 
##D # Plot results
##D plot(y, type="l", col="gray", main="Graph Spectral LOWESS")
##D lines(result$predictions, col="red", lwd=2)
## End(Not run)




### Name: graph.spectrum
### Title: Compute Graph Spectrum
### Aliases: graph.spectrum

### ** Examples

# Create a simple graph adjacency list
graph <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))

# Compute spectrum using R implementation
spec_r <- graph.spectrum(graph, nev = 3, use.R = TRUE)

# Compute spectrum using C implementation
spec_c <- graph.spectrum(graph, nev = 3, use.R = FALSE)

# Get spectrum with Laplacian matrix
spec_lap <- graph.spectrum(graph, return.Laplacian = TRUE)




### Name: grid.critical.points.plot
### Title: Plot Critical Points on Grid
### Aliases: grid.critical.points.plot

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D critical <- find.critical.points(f.grid)
##D grid.critical.points.plot(grid, critical)
## End(Not run)



### Name: harmonic.smoother
### Title: Perform Harmonic Smoothing with Topology Tracking
### Aliases: harmonic.smoother

### ** Examples

## Not run: 
##D # Create a simple grid graph
##D grid.graph <- create.graph.from.grid(10, 10)
##D 
##D # Create noisy function values
##D values <- sin(0.1 * seq_len(100)) + rnorm(100, 0, 0.1)
##D 
##D # Define a region for smoothing (center of the grid)
##D region <- 35:65
##D 
##D # Apply harmonic smoothing with topology tracking
##D result <- harmonic.smoother(
##D   grid.graph$adj.list,
##D   grid.graph$weight.list,
##D   values,
##D   region,
##D   max.iterations = 200,
##D   tolerance = 1e-8,
##D   record.frequency = 5,  # Record every 5 iterations
##D   stability.window = 3,
##D   stability.threshold = 0.05
##D )
##D 
##D # Get the smoothed values
##D smoothed.values <- result$harmonic_predictions
##D 
##D # Plot original vs smoothed values
##D plot(values, type = "l", col = "gray")
##D lines(smoothed.values, col = "red")
##D 
##D # Plot the evolution of topology differences
##D plot(result$topology_differences, type = "l",
##D      xlab = "Iteration", ylab = "Topology Difference")
##D abline(v = result$stable_iteration, col = "blue", lty = 2)
## End(Not run)




### Name: hdbscan.cltr
### Title: Apply HDBSCAN Clustering
### Aliases: hdbscan.cltr

### ** Examples

## Not run: 
##D # Generate sample data
##D set.seed(123)
##D X <- rbind(
##D   matrix(rnorm(100, mean = 0), ncol = 2),
##D   matrix(rnorm(100, mean = 5), ncol = 2)
##D )
##D 
##D # Apply HDBSCAN clustering
##D result <- hdbscan.cltr(X, min.pts = 5:50, soft.K = 20, verbose = TRUE)
##D 
##D # View optimal clustering for Dunn index
##D table(result$dunn.cltr)
## End(Not run)




### Name: hgrid
### Title: Construct a Hierarchical Uniform Grid Around a State Space
### Aliases: hgrid

### ** Examples

## Not run: 
##D # Generate a synthetic 2D state space
##D X <- matrix(runif(200), ncol = 2)
##D # Construct the hierarchical grid with edge length 0.1 and epsilon 0.2
##D grid <- hgrid(X, w = 0.1, epsilon = 0.2)
##D plot(X, col = "red")
##D points(grid, col = "blue", pch = 3)
## End(Not run)



### Name: hist1
### Title: Histogram with customized default settings
### Aliases: hist1

### ** Examples

# Basic usage
data <- rnorm(1000)
hist1(data)

# With custom settings
hist1(data, main = "Normal Distribution", n.breaks = 50, col = "blue")

# Using additional parameters
hist1(data, main = "Custom Histogram", n.breaks = 30, col = "green",
      xlab = "Values", ylab = "Frequency")




### Name: hist2
### Title: Plot two overlapping histograms
### Aliases: hist2

### ** Examples

# Basic usage
x1 <- rnorm(1000, mean = 0, sd = 1)
x2 <- rnorm(1000, mean = 2, sd = 1.5)
hist2(x1, x2, x1.lab = "Group A", x2.lab = "Group B")

# With custom settings
hist2(x1, x2, x1.lab = "Control", x2.lab = "Treatment",
      xlab = "Values", main = "Comparison of Distributions",
      n.x1.breaks = 50, n.x2.breaks = 50)

# Custom colors and no legend
hist2(x1, x2, x1.lab = "A", x2.lab = "B",
      col.x1 = "lightgreen", col.x2 = "lightcoral",
      legend.pos = "none")




### Name: hist3
### Title: Plot three overlapping histograms
### Aliases: hist3

### ** Examples

# Basic usage
x1 <- rnorm(1000, mean = 0, sd = 1)
x2 <- rnorm(1000, mean = 2, sd = 1.2)
x3 <- rnorm(1000, mean = 1, sd = 0.8)
hist3(x1, x2, x3,
      x1.lab = "Group A", x2.lab = "Group B", x3.lab = "Group C")

# With custom settings
hist3(x1, x2, x3,
      x1.lab = "Control", x2.lab = "Treatment 1", x3.lab = "Treatment 2",
      xlab = "Measurements", main = "Three-way Comparison",
      n.x1.breaks = 50, n.x2.breaks = 50, n.x3.breaks = 50)

# Fixed y-axis limit and custom colors
hist3(x1, x2, x3,
      x1.lab = "A", x2.lab = "B", x3.lab = "C",
      ylim.max = 0.5,
      col.x1 = "lightgreen", col.x2 = "lightblue", col.x3 = "lightcoral")




### Name: hor.error.bar
### Title: Add Horizontal Error Bar to Plot
### Aliases: hor.error.bar

### ** Examples

plot(rnorm(10), 1:10, xlim = c(-3, 3))
hor.error.bar(-1, 1, 5, col = "blue")




### Name: identical.vertex.set.weighted.graph.similarity
### Title: Weighted Graph Matching between two graphs with identical vertex
###   sets
### Aliases: identical.vertex.set.weighted.graph.similarity

### ** Examples

# Create adjacency lists (as lists, not matrices)
g1.adj.list <- list(c(2), c(3), c(1))  # vertex 1->2, vertex 2->3, vertex 3->1
g1.weights <- list(c(1), c(2), c(3))   # corresponding weights

g2.adj.list <- list(c(2), c(3), c(1))  # same structure
g2.weights <- list(c(1), c(3), c(2))   # different weights

similarity <- identical.vertex.set.weighted.graph.similarity(
  g1.adj.list, g1.weights, g2.adj.list, g2.weights
)
print(similarity)




### Name: identify_points_within_box
### Title: Identify Points Within a Box
### Aliases: identify_points_within_box

### ** Examples

## Not run: 
##D X <- matrix(runif(100), ncol = 2)
##D L <- c(0.2, 0.2)
##D R <- c(0.5, 0.5)
##D points_within_box <- identify_points_within_box(X, L, R)
## End(Not run)




### Name: init.version.nn.distance.ratio.estimator
### Title: Estimate Relative Entropy Using Nearest Neighbor Distance Ratio
### Aliases: init.version.nn.distance.ratio.estimator

### ** Examples

X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- nn.distance.ratio.estimator(X, Y)
print(result)




### Name: instrumented.gds
### Title: Instrumented Graph Diffusion Smoother
### Aliases: instrumented.gds

### ** Examples

## Not run: 
##D # Example with simple graph
##D graph <- list(c(2), c(1,3), c(2))
##D edge.lengths <- list(c(1), c(1,1), c(1))
##D y.true <- c(0, 1, 0)
##D y.noisy <- y.true + rnorm(3, 0, 0.1)
##D 
##D result <- instrumented.gds(
##D   graph = graph,
##D   edge.lengths = edge.lengths,
##D   y = y.noisy,
##D   y.true = y.true,
##D   n.time.steps = 50,
##D   base.step.factor = 0.5
##D )
## End(Not run)




### Name: itriangle.plot
### Title: itriangle.plot
### Aliases: itriangle.plot

### ** Examples

## Not run: 
##D   library(igraph)
##D 
##D   # Define the inverted triangle plot function
##D   itriangle.plot <- function(coords, v = NULL, params) {
##D       vertex.color <- params("vertex", "color")
##D       if (length(vertex.color) != 1 && !is.null(v)) {
##D           vertex.color <- vertex.color\[v\]
##D       }
##D       vertex.frame.color <- params("vertex", "frame.color")
##D       if (length(vertex.frame.color) != 1 && !is.null(v)) {
##D           vertex.frame.color <- vertex.frame.color\[v\]
##D       }
##D       vertex.frame.width <- params("vertex", "frame.width")
##D       if (length(vertex.frame.width) != 1 && !is.null(v)) {
##D           vertex.frame.width <- vertex.frame.width\[v\]
##D       }
##D       vertex.size <- 1/200 * params("vertex", "size")
##D       if (length(vertex.size) != 1 && !is.null(v)) {
##D           vertex.size <- vertex.size\[v\]
##D       }
##D       vertex.size <- rep(vertex.size, length.out = nrow(coords))
##D 
##D       # Adjust the size for the inverted triangle
##D       side.length <- sqrt(4 * pi / sqrt(3)) * vertex.size / 2
##D 
##D       vertex.frame.color\[vertex.frame.width <= 0\] <- NA
##D       vertex.frame.width\[vertex.frame.width <= 0\] <- 1
##D 
##D       for (i in 1:nrow(coords)) {
##D           x <- coords\[i, 1\]
##D           y <- coords\[i, 2\]
##D           size <- side.length\[i\]
##D           polygon(x + size * c(cos(3*pi/2), cos(5*pi/6), cos(pi/6)),
##D                   y + size * c(sin(3*pi/2), sin(5*pi/6), sin(pi/6)),
##D                   col = vertex.color\[i\],
##D                   border = vertex.frame.color\[i\],
##D                   lwd = vertex.frame.width\[i\])
##D       }
##D   }
##D 
##D   # Add the inverted triangle shape to igraph
##D   add_shape("itriangle", clip = shapes(shape = "circle")$clip,
##D             plot = itriangle.plot)
##D 
##D   # Example graph
##D   g <- make_ring(10)
##D   V(g)$shape <- rep(c("circle", "itriangle"), length.out = vcount(g))
##D   plot(g, vertex.size = 15, vertex.color = "skyblue")
## End(Not run)




### Name: jaccard.index
### Title: Calculate Jaccard index between two sets represented as numeric
###   vectors
### Aliases: jaccard.index

### ** Examples

# Example 1: Basic usage with duplicates
samples.1 <- c(1, 2, 5, 5)
samples.2 <- c(2, 2, 5, 7)
jaccard.index(samples.1, samples.2)  # Returns 0.5 (intersection: {2,5}, union: {1,2,5,7})

# Example 2: Identical sets
jaccard.index(c(1, 2, 3), c(3, 2, 1))  # Returns 1

# Example 3: No overlap
jaccard.index(c(1, 2, 3), c(4, 5, 6))  # Returns 0

# Example 4: One empty set
jaccard.index(c(1, 2, 3), numeric(0))  # Returns 0

# Example 5: Comparing sample indices or cluster memberships
cluster1_samples <- c(1, 3, 5, 7, 9)
cluster2_samples <- c(2, 4, 5, 7, 8)
similarity <- jaccard.index(cluster1_samples, cluster2_samples)
print(paste("Cluster overlap:", round(similarity * 100, 1), "%"))



### Name: jensen.shannon.divergence
### Title: Calculate Jensen-Shannon Divergence Between Two Probability Mass
###   Functions
### Aliases: jensen.shannon.divergence

### ** Examples

# Calculate JSD between two simple distributions
p <- c(0.4, 0.6)
q <- c(0.5, 0.5)
jensen.shannon.divergence(p, q)

# Calculate JSD between distributions of different lengths
p <- c(0.3, 0.7)
q <- c(0.2, 0.3, 0.5)
jensen.shannon.divergence(p, q)  # p will be padded with 0 to match q's length

# Calculate JSD between two categorical distributions
p <- c(0.2, 0.3, 0.5)
q <- c(0.1, 0.4, 0.5)
jensen.shannon.divergence(p, q)




### Name: join.graphs
### Title: Join two graphs at specified vertices
### Aliases: join.graphs

### ** Examples

graph1 <- list(c(2, 3), c(1), c(1))
graph2 <- list(c(2), c(1, 3), c(2))
joined_graph <- join.graphs(graph1, graph2, 2, 1)
print(joined_graph)




### Name: knn_adaptive_mean_shift_gfa
### Title: kNN-Adaptive Mean Shift with Gradient Field Averaging
### Aliases: knn_adaptive_mean_shift_gfa

### ** Examples

set.seed(1)
X <- matrix(rnorm(200), ncol = 2)
out <- knn_adaptive_mean_shift_gfa(
  X, k = 10, density_k = 10, n_steps = 5, step_size = 0.1
)
sapply(out$X_traj, dim)



### Name: kNN.cltr.imputation
### Title: kNN-Based Cluster Imputation
### Aliases: kNN.cltr.imputation

### ** Examples

## Not run: 
##D # Generate sample data with 3 clusters and some noise
##D set.seed(123)
##D X <- rbind(
##D   matrix(rnorm(100, mean = 0), ncol = 2),
##D   matrix(rnorm(100, mean = 5), ncol = 2),
##D   matrix(rnorm(100, mean = c(5, 0)), ncol = 2),
##D   matrix(runif(20, -2, 7), ncol = 2)  # noise points
##D )
##D 
##D # Initial clustering with some points labeled as 0 (noise)
##D cltr <- c(rep(1, 50), rep(2, 50), rep(3, 50), rep(0, 20))
##D names(cltr) <- rownames(X) <- paste0("point_", 1:nrow(X))
##D 
##D # Apply kNN imputation
##D new_cltr <- kNN.cltr.imputation(X, cltr, ref.cltr = 0, K = 10)
##D 
##D # Check results
##D table(Original = cltr, Imputed = new_cltr)
## End(Not run)




### Name: knn.weighted.mean
### Title: Calculates K-Nearest Neighbor Weighted Mean
### Aliases: knn.weighted.mean

### ** Examples

X <- matrix(rnorm(100), ncol = 2)
y <- rnorm(50)
weighted_means <- knn.weighted.mean(X, y, k = 5)




### Name: kullback.leibler.divergence
### Title: Calculate Kullback-Leibler Divergence Between Two Probability
###   Mass Functions
### Aliases: kullback.leibler.divergence

### ** Examples

# Calculate KL divergence between two simple distributions
p <- c(0.4, 0.6)
q <- c(0.5, 0.5)
kullback.leibler.divergence(p, q)

# Using different zero handling methods
p <- c(0.2, 0.8, 0)
q <- c(0.3, 0.6, 0.1)
kullback.leibler.divergence(p, q, zero.handling = "exclude")
kullback.leibler.divergence(p, q, zero.handling = "pseudo")

# Calculate KL divergence between distributions of different lengths
p <- c(0.3, 0.7)
q <- c(0.2, 0.3, 0.5)
kullback.leibler.divergence(p, q)  # p will be padded with 0 to match q's length




### Name: left.asymmetric.bump.fn
### Title: Creates Left Asymmetric Bump Function
### Aliases: left.asymmetric.bump.fn

### ** Examples

## Not run: 
##D x <- seq(-2, 2, length.out = 100)
##D plot(x, AsymmetricBumpFunction(x), type = "l")
## End(Not run)



### Name: lmax.basins
### Title: Identify Basins of Attraction of Local Maxima on a Weighted
###   Graph
### Aliases: lmax.basins

### ** Examples

## Not run: 
##D # Create a weighted graph
##D adj.list <- list(
##D   c(2, 3),      # vertex 1 neighbors
##D   c(1, 3, 4),   # vertex 2 neighbors
##D   c(1, 2, 5),   # vertex 3 neighbors
##D   c(2, 5),      # vertex 4 neighbors
##D   c(3, 4)       # vertex 5 neighbors
##D )
##D 
##D # Edge weights (optional)
##D weight.list <- list(
##D   c(1.0, 0.5),     # weights for edges from vertex 1
##D   c(1.0, 0.8, 1.2), # weights for edges from vertex 2
##D   c(0.5, 0.8, 1.0), # weights for edges from vertex 3
##D   c(1.2, 1.5),      # weights for edges from vertex 4
##D   c(1.0, 1.5)       # weights for edges from vertex 5
##D )
##D 
##D # Values at vertices
##D y <- c(0.5, 0.8, 0.3, 1.2, 1.0)
##D 
##D # Define local maxima
##D lmax.list <- list(
##D   list(lmax = 4, vertices = c(4), label = "max1"),
##D   list(lmax = 5, vertices = c(5), label = "max2")
##D )
##D 
##D # Compute basins
##D basins <- lmax.basins(adj.list, weight.list, y, lmax.list, verbose = TRUE)
## End(Not run)




### Name: lmin.basins
### Title: Identify Basins of Attraction of Local Minima on a Weighted
###   Graph
### Aliases: lmin.basins

### ** Examples

## Not run: 
##D # Create a weighted graph
##D adj.list <- list(
##D   c(2, 3),      # vertex 1 neighbors
##D   c(1, 3, 4),   # vertex 2 neighbors
##D   c(1, 2, 5),   # vertex 3 neighbors
##D   c(2, 5),      # vertex 4 neighbors
##D   c(3, 4)       # vertex 5 neighbors
##D )
##D 
##D # Values at vertices
##D y <- c(0.5, 0.8, 0.3, 1.2, 0.2)
##D 
##D # Define local minima
##D lmin.list <- list(
##D   list(lmin = 5, vertices = c(5), label = "min1"),
##D   list(lmin = 3, vertices = c(3), label = "min2")
##D )
##D 
##D # Compute basins
##D basins <- lmin.basins(adj.list, weight.list = NULL, y, lmin.list)
## End(Not run)




### Name: load.graph.data
### Title: Load Graph Data from RDA Files
### Aliases: load.graph.data

### ** Examples

## Not run: 
##D # Load graphs for k = 5, 10, 15
##D k.values <- c(5, 10, 15)
##D graph.data <- load.graph.data(
##D   k.values,
##D   prefix = "path/to/graphs/graph_k",
##D   suffix = "_NN_graph.rda",
##D   graph.object.name = "S.graph"
##D )
## End(Not run)




### Name: mabilo.plus
### Title: Model-Averaged Locally Weighted Scatterplot Smoothing Plus
###   (MABILO Plus)
### Aliases: mabilo.plus

### ** Examples

# Basic usage
x <- seq(0, 10, length.out = 100)
y <- sin(x) + rnorm(100, 0, 0.1)
fit <- mabilo.plus(x, y)

# With custom parameters
fit2 <- mabilo.plus(x, y, k.min = 5, k.max = 20,
                    model.averaging.strategy = "kernel.weights.only")

# Plot the results
plot(x, y)
lines(fit$x_sorted, fit$ma_predictions, col = "red", lwd = 2)




### Name: mabilo
### Title: Model-Averaged Locally Weighted Scatterplot Smoothing (MABILO)
### Aliases: mabilo

### ** Examples

# Basic usage
x <- seq(0, 10, length.out = 100)
y <- sin(x) + rnorm(100, 0, 0.1)
fit <- mabilo(x, y, k.min = 3, k.max = 10)
plot(x, y)
lines(x, fit$predictions, col = "red")

# With Bayesian bootstrap
fit_bb <- mabilo(x, y, k.min = 3, k.max = 10, n.bb = 100, p = 0.95)
lines(x, fit_bb$cri_L, col = "blue", lty = 2)
lines(x, fit_bb$cri_U, col = "blue", lty = 2)




### Name: mabilog
### Title: Model-Averaged Binary Locally-Weighted Logistic Smoothing
###   (MABILOG)
### Aliases: mabilog

### ** Examples

# Generate binary data with smooth probability structure
set.seed(42)
x <- seq(0, 10, length.out = 200)
true_prob <- plogis(sin(x) - 0.5)
y <- rbinom(length(x), 1, true_prob)

# Fit MABILOG model
fit <- mabilog(x, y, y.true = true_prob, k.min = 5, k.max = 20)

# Plot results
plot(x, y, col = c("red", "blue")[y + 1], pch = 19, cex = 0.5)
lines(fit$x_sorted, fit$predictions, lwd = 2)
lines(x, true_prob, col = "green", lty = 2)
legend("topright", c("Data (y=0)", "Data (y=1)", "Fitted", "True"),
       col = c("red", "blue", "black", "green"),
       pch = c(19, 19, NA, NA), lty = c(NA, NA, 1, 2))

# With bootstrap confidence intervals
## No test: 
fit_bb <- mabilog(x, y, k.min = 5, k.max = 20, n.bb = 100)
plot(x, y, col = c("red", "blue")[y + 1], pch = 19, cex = 0.5)
polygon(c(fit_bb$x_sorted, rev(fit_bb$x_sorted)),
        c(fit_bb$cri_L, rev(fit_bb$cri_U)),
        col = "gray80", border = NA)
lines(fit_bb$x_sorted, fit_bb$predictions, lwd = 2)
## End(No test)




### Name: mae.plot
### Title: Plot Median Absolute Error
### Aliases: mae.plot

### ** Examples

## Not run: 
##D mae.mean <- c(NA, 0.5, 0.4, 0.35, 0.33, 0.32)
##D mae.mad <- c(NA, 0.1, 0.08, 0.07, 0.06, 0.06)
##D mae.plot(mae.mean, mae.mad)
## End(Not run)




### Name: maelog
### Title: Model-Averaged Local Logistic Regression
### Aliases: maelog

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 200
##D x <- seq(0, 1, length.out = n)
##D true.prob <- plogis(10 * (x - 0.5))  # Logistic function
##D y <- rbinom(n, 1, true.prob)
##D 
##D # Fit model with automatic bandwidth selection
##D fit1 <- maelog(x, y)
##D 
##D # Plot results
##D plot(x, y, col = rgb(0, 0, 0, 0.5), pch = 16,
##D      main = "Local Logistic Regression",
##D      xlab = "x", ylab = "Probability")
##D lines(x, fit1$predictions, col = "blue", lwd = 2)
##D lines(x, true.prob, col = "red", lwd = 2, lty = 2)
##D legend("topleft", c("Fitted", "True"),
##D        col = c("blue", "red"), lwd = 2, lty = c(1, 2))
##D 
##D # Fit with quadratic local models
##D fit2 <- maelog(x, y, fit.quadratic = TRUE, cv.folds = 5)
##D 
##D # Fit with fixed bandwidth
##D fit3 <- maelog(x, y, pilot.bandwidth = 0.1)
##D 
##D # Compare different kernels
##D fit.tricube <- maelog(x, y, kernel = 7)
##D fit.gauss <- maelog(x, y, kernel = 5)
##D 
##D # Larger example with cross-validation
##D n <- 1000
##D x <- runif(n, -2, 2)
##D prob <- plogis(2 * sin(pi * x))
##D y <- rbinom(n, 1, prob)
##D 
##D # Compare linear vs quadratic with 10-fold CV
##D fit.linear <- maelog(x, y, fit.quadratic = FALSE, cv.folds = 10)
##D fit.quad <- maelog(x, y, fit.quadratic = TRUE, cv.folds = 10)
##D 
##D # Plot bandwidth selection results
##D par(mfrow = c(1, 2))
##D plot(fit.linear$candidate.bandwidths, fit.linear$mean.errors,
##D      type = "l", xlab = "Bandwidth", ylab = "CV Error",
##D      main = "Linear Model")
##D abline(v = fit.linear$opt.bw, col = "red", lty = 2)
##D 
##D plot(fit.quad$candidate.bandwidths, fit.quad$mean.errors,
##D      type = "l", xlab = "Bandwidth", ylab = "CV Error",
##D      main = "Quadratic Model")
##D abline(v = fit.quad$opt.bw, col = "red", lty = 2)
## End(Not run)




### Name: magelo
### Title: Model Averaged Grid-based Epsilon LOwess (MAGELO)
### Aliases: magelo

### ** Examples

x <- seq(0, 10, length.out = 100)
y <- sin(x) + rnorm(100, 0, 0.1)
fit <- magelo(x, y, degree = 1)
plot(x, y)
lines(fit$xgrid, fit$gpredictions, col = "red")




### Name: magelog
### Title: Grid-based Model Averaged Bandwidth Logistic Regression
### Aliases: magelog

### ** Examples

## Not run: 
##D # Generate example data with logistic relationship
##D set.seed(123)
##D n <- 200
##D x <- seq(0, 1, length.out = n)
##D # True probability function
##D p_true <- 1/(1 + exp(-(x - 0.5)*10))
##D y <- rbinom(n, 1, p_true)
##D 
##D # Fit model with automatic bandwidth selection
##D fit <- magelog(x, y, grid.size = 100, fit.quadratic = FALSE, cv.folds = 5)
##D 
##D # Plot results
##D plot(x, y, pch = 19, col = adjustcolor("black", 0.5),
##D      xlab = "x", ylab = "Probability",
##D      main = "Local Logistic Regression")
##D # Add true probability curve
##D lines(x, p_true, col = "gray", lwd = 2, lty = 2)
##D # Add fitted curve
##D lines(fit$x.grid, fit$bw.grid.predictions[, fit$opt.brier.bw.idx],
##D       col = "red", lwd = 2)
##D legend("topleft", c("True probability", "Fitted curve", "Data"),
##D        col = c("gray", "red", "black"), lty = c(2, 1, NA),
##D        pch = c(NA, NA, 19), lwd = c(2, 2, NA))
##D 
##D # Examine bandwidth selection
##D plot(fit$bws, fit$mean.brier.errors, type = "b",
##D      xlab = "Bandwidth", ylab = "Cross-validation Brier Score",
##D      main = "Bandwidth Selection")
##D abline(v = fit$bws[fit$opt.brier.bw.idx], col = "red", lty = 2)
## End(Not run)




### Name: map.S.to.X
### Title: Map Sample IDs to 3D Space
### Aliases: map.S.to.X

### ** Examples

## Not run: 
##D X <- matrix(rnorm(300), ncol = 3)
##D rownames(X) <- paste0("Sample", 1:100)
##D S <- paste0("Sample", c(1, 5, 10, 15, 20))
##D 
##D map.S.to.X(S, X, radius = 0.075, col = 'red')
## End(Not run)




### Name: matrix.format
### Title: Format Numbers in a Matrix with Specified Decimal Places
### Aliases: matrix.format

### ** Examples

m <- matrix(c(0.7, 1.235, 0.1, 0.8876), nrow = 2)
matrix.format(m)               # default: 2 digits
matrix.format(m, digits = 3)   # 3 digits
matrix.format(m, nsmall = 3)   # at least 3 decimal places




### Name: meanshift.data.smoother
### Title: Mean Shift Data Smoother
### Aliases: meanshift.data.smoother

### ** Examples

## Not run: 
##D # Generate sample data
##D set.seed(123)
##D X <- matrix(rnorm(200), ncol = 2)
##D 
##D # Basic mean shift smoothing
##D result1 <- meanshift.data.smoother(X, k = 5)
##D 
##D # Adaptive mean shift with gradient field averaging
##D result2 <- meanshift.data.smoother(X, k = 5, method = "adaptive")
##D 
##D # With momentum
##D result3 <- meanshift.data.smoother(X, k = 5, method = "adaptive_momentum",
##D                                   momentum = 0.95)
##D 
##D # Plot the trajectory
##D plot(result1$median.kdistances, type = 'l',
##D      xlab = 'Step', ylab = 'Median k-distance')
##D abline(v = result1$opt.step, col = 'red', lty = 2)
##D 
##D # Compare original and smoothed data
##D par(mfrow = c(1, 2))
##D plot(X, main = "Original", pch = 19, cex = 0.5)
##D plot(result1$dX, main = "Smoothed", pch = 19, cex = 0.5)
## End(Not run)




### Name: minh.limit
### Title: Find Minimum Hop Limit for Path Existence
### Aliases: minh.limit

### ** Examples

## Not run: 
##D pgs <- create.path.graph.series(graph, edge.lengths, h.values = 1:5)
##D min.h <- minh.limit(pgs, from = 1, to = 5)
##D if (!is.null(min.h)) {
##D   cat("Minimum hops needed:", min.h, "\n")
##D }
## End(Not run)




### Name: minmax.normalize
### Title: Min-Max Normalization
### Aliases: minmax.normalize

### ** Examples

x <- c(1, 2, 3, 4, 5)
minmax.normalize(x, 0, 1)  # Normalizes x to the range \eqn{[0, 1]}




### Name: morse.analysis.plot
### Title: Plot Morse Analysis
### Aliases: morse.analysis.plot

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D critical <- find.critical.points(f.grid)
##D traj <- compute.gradient.trajectory(15, 20, f.grid)
##D morse.analysis.plot(grid, f.grid, traj, critical)
## End(Not run)



### Name: morse.smale.cells.plot
### Title: Plot Morse-Smale Cells
### Aliases: morse.smale.cells.plot

### ** Examples

## Not run: 
##D grid <- create.grid(30)
##D f <- function(x, y) sin(3*x) * cos(3*y)
##D f.grid <- evaluate.function.on.grid(f, grid)
##D complex <- create.morse.smale.complex(f.grid)
##D morse.smale.cells.plot(grid, complex, f.grid)
## End(Not run)



### Name: morse.smale.complex.lmin.lmax.summary
### Title: Process Complete Morse-Smale Complex Analysis Pipeline
### Aliases: morse.smale.complex.lmin.lmax.summary

### ** Examples

## Not run: 
##D results <- morse.smale.complex.lmin.lmax.summary(
##D     MS.res = smoothed.condE.MS.res,
##D     state.space = asv.matrix,
##D     taxonomy = asv.taxonomy,
##D     rel.condE = rel.smoothed.condE,
##D     outcome.name = "smoothed conditional expectation"
##D )
## End(Not run)



### Name: morse.smale.complex.plot.from.critical
### Title: Plot Morse-Smale Complex from Critical Points
### Aliases: morse.smale.complex.plot.from.critical

### ** Examples

## Not run: 
##D # Example with known critical points
##D grid <- create.grid(50)
##D mixture <- list(
##D   f = function(x, y) x^2 - y^2,  # Simple saddle
##D   gradient = function(x, y) c(2*x, -2*y)
##D )
##D critical_points <- list(
##D   maxima = matrix(c(0.8, 0.2), nrow = 1),
##D   minima = matrix(c(0.2, 0.8), nrow = 1),
##D   saddles = matrix(c(0.5, 0.5), nrow = 1)
##D )
##D morse.smale.complex.plot.from.critical(critical_points, mixture, grid)
## End(Not run)



### Name: morse.smale.complex.plot
### Title: Plot Morse-Smale Complex from Function
### Aliases: morse.smale.complex.plot

### ** Examples

## Not run: 
##D # Create a simple function with two critical points
##D grid <- create.grid(50)
##D mixture <- list(
##D   f = function(x, y) sin(3*pi*x) * cos(3*pi*y),
##D   gradient = function(x, y) {
##D     c(3*pi*cos(3*pi*x)*cos(3*pi*y),
##D       -3*pi*sin(3*pi*x)*sin(3*pi*y))
##D   }
##D )
##D morse.smale.complex.plot(grid, mixture, n_sample_points = 10)
## End(Not run)



### Name: mutual.information
### Title: Estimate Mutual Information Between Dataset Membership and
###   Features
### Aliases: mutual.information

### ** Examples

# Example 1: Datasets with different means
set.seed(123)
X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- mutual.information(X, Y)
print(paste("MI between datasets:", round(result, 4), "nats"))

# Example 2: Identical datasets (MI should be ~0)
X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000), ncol = 2)
result <- mutual.information(X, Y)
print(paste("MI for identical distributions:", round(result, 4)))

# Example 3: Feature-wise contribution
X <- matrix(rnorm(500), ncol = 5)
Y <- X
Y[,3] <- Y[,3] + 2  # Only change feature 3
# Compute MI for each feature separately
## Not run: 
##D feature_mi <- sapply(1:5, function(i) {
##D   mutual.information(X[,i,drop=FALSE], Y[,i,drop=FALSE])
##D })
##D barplot(feature_mi, names.arg = paste("Feature", 1:5),
##D         main = "Feature-wise MI Contribution")
## End(Not run)



### Name: nada.graph.spectral.lowess
### Title: Non-adaptive Local Regression on Graphs Using Spectral Embedding
### Aliases: nada.graph.spectral.lowess

### ** Examples

## Not run: 
##D # Create a simple graph with 100 vertices
##D n <- 100
##D set.seed(123)
##D 
##D # Create a ring graph
##D adj.list <- vector("list", n)
##D weight.list <- vector("list", n)
##D 
##D for (i in 1:n) {
##D   neighbors <- c(i-1, i+1)
##D   # Handle wrap-around for ring structure
##D   neighbors\[neighbors == 0\] <- n
##D   neighbors\[neighbors == n+1\] <- 1
##D 
##D   adj.list\[\[i\]\] <- neighbors
##D   weight.list\[\[i\]\] <- rep(1, length(neighbors))
##D }
##D 
##D # Generate response values with spatial pattern
##D y <- sin(2*pi*(1:n)/n) + rnorm(n, 0, 0.2)
##D 
##D # Apply spectral LOWESS
##D result <- nada.graph.spectral.lowess(
##D   adj.list = adj.list,
##D   weight.list = weight.list,
##D   y = y,
##D   n.evectors = 5,
##D   verbose = TRUE
##D )
##D 
##D # Plot results
##D plot(y, type="l", col="gray", main="Graph Spectral LOWESS")
##D lines(result$predictions, col="red", lwd=2)
## End(Not run)




### Name: nerve.cx.spectral.filter
### Title: Apply Spectral Filtering to a Nerve Complex
### Aliases: nerve.cx.spectral.filter

### ** Examples

## Not run: 
##D # Generate 2D points
##D set.seed(123)
##D coords <- matrix(runif(200), ncol = 2)
##D 
##D # Create a smooth function with noise
##D f_true <- function(x) sin(2*pi*x[1]) * cos(2*pi*x[2])
##D y_true <- apply(coords, 1, f_true)
##D y_noisy <- y_true + rnorm(length(y_true), 0, 0.2)
##D 
##D # Create nerve complex
##D complex <- create.nerve.complex(coords, k = 8, max.dim = 2)
##D complex <- set.complex.function.values(complex, y_noisy)
##D 
##D # Apply spectral filtering with default parameters
##D result <- nerve.cx.spectral.filter(complex, y_noisy)
##D 
##D # Apply with custom parameters emphasizing higher dimensions
##D result2 <- nerve.cx.spectral.filter(
##D   complex, y_noisy,
##D   laplacian.type = "NORMALIZED",
##D   filter.type = "GAUSSIAN",
##D   dim.weights = c(1.0, 0.5, 0.25),
##D   n.candidates = 150,
##D   verbose = TRUE
##D )
##D 
##D # Compare results
##D mse1 <- mean((result$predictions - y_true)^2)
##D mse2 <- mean((result2$predictions - y_true)^2)
##D cat("MSE (default):", mse1, "\n")
##D cat("MSE (custom):", mse2, "\n")
## End(Not run)




### Name: nerve.graph
### Title: Constructs the 1-skeleton (nerve graph) of the nerve simplicial
###   complex associate with a covering of some finite set
### Aliases: nerve.graph

### ** Examples

covering.list <- list(c(1, 2), c(2, 3), c(1, 3, 4))
nerve.graph(covering.list, n.cores = 2)



### Name: nn.distance.ratio.estimator
### Title: Estimate Relative Entropy Using Nearest Neighbor Distance Ratio
### Aliases: nn.distance.ratio.estimator

### ** Examples

X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- nn.distance.ratio.estimator(X, Y)
print(result)




### Name: normalize.and.inv.logit.fn
### Title: Transform Continuous Values to Probabilities via Normalization
###   and Logistic Function
### Aliases: normalize.and.inv.logit.fn

### ** Examples

# Transform smoothed values into probabilities
y.smooth <- c(-1, 0, 0.5, 1, 2)
probs <- normalize.and.inv.logit.fn(y.smooth)

# Use custom range for normalization
probs.wide <- normalize.and.inv.logit.fn(y.smooth, y.min = -5, y.max = 5)




### Name: normalize.and.inv.logit
### Title: Normalize and Apply Inverse Logit Transformation
### Aliases: normalize.and.inv.logit

### ** Examples

## Not run: 
##D x <- rnorm(10)
##D normalize.and.inv.logit(x)  # Normalizes x and applies inverse logit
## End(Not run)



### Name: o.inv.fn
### Title: Computes the inverse of a permutation vector generated by
###   order()
### Aliases: o.inv.fn

### ** Examples

x <- c(3, 1, 4, 2)
o <- order(x)  # Returns c(2, 4, 1, 3)
o.inv <- o.inv.fn(o)  # Returns c(3, 1, 4, 2)
# Verify: x[o][o.inv] equals x




### Name: overlap.coefficient
### Title: Calculate the Overlap Coefficient Between Two Numeric Vectors
### Aliases: overlap.coefficient

### ** Examples

overlap.coefficient(c(1, 2, 3), c(1, 2, 3))     # Returns 1
overlap.coefficient(c(1, 2), c(1, 2, 3, 4))     # Returns 1
overlap.coefficient(c(1, 2, 2, 3), c(1, 2, 4))  # Returns 0.67
overlap.coefficient(c(1, 2), c(3, 4))           # Returns 0
overlap.coefficient(numeric(0), c(1, 2))        # Returns 0




### Name: parameterize.circular.graph
### Title: Parameterize a Circular Graph Structure for Biological Cycle
###   Analysis
### Aliases: parameterize.circular.graph

### ** Examples

# Example 1: Simple circular graph (like cell cycle progression)
n <- 8
adj.list <- lapply(seq_len(n), function(i) {
  c(if (i == n) 1 else i + 1,  # next vertex
    if (i == 1) n else i - 1)  # previous vertex
})
weight.list <- lapply(adj.list, function(adj) rep(1.0, length(adj)))

# Get circular parameterization
result <- parameterize.circular.graph(adj.list, weight.list, TRUE)

# Display angles (positions along the cycle)
print(round(result$angles, 2))

# Example 2: Biological application - genes with circular expression pattern
# Simulate a gene similarity network where genes are connected if their
# expression patterns are similar (simplified example)
n_genes <- 12
# Create connections based on proximity in the cycle
adj.list <- lapply(seq_len(n_genes), function(i) {
  # Connect to 2 neighbors on each side to simulate local similarity
  neighbors <- c((i - 2):(i - 1), (i + 1):(i + 2))
  neighbors <- ((neighbors - 1) %% n_genes) + 1
  neighbors[neighbors != i]  # Remove self-connections
})
# Weights decrease with distance in the cycle
weight.list <- lapply(seq_len(n_genes), function(i) {
  neighbors <- adj.list[[i]]
  weights <- sapply(neighbors, function(j) {
    dist <- min(abs(i - j), n_genes - abs(i - j))
    exp(-dist/2)  # Exponential decay
  })
  weights
})

result <- parameterize.circular.graph(adj.list, weight.list)

# Plot genes positioned by their inferred cycle position
plot(cos(result$angles), sin(result$angles),
     xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
     pch = 19, xlab = "x", ylab = "y", asp = 1,
     main = "Gene Positions in Expression Cycle")

# Add gene labels
text(1.1 * cos(result$angles), 1.1 * sin(result$angles),
     labels = paste0("G", seq_len(n_genes)), cex = 0.8)

# The angles can be interpreted as positions in the biological cycle
# For cell cycle: 0 = G1/S, pi/2 = S, pi = G2, 3*pi/2 = M phase




### Name: partition.of.unity.1D
### Title: Creates a Partition of Unity in 1D
### Aliases: partition.of.unity.1D

### ** Examples

## Not run: 
##D x.center <- seq(1, 9, by = 2)
##D partition <- partition.of.unity.1D(x.center)
##D plot(seq(0, 10, length.out = 400), partition[,1], type = "l")
## End(Not run)



### Name: path.dist
### Title: Distances from the first vertex of a path to each consecutive
###   vertex along the given path normalized so that the total distance is
###   1
### Aliases: path.dist

### ** Examples

# Create a simple vertex matrix
V <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol=2, byrow=TRUE)
path <- c(1, 2, 3, 4)
distances <- path.dist(path, V)
# distances will be c(0.0, 0.333..., 0.666..., 1.0)




### Name: path.fn
### Title: Animate Path Along Trajectory
### Aliases: path.fn

### ** Examples

## Not run: 
##D # This function is typically used within an animation loop
##D Epath <- matrix(rnorm(30), ncol = 3)
##D sphere.id <- NULL
##D for (i in 2:nrow(Epath)) {
##D   path.fn(i, Epath)
##D   Sys.sleep(0.1)
##D }
## End(Not run)




### Name: path.length
### Title: Computes the length of a path specified by a matrix X
### Aliases: path.length

### ** Examples

# Create a simple path in 2D
path_points <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol=2, byrow=TRUE)
length <- path.length(path_points)
# length will be 3 (three unit-length segments)

# 3D path example
path_3d <- matrix(c(0,0,0, 1,0,0, 1,1,0, 1,1,1), ncol=3, byrow=TRUE)
length_3d <- path.length(path_3d)




### Name: pca.variance.analysis
### Title: Perform PCA Variance Analysis
### Aliases: pca.variance.analysis

### ** Examples

data <- matrix(rnorm(1000*50), ncol=50)
result <- pca.variance.analysis(data)
print(result$pc.95)




### Name: pearson.wcor
### Title: Pearson Weighted Correlation Coefficient
### Aliases: pearson.wcor

### ** Examples

x <- rnorm(100)
y <- rnorm(100)
w <- runif(100)
pearson.wcor(x, y, w)




### Name: perform.harmonic.smoothing
### Title: Perform Harmonic Smoothing on Graph Function Values
### Aliases: perform.harmonic.smoothing

### ** Examples

## Not run: 
##D # Create a simple grid graph
##D grid.graph <- create.graph.from.grid(10, 10)
##D 
##D # Create noisy function values
##D values <- sin(0.1 * seq_len(100)) + rnorm(100, 0, 0.1)
##D 
##D # Define a region for smoothing (center of the grid)
##D region <- 35:65
##D 
##D # Apply harmonic smoothing
##D smoothed.values <- perform.harmonic.smoothing(
##D   grid.graph$adj.list,
##D   grid.graph$weight.list,
##D   values,
##D   region,
##D   max.iterations = 200,
##D   tolerance = 1e-8
##D )
##D 
##D # Plot original vs smoothed values
##D plot(values, type = "l", col = "gray")
##D lines(smoothed.values, col = "red")
## End(Not run)




### Name: pgmalo
### Title: Path Graph Model Averaging Local Linear Model
### Aliases: pgmalo

### ** Examples

## No test: 
# Create a simple chain graph
n <- 50
neighbors <- vector("list", n)
edge_lengths <- vector("list", n)

# Build chain structure
for(i in 1:n) {
  if(i == 1) {
    neighbors[[i]] <- 2L
    edge_lengths[[i]] <- 1.0
  } else if(i == n) {
    neighbors[[i]] <- (n-1L)
    edge_lengths[[i]] <- 1.0
  } else {
    neighbors[[i]] <- c(i-1L, i+1L)
    edge_lengths[[i]] <- c(1.0, 1.0)
  }
}

# Generate smooth response with noise
x_pos <- seq(0, 1, length.out = n)
y <- sin(2 * pi * x_pos) + rnorm(n, 0, 0.1)

# Fit model with small h range for speed
  fit <- pgmalo(neighbors, edge_lengths, y, h.max = 10, n.CVs = 10)
  summary(fit)
  print(fit)
## End(No test)




### Name: plaplace
### Title: Laplace Distribution Function
### Aliases: plaplace

### ** Examples

x <- seq(-5, 5, by = 0.1)
y <- plaplace(x, location = 0, scale = 1)
plot(x, y, type = "l", main = "Laplace CDF")




### Name: plot.assoc0
### Title: Plot Method for assoc0 Objects
### Aliases: plot.assoc0

### ** Examples

## Not run: 
##D # Generate example data and run test
##D set.seed(123)
##D n <- 200
##D x <- runif(n)
##D y <- sin(2*pi*x) + rnorm(n, sd = 0.3)
##D result <- fassoc0.test(x, y, n.cores = 2, n.perms = 1000, plot.it = FALSE)
##D 
##D # Create different plots
##D plot(result, plot = "Exy")
##D plot(result, plot = "d1hist")
##D plot(result, plot = "bc.d1hist")
## End(Not run)




### Name: plot.assoc1
### Title: Plot Method for assoc1 Objects
### Aliases: plot.assoc1

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D n <- 200
##D x <- runif(n)
##D y <- sin(4*pi*x) + rnorm(n, sd = 0.2)
##D result <- fassoc1.test(x, y, n.cores = 2, n.perms = 500, plot.it = FALSE)
##D 
##D # Create different plots
##D plot(result, plot = "Exy")
##D plot(result, plot = "dExy")
##D plot(result, plot = "d1hist")
## End(Not run)




### Name: plot.box.tiling
### Title: Plots boxes of a box tiling of a 2D state space
### Aliases: plot.box.tiling

### ** Examples

## Not run: 
##D w <- 1
##D L <- c(0, 0)
##D R <- c(2, 3)
##D grid_boxes <- create.ED.boxes(w, L, R)
##D plot.boxes(grid_boxes, L, R)
## End(Not run)




### Name: plot.chain.with.path
### Title: Plot a Chain Graph with Optional Highlighted Paths
### Aliases: plot.chain.with.path

### ** Examples

# Create a simple chain graph
adj_list <- list(
  c(2),    # Vertex 1 connected to 2
  c(1, 3), # Vertex 2 connected to 1 and 3
  c(2)     # Vertex 3 connected to 2
)
weight_list <- list(
  c(1),    # Weight for edge 1-2
  c(1, 2), # Weights for edges 2-1 and 2-3
  c(2)     # Weight for edge 3-2
)
path_data <- list(
  list(vertices = c(1, 2, 3), ref_vertex = 1)
)

chain.with.path.obj <- chain.with.path(adj_list, weight_list, path_data)
plot(chain.with.path.obj)




### Name: plot.circular_parameterization
### Title: Plot Method for circular_parameterization Objects
### Aliases: plot.circular_parameterization

### ** Examples

# Create example graph
n <- 8
adj.list <- lapply(seq_len(n), function(i) {
  c(if (i == n) 1 else i + 1, if (i == 1) n else i - 1)
})
weight.list <- lapply(adj.list, function(adj) rep(1.0, length(adj)))
result <- parameterize.circular.graph(adj.list, weight.list)

# Basic plot
plot(result)

# Plot with edges
plot(result, adj.list = adj.list,
     main = "Circular Graph Visualization")




### Name: plot.gaussian_mixture_data
### Title: Plot Method for Gaussian Mixture Data Objects
### Aliases: plot.gaussian_mixture_data

### ** Examples

## Not run: 
##D # Create a 2D Gaussian mixture
##D gm <- get.gaussian.mixture.2d(n.components = 3, sd.value = 0.1)
##D 
##D # Sample points from the mixture
##D data <- sample(gm, n = 500, sampling.method = "random", noise.sd = 0.05)
##D 
##D # Basic 2D scatter plot
##D plot(data)
##D 
##D # Contour plot of underlying function
##D plot(data, type = "contour", nlevels = 20)
##D 
##D # View individual components
##D plot(data, type = "components")
##D 
##D # Compare true vs noisy values
##D plot(data, type = "compare")
##D 
##D # 3D visualization (requires rgl)
##D plot(data, type = "surface", grid_size = 75)
## End(Not run)




### Name: plot.geodesic_stats
### Title: Plot Geodesic Statistics
### Aliases: plot.geodesic_stats

### ** Examples

## Not run: 
##D stats <- compute.geodesic.stats(adj.list, weight.list)
##D plot.geodesic.stats(stats, "all")
## End(Not run)




### Name: plot.ggraph
### Title: Plot a ggraph Object
### Aliases: plot.ggraph

### ** Examples

## Not run: 
##D adj_list <- list(c(2, 3), c(1, 3), c(1, 2))
##D weight_list <- list(c(1, 1), c(1, 1), c(1, 1))
##D g <- ggraph(adj_list, weight_list)
##D plot.res <- plot(g)
##D plot(g, vertex.size = 5, vertex.label = c("A", "B", "C"))
##D plot(g, layout = "in_circle")
##D plot(g, vertex.size = 2, use.saved.layout = plot.res$layout)
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
##D   plot(g, dim = 3)
##D }
## End(Not run)



### Name: plot.graph_kernel_smoother
### Title: Plot Method for graph_kernel_smoother Objects
### Aliases: plot.graph_kernel_smoother

### ** Examples

# See examples in graph.kernel.smoother()




### Name: plot.graphMScx
### Title: Plot Method for graphMScx Objects
### Aliases: plot.graphMScx

### ** Examples

## Not run: 
##D # Create a simple graph
##D adj_list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
##D weight_list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
##D g <- ggraph(adj_list, weight_list)
##D 
##D # Define function values on vertices
##D Ey <- c(0.5, 0.3, 0.8, 0.2)
##D 
##D # Compute Morse-Smale complex
##D mscx <- graph.MS.cx(g, Ey)
##D 
##D # Plot the Morse-Smale complex graph
##D plot(mscx)
##D 
##D # Show branches of attraction for local minima
##D plot(mscx, type = "lmin_BoA")
##D 
##D # Plot a specific cell
##D plot(mscx, type = "MS_graphs", i = 2)
##D 
##D # Show the entire graph with first cell highlighted
##D plot(mscx, type = "MS_cell", i = 1)
## End(Not run)




### Name: plot.harmonic_smoother
### Title: Plot Method for Harmonic Smoother Results
### Aliases: plot.harmonic_smoother

### ** Examples

## Not run: 
##D # After running harmonic.smoother()
##D result <- harmonic.smoother(adj.list, weight.list, values, region)
##D 
##D # Plot topology evolution
##D plot(result, type = "topology")
##D 
##D # Plot extrema counts
##D plot(result, type = "extrema")
##D 
##D # Plot smoothed values
##D plot(result, type = "values")
## End(Not run)




### Name: plot.IkNNgraphs
### Title: Plot Diagnostics for Intersection k-NN Graph Analysis
### Aliases: plot.IkNNgraphs

### ** Examples

## Not run: 
##D # Create sample data for IkNNgraphs analysis
##D set.seed(123)
##D res <- list(
##D   k.values = 3:10,
##D   edit.distances = c(25, 18, 14, 12, 10, 9, 8, 7),
##D   n.edges.in.pruned.graph = c(50, 80, 100, 115, 125, 132, 138, 142),
##D   js.div = c(0.5, 0.3, 0.2, 0.15, 0.12, 0.1, 0.09, 0.08),
##D   edit.distances.breakpoint = 5,
##D   n.edges.in.pruned.graph.breakpoint = 6,
##D   js.div.breakpoint = 5.5
##D )
##D 
##D # Plot diagnostics with default settings
##D   plot.IkNNgraphs(res)
##D 
##D   # Plot only edit distances and edge counts
##D   plot.IkNNgraphs(res, diags = c("edist", "edge"))
## End(Not run)




### Name: plot.instrumented.gds
### Title: Plot Instrumented GDS Results
### Aliases: plot.instrumented.gds

### ** Examples

## Not run: 
##D # Create example data
##D graph <- list(c(2), c(1,3), c(2))
##D edge.lengths <- list(c(1), c(1,1), c(1))
##D y.true <- c(0, 1, 0)
##D y.noisy <- y.true + rnorm(3, 0, 0.1)
##D 
##D # Run instrumented GDS
##D result <- instrumented.gds(
##D   graph = graph,
##D   edge.lengths = edge.lengths,
##D   y = y.noisy,
##D   y.true = y.true,
##D   n.time.steps = 50,
##D   base.step.factor = 0.5
##D )
##D 
##D # Default multi-panel plot
##D plot(result)
##D 
##D # Single panel with all metrics
##D plot(result, layout = "single")
##D 
##D # Single panel with normalization
##D plot(result, layout = "single", normalize = TRUE)
##D 
##D # Custom colors for multi-panel
##D plot(result, col = c("blue", "red", "darkgreen"), lwd = 3)
##D 
##D # Plot only specific metrics in single panel
##D plot(result, metrics = c("snr", "mad"), layout = "single",
##D      col = c("navy", "darkred"), lty = c(1, 2))
## End(Not run)




### Name: plot.kh.matrix
### Title: Plot a k-h Matrix
### Aliases: plot.kh.matrix

### ** Examples

## Not run: 
##D # Create a kh.matrix object
##D kh.mat <- matrix(c(0,1,1,0,1,0,1,1,0), nrow = 3)
##D existing.k <- c(10, 20, 30)
##D h.values <- c(1, 2, 3)
##D khm <- kh.matrix(kh.mat, existing.k, h.values, id = "LM1")
##D 
##D # Basic plot
##D plot(khm)
##D 
##D # With custom colors and a legend
##D plot(khm,
##D      color.palette = c("lightblue", "darkred"),
##D      with.legend = TRUE)
##D 
##D # You can also override parameters from the object
##D plot(khm, main = "Custom Title", xlab = "k values")
## End(Not run)




### Name: plot.local_extrema
### Title: Plot Local Extrema Detection Results
### Aliases: plot.local_extrema

### ** Examples

# Create example data
adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
y <- c(1, 3, 2, 5, 1)

# Detect and plot maxima
maxima <- detect.local.extrema(adj.list, weight.list, y, 2, 2)
plot(maxima)




### Name: plot.mabilo_plus
### Title: Plot Method for Mabilo Plus Objects
### Aliases: plot.mabilo_plus

### ** Examples

## Not run: 
##D # Generate example data
##D x <- seq(0, 10, length.out = 100)
##D y <- sin(x) + rnorm(100, 0, 0.1)
##D 
##D # Fit mabilo.plus model
##D fit <- mabilo.plus(x, y)
##D 
##D # Basic fit plot
##D plot(fit)
##D 
##D # Plot both SM and MA predictions
##D plot(fit, type = "fit", predictions.type = "both")
##D 
##D # Diagnostic plot
##D plot(fit, type = "diagnostic", diagnostic.type = "both")
##D 
##D # Residual plots
##D plot(fit, type = "residuals", predictions.type = "ma")
##D plot(fit, type = "residuals.hist", predictions.type = "both")
## End(Not run)




### Name: plot.mabilo
### Title: Plot Method for Mabilo Objects
### Aliases: plot.mabilo

### ** Examples

# Generate example data
x <- seq(0, 10, length.out = 100)
y <- sin(x) + rnorm(100, 0, 0.1)

# Fit mabilo model
fit <- mabilo(x, y, n.bb = 100)

# Basic fit plot with credible intervals
plot(fit)

# Diagnostic plot showing errors
plot(fit, type = "diagnostic")

# Residual plot
plot(fit, type = "residuals")




### Name: plot.maelog
### Title: Plot Method for maelog Objects
### Aliases: plot.maelog

### ** Examples

# Generate example data
set.seed(123)
n <- 200
x <- runif(n)
y <- rbinom(n, 1, plogis(10 * (x - 0.5)))

# Fit with fixed bandwidth
fit1 <- maelog(x, y, pilot.bandwidth = 0.1)

# Plot fitted values
plot(fit1)

## Not run: 
##D # Plot residuals
##D plot(fit1, which = 3)
##D 
##D # Plot QQ plot
##D plot(fit1, which = 4)
##D 
##D # Fit with automatic bandwidth selection
##D fit2 <- maelog(x, y, pilot.bandwidth = -1, n.bws = 30)
##D 
##D # Show all plots
##D par(mfrow = c(2, 2))
##D plot(fit2, which = 1)  # Fitted values
##D plot(fit2, which = 2)  # Bandwidth selection
##D plot(fit2, which = 3)  # Residuals
##D plot(fit2, which = 4)  # QQ plot
##D par(mfrow = c(1, 1))
## End(Not run)




### Name: plot.morse.smale.trajectories.from.point
### Title: Plot Morse-Smale Trajectories from Point
### Aliases: plot.morse.smale.trajectories.from.point

### ** Examples

## Not run: 
##D grid <- create.grid(50)
##D mixture <- list(
##D   f = function(x, y) -((x-0.3)^2 + (y-0.7)^2) - 2*((x-0.7)^2 + (y-0.3)^2),
##D   gradient = function(x, y) {
##D     c(-2*(x-0.3) - 4*(x-0.7), -2*(y-0.7) - 4*(y-0.3))
##D   }
##D )
##D # Plot trajectories from a point
##D result <- plot.morse.smale.trajectories.from.point(0.5, 0.5, mixture, grid)
## End(Not run)



### Name: plot.MSD
### Title: Plot Results of Mean Shift Data Denoising
### Aliases: plot.MSD

### ** Examples

## Not run: 
##D # Generate sample data
##D X <- matrix(rnorm(200), ncol = 2)
##D 
##D # Run mean shift smoothing
##D res <- meanshift.data.smoother(X, k = 5)
##D 
##D # Plot original vs denoised data
##D plot(res, type = "dX")
##D 
##D # Plot trajectory at step 3
##D plot(res, type = "dXi", i = 3)
##D 
##D # Plot median k-distances
##D plot(res, type = "kdists")
## End(Not run)




### Name: plot.nerve_cx_spectral_filter
### Title: Plot Nerve Complex Spectral Filter Results
### Aliases: plot.nerve_cx_spectral_filter

### ** Examples

## Not run: 
##D # Create example data and apply filtering
##D coords <- matrix(runif(200), ncol = 2)
##D complex <- create.nerve.complex(coords, k = 5, max.dim = 2)
##D y <- sin(coords[,1] * 2*pi) * cos(coords[,2] * 2*pi) + rnorm(100, 0, 0.1)
##D complex <- set.complex.function.values(complex, y)
##D result <- nerve.cx.spectral.filter(complex, y)
##D 
##D # Create different plot types
##D par(mfrow = c(1, 3))
##D plot(result, type = "parameters")
##D plot(result, type = "predictions", y = y)
##D plot(result, type = "residuals", y = y)
## End(Not run)




### Name: plot.pgmalo
### Title: Plot Method for pgmalo Objects
### Aliases: plot.pgmalo

### ** Examples

# See examples in pgmalo()




### Name: plot.prediction.errors
### Title: Plot Prediction Errors for Multiple Smoothing Models
### Aliases: plot.prediction.errors

### ** Examples

## Not run: 
##D # Using a named list
##D errors_list <- list(
##D   model1 = c(1,2,3,4,5),
##D   model2 = c(2,3,4,5,6)
##D )
##D xvals <- c(1,2,3,4,5)
##D 
##D # Create prediction.errors object
##D pe <- prediction.errors(errors_list, xvals)
##D 
##D # Basic plot
##D plot(pe)
##D 
##D # Using a data frame
##D errors_df <- data.frame(
##D   model1 = c(1,2,3,4,5),
##D   model2 = c(2,3,4,5,6)
##D )
##D 
##D # Create prediction.errors object from data frame
##D pe_df <- prediction.errors(errors_df, xvals = 1:5)
##D 
##D # Plot with custom parameters
##D plot(pe_df,
##D      cols = c("red", "blue"),
##D      ltys = c(1, 2),
##D      main = "Prediction Errors",
##D      legend.pos = "topright")
## End(Not run)




### Name: plot.pwlm
### Title: Plot a Piecewise Linear Model
### Aliases: plot.pwlm

### ** Examples

# Generate sample data
set.seed(123)
x <- 1:20
y <- c(1:10, 20:11) + rnorm(20, 0, 2)

# Fit piecewise linear model with one breakpoint
pwlm_fit <- fit.pwlm(x, y)

# Plot the model
## Not run: 
##D   plot(pwlm_fit)
##D 
##D   # Customize plot appearance
##D   plot(pwlm_fit, main = "My Custom PWLM Plot",
##D        point_color = "blue", line_color = "red")
## End(Not run)




### Name: plot.star_object
### Title: Visualize 2D Star Graph with Smooth Function Values
### Aliases: plot.star_object

### ** Examples

# Generate basic star graph data
data <- generate.star.dataset(
  n.points = 20,
  n.arms = 5,
  min.arm.length = 1,
  max.arm.length = 3,
  dim = 2,
  fn = "exp",
  noise = "norm",
  noise.sd = 0.1
)

# Basic visualization
plot(data)

# Visualization with custom styling and path highlighting
path_data <- list(
  list(vertices = c(1, 2, 3), ref_vertex = 1)
)
plot(data,
  point.size = 2,
  color.edges = TRUE,
  show.vertex.labs = TRUE,
  gpd.obj = path_data
)




### Name: plot.ugkmm
### Title: Plot ugkmm Objects
### Aliases: plot.ugkmm

### ** Examples

# See examples in univariate.gkmm




### Name: plot.upgmalo
### Title: Plot Method for upgmalo Objects
### Aliases: plot.upgmalo

### ** Examples

## Not run: 
##D # Generate example data
##D n <- 100
##D x <- seq(0, 2*pi, length.out = n)
##D y <- sin(x) + rnorm(n, 0, 0.2)
##D 
##D # Fit model
##D fit <- upgmalo(x, y, n.bb = 100)
##D 
##D # Create various plots
##D plot(fit, type = "fit", with.pts = TRUE)
##D plot(fit, type = "diagnostic")
##D plot(fit, type = "residuals")
## End(Not run)




### Name: plot2D.cltr.profiles
### Title: Plot Cluster Profile Lines
### Aliases: plot2D.cltr.profiles

### ** Examples

## Not run: 
##D X <- matrix(rnorm(100), ncol = 10)
##D cltr <- sample(1:3, 10, replace = TRUE)
##D plot2D.cltr.profiles(X, cltr, id = 2, xlab = "Variables", ylab = "Values")
## End(Not run)




### Name: plot2D.colored.graph
### Title: Plot a Graph with Colored Vertices
### Aliases: plot2D.colored.graph

### ** Examples

## Not run: 
##D # Basic usage with simulated data
##D n <- 100  # number of vertices
##D embedding <- matrix(rnorm(2*n), ncol=2)
##D adj.list <- lapply(1:n, function(i) sample(1:n, 5))
##D vertex.colors <- rnorm(n)
##D plot.colored.graph(embedding, adj.list, vertex.colors)
##D 
##D # Custom styling
##D plot2D.colored.graph(
##D     embedding,
##D     adj.list,
##D     vertex.colors,
##D     vertex.size = 1.5,
##D     edge.alpha = 0.3,
##D     color.palette = colorRampPalette(c("navy", "white", "darkred"))(100),
##D     main = "My Graph",
##D     add.legend = TRUE
##D )
## End(Not run)




### Name: plot2D.cont
### Title: Create 2D Plot with Continuous Color Coding
### Aliases: plot2D.cont

### ** Examples

## Not run: 
##D # Generate sample data
##D X <- matrix(rnorm(200), ncol = 2)
##D y <- runif(100)
##D 
##D # Basic plot
##D plot2D.cont(X, y)
##D 
##D # Plot with custom legend position and no axes
##D plot2D.cont(X, y, legend.title = "Values", axes = FALSE)
## End(Not run)




### Name: plot2D.node.level.props
### Title: Plot Node Level Proportions on Graph
### Aliases: plot2D.node.level.props

### ** Examples

## Not run: 
##D # Create example adjacency matrix
##D adj.mat <- matrix(0, 5, 5)
##D adj.mat\code{[1,2]} <- adj.mat\code{[2,1]} <- 1
##D adj.mat\code{[2,3]} <- adj.mat\code{[3,2]} <- 1
##D 
##D # Create proportions matrix
##D props <- matrix(runif(15), ncol = 3)
##D props <- props / rowSums(props)
##D colnames(props) <- c("Type A", "Type B", "Type C")
##D 
##D plot2D.node.level.props(adj.mat, props)
## End(Not run)




### Name: plot2D.tree
### Title: Add Minimal Spanning Tree to 2D Plot
### Aliases: plot2D.tree

### ** Examples

## Not run: 
##D X <- matrix(rnorm(20), ncol = 2)
##D T.edges <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol = 2, byrow = TRUE)
##D 
##D plot(X, pch = 19)
##D plot2D.tree(X, T.edges)
## End(Not run)




### Name: plot3D.cl
### Title: Plot a Specific Cluster with Custom Colors
### Aliases: plot3D.cl

### ** Examples

## Not run: 
##D X <- matrix(rnorm(300), ncol = 3)
##D cltr <- sample(1:5, 100, replace = TRUE)
##D 
##D # Highlight single cluster
##D plot3D.cl(2, cltr, X)
##D 
##D # Highlight multiple clusters
##D plot3D.cl(c(2, 4), cltr, X)
## End(Not run)




### Name: plot3D.cltrs
### Title: Create 3D Plot with Cluster Visualization
### Aliases: plot3D.cltrs

### ** Examples

## Not run: 
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
##D   set.seed(123)
##D   X <- matrix(rnorm(300), ncol = 3)
##D   cltr <- sample(c("A", "B", "C"), 100, replace = TRUE)
##D   plot3D.cltrs(X, cltr, show.cltr.labels = TRUE)
##D }
## End(Not run)



### Name: plot3D.diskEmbdg
### Title: Plot 3d Disk Embedding Object
### Aliases: plot3D.diskEmbdg

### ** Examples

## Not run: 
##D # Assuming ebdg.obj is a properly formatted disk embedding object
##D plot3D.diskEmbdg(ebdg.obj, col = "blue", vertex.col = "red")
## End(Not run)




### Name: plot3D.gaussian_mixture
### Title: 3D Plot of Gaussian Mixture
### Aliases: plot3D.gaussian_mixture

### ** Examples

## Not run: 
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
##D   gm <- get.gaussian.mixture.2d(n.components = 2)
##D   plot3D(gm)
##D }
## End(Not run)



### Name: plot3D.geodesic
### Title: Plot Geodesic Path in 3D
### Aliases: plot3D.geodesic

### ** Examples

## Not run: 
##D library(igraph)
##D # Create a simple graph
##D G <- make_ring(10)
##D S.3d <- matrix(rnorm(30), ncol = 3)
##D 
##D plot3D.plain(S.3d)
##D plot3D.geodesic(S.3d, 1, 5, G, S.3d)
## End(Not run)




### Name: plot3D.graph
### Title: 3D Graph Visualization with Height Values
### Aliases: plot3D.graph

### ** Examples

## Not run: 
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
##D   adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2), c(2))
##D   weight.list <- list(c(1, 1), c(1, 1, 1), c(1, 1), c(1))
##D   z.values <- c(0.5, 1.2, 0.8, 1.5)
##D   g <- ggraph(adj.list, weight.list)
##D   plot3D.graph(g, z.values)
##D }
## End(Not run)



### Name: plot3D.path
### Title: Plot Path in 3D Graph
### Aliases: plot3D.path

### ** Examples

## Not run: 
##D V <- matrix(rnorm(15), ncol = 3)
##D s <- c(1, 3, 2, 5, 4)
##D 
##D plot3D.plain(V)
##D plot3D.path(s, V, edge.col = "red")
## End(Not run)




### Name: plot3D.tree
### Title: Add Minimal Spanning Tree to 3D Plot
### Aliases: plot3D.tree

### ** Examples

## Not run: 
##D X <- matrix(rnorm(30), ncol = 3)
##D # Assuming T.edges is computed from a minimal spanning tree algorithm
##D T.edges <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol = 2, byrow = TRUE)
##D 
##D plot3D.plain(X)
##D plot3D.tree(X, T.edges)
## End(Not run)




### Name: plotlcor.1D
### Title: Plot Output from lcor.1D()
### Aliases: plotlcor.1D

### ** Examples

## Not run: 
##D # Assuming 'result' is output from lcor.1D()
##D plotlcor.1D(result, with.CrI = TRUE, title = "Local Correlations")
##D 
##D # Plot with credible intervals as lines instead of polygon
##D plotlcor.1D(result, with.CrI = TRUE, CrI.as.polygon = FALSE)
## End(Not run)




### Name: predict.maelog
### Title: Predict Method for maelog Objects
### Aliases: predict.maelog

### ** Examples

## Not run: 
##D # Fit model
##D set.seed(123)
##D x <- seq(0, 1, length.out = 100)
##D y <- rbinom(100, 1, plogis(10 * (x - 0.5)))
##D fit <- maelog(x, y, with.errors = TRUE)
##D 
##D # Predictions at new points
##D x.new <- seq(0.2, 0.8, by = 0.1)
##D pred <- predict(fit, x.new)
##D 
##D # Predictions with standard errors
##D pred.se <- predict(fit, x.new, se.fit = TRUE)
## End(Not run)



### Name: predict.magelo
### Title: Predicting values from a magelo model
### Aliases: predict.magelo

### ** Examples

## Not run: 
##D model <- magelo(x, y)
##D predictions <- predict(model, newdata = c(1, 2, 3))
## End(Not run)



### Name: predict.magelog
### Title: Predict Method for magelog Objects
### Aliases: predict.magelog

### ** Examples

## Not run: 
##D # Fit model
##D fit <- magelog(x, y)
##D 
##D # Predictions at original points
##D pred <- predict(fit)
##D 
##D # Predictions at new points
##D x_new <- seq(min(x), max(x), length.out = 50)
##D pred_new <- predict(fit, newdata = x_new)
## End(Not run)




### Name: print.basin_cx
### Title: Print Method for Basin Complex Objects
### Aliases: print.basin_cx

### ** Examples

## Not run: 
##D # Create a basin complex
##D adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
##D weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
##D y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
##D basin_cx <- create.basin.cx(adj_list, weight_list, y)
##D 
##D # Print basic information
##D print(basin_cx)
##D # or simply:
##D basin_cx
## End(Not run)




### Name: print.circular_parameterization
### Title: Print Method for circular_parameterization Objects
### Aliases: print.circular_parameterization

### ** Examples

# Create example graph
n <- 6
adj.list <- lapply(seq_len(n), function(i) {
  c(if (i == n) 1 else i + 1, if (i == 1) n else i - 1)
})
weight.list <- lapply(adj.list, function(adj) rep(1.0, length(adj)))
result <- parameterize.circular.graph(adj.list, weight.list)
print(result)




### Name: print.gflow_cx
### Title: Print Method for gflow_cx Objects
### Aliases: print.gflow_cx

### ** Examples

## Not run: 
##D # Create simple example
##D adj.list <- list(c(2,3), c(1,3), c(1,2))
##D weight.list <- list(c(1,1), c(1,1), c(1,1))
##D y <- c(0.5, 1.0, 0.7)
##D 
##D result <- create.gflow.cx(adj.list, weight.list, y, verbose = FALSE)
##D print(result)
## End(Not run)



### Name: print.mknn_graphs
### Title: Print Method for mknn_graphs Objects
### Aliases: print.mknn_graphs

### ** Examples

## Not run: 
##D # Create multiple graphs
##D X <- matrix(rnorm(100 * 3), ncol = 3)
##D graphs <- create.mknn.graphs(X, kmin = 5, kmax = 10)
##D print(graphs)
## End(Not run)



### Name: print.mst_completion_graph
### Title: Print Method for MST Completion Graph Objects
### Aliases: print.mst_completion_graph

### ** Examples

## Not run: 
##D # Create and print a graph
##D X <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
##D graph <- create.cmst.graph(X, q.thld = 0.8)
##D print(graph)
## End(Not run)




### Name: print.path.graph
### Title: Print Method for path.graph Objects
### Aliases: print.path.graph

### ** Examples

## Not run: 
##D pg <- create.path.graph(graph, edge.lengths, h = 2)
##D print(pg)
## End(Not run)




### Name: print.summary.basin_cx
### Title: Print Method for Basin Complex Summary Objects
### Aliases: print.summary.basin_cx

### ** Examples

## Not run: 
##D # Create a basin complex
##D adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
##D weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
##D y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
##D basin_cx <- create.basin.cx(adj_list, weight_list, y)
##D 
##D # Get and print detailed summary
##D basin_summary <- summary(basin_cx)
##D print(basin_summary)
##D # or simply:
##D basin_summary
## End(Not run)




### Name: print.summary.ggflow
### Title: Print Gradient Flow Structure Summary
### Aliases: print.summary.ggflow

### ** Examples

## Not run: 
##D flow <- construct.graph.gradient.flow(adj.list, weight.list, y, scale)
##D summary(flow)
## End(Not run)




### Name: print.summary.mst_completion_graph
### Title: Print Summary of MST Completion Graph
### Aliases: print.summary.mst_completion_graph

### ** Examples

## Not run: 
##D X <- matrix(rnorm(75 * 4), nrow = 75, ncol = 4)
##D graph <- create.cmst.graph(X)
##D graph_summary <- summary(graph)
##D print(graph_summary, digits = 4)
## End(Not run)




### Name: print.summary.pgmalo
### Title: Print Method for summary.pgmalo Objects
### Aliases: print.summary.pgmalo

### ** Examples

## Not run: 
##D fit <- pgmalo(neighbors, edge_lengths, y)
##D summary(fit)
##D 
##D # Show more CV details
##D print(summary(fit), show_cv_details = TRUE)
## End(Not run)




### Name: print.summary.upgmalo
### Title: Print Method for summary.upgmalo Objects
### Aliases: print.summary.upgmalo

### ** Examples

# See examples in upgmalo() and summary.upgmalo()




### Name: prof.fn
### Title: Extract Most Abundant ASVs with Taxonomy Information
### Aliases: prof.fn

### ** Examples

# Create example data
S <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(S) <- paste0("Sample", 1:10)
colnames(S) <- paste0("ASV", 1:10)

# Extract profile for first sample
prof <- prof.fn("Sample1", S, n.prof = 3)

# With taxonomy
taxonomy <- paste0("Taxon", 1:10)
names(taxonomy) <- colnames(S)
prof_with_tax <- prof.fn(1, S, bm.tx = taxonomy, n.prof = 3, verbose = TRUE)




### Name: project.onto.subspace
### Title: Projects a Vector onto a Subspace in R^n
### Aliases: project.onto.subspace

### ** Examples

u1 <- c(1, 0, 0)
u2 <- c(0, 1, 0)
u3 <- c(0, 0, 1)
U <- matrix(c(u1, u2, u3), nrow=3)
x <- c(2, 3, 4)
project.onto.subspace(x, U)




### Name: qexp.gaussian
### Title: q-Exponential Gaussian Function
### Aliases: qexp.gaussian

### ** Examples

## Not run: 
##D x <- seq(-5, 5, length.out = 100)
##D plot(x, qgaussian(x, q = 4), type = "l")
## End(Not run)



### Name: qlaplace
### Title: Laplace Distribution Quantile Function
### Aliases: qlaplace

### ** Examples

p <- seq(0, 1, by = 0.1)
q <- qlaplace(p, location = 0, scale = 1)
plot(p, q, type = "l", main = "Laplace Quantile Function")




### Name: quantize.for.legend
### Title: Quantize Continuous Variable for Legend Display
### Aliases: quantize.for.legend

### ** Examples

## Not run: 
##D y <- rnorm(100)
##D q <- quantize.for.legend(y, n.levels = 5)
##D plot(1:100, y, col = q$y.cols, pch = 19)
##D legend("topright", legend = q$legend.labs, fill = q$y.col.tbl)
## End(Not run)




### Name: radial.Lp.transform
### Title: Perform Lp Radial Transformation
### Aliases: radial.Lp.transform

### ** Examples

X <- matrix(runif(20), ncol = 2)
transformed_X <- radial.Lp.transform(X, p = 0.5)




### Name: Rdensity.distance
### Title: Calculates a Symmetric Density-Associated Distance
### Aliases: Rdensity.distance

### ** Examples

density.func <- function(x) { return(dnorm(x, mean=0, sd=1)) }
p <- c(0, 0)
q <- c(1, 1)
p.density <- density.func(p)
q.density <- density.func(q)
Rdensity.distance(p, q, p.density, q.density, density.func)




### Name: remove.close.neighbors
### Title: Removes Close Neighbors
### Aliases: remove.close.neighbors

### ** Examples

## Not run: 
##D x <- c(0.1, 0.2, 0.4, 0.7, 1.0)
##D min.dist <- 0.25
##D remove.close.neighbors(x, min.dist)
## End(Not run)



### Name: remove.knn.outliers
### Title: Remove Outliers from a State Space Using k-Nearest Neighbors
### Aliases: remove.knn.outliers

### ** Examples

# Create a sample dataset with outliers
set.seed(123)
n_normal <- 1000
n_outliers <- 10

# Generate normal data
normal_data <- matrix(rnorm(n_normal * 2), ncol = 2)

# Generate outliers far from the main cluster
outliers <- cbind(rnorm(n_outliers, mean = 10, sd = 0.5),
                  rnorm(n_outliers, mean = 10, sd = 0.5))

# Combine data
S <- rbind(normal_data, outliers)
y <- c(rnorm(n_normal), rnorm(n_outliers, mean = 5))

# Remove outliers using the default method
result <- remove.knn.outliers(S, y, K = 10)

# Print summary
cat("Number of outliers detected:", result$n.outliers, "\n")
cat("Percentage of data retained:",
    round(100 * nrow(result$S.q) / nrow(S), 2), "%\n")

# Use a different method with more conservative threshold
result2 <- remove.knn.outliers(S, y, p = 0.95, method = "dist.factor",
                               dist.factor = 5, K = 10)

## Not run: 
##D # Plot the original and filtered state space
##D par(mfrow = c(1, 2))
##D plot(S, col = "gray", main = "Original State Space",
##D      xlab = "Dimension 1", ylab = "Dimension 2")
##D points(S[!result$idx, ], col = "red", pch = 16, cex = 1.2)
##D 
##D plot(result$S.q, col = "blue", pch = 16,
##D      main = "Filtered State Space",
##D      xlab = "Dimension 1", ylab = "Dimension 2")
##D legend("topright", legend = c("Retained", "Removed"),
##D        col = c("blue", "red"), pch = 16)
## End(Not run)




### Name: replace.basin.label
### Title: Replace a Basin Label in Gradient Flow Object
### Aliases: replace.basin.label

### ** Examples

## Not run: 
##D # Replace a generic label with a meaningful one
##D flow_updated <- replace.basin.label(flow, "M1", "HighState")
## End(Not run)




### Name: restric.to.Linf.unit.sphere
### Title: Restrict Points to L-infinity Unit Sphere
### Aliases: restric.to.Linf.unit.sphere

### ** Examples

X <- matrix(runif(20, -2, 2), ncol = 2)
restricted_X <- restric.to.Linf.unit.sphere(X)




### Name: right.asymmetric.bump.fn
### Title: Creates Right Asymmetric Bump Function
### Aliases: right.asymmetric.bump.fn

### ** Examples

## Not run: 
##D x <- seq(-2, 2, length.out = 100)
##D plot(x, AsymmetricBumpFunction(x), type = "l")
## End(Not run)



### Name: right.winsorize
### Title: Right Winsorize a numeric vector
### Aliases: right.winsorize

### ** Examples

y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
right.winsorize(y, p = 0.1)




### Name: rlaplace
### Title: Generate Random Variates from the Laplace Distribution
### Aliases: rlaplace

### ** Examples

# Generate 1000 random variates from the standard Laplace distribution
x <- rlaplace(1000)

# Generate 500 random variates with location 2 and scale 3, using a specific seed
y <- rlaplace(500, location = 2, scale = 3, seed = 12345)




### Name: rm.self.loops
### Title: Removes self-loops from a graph
### Aliases: rm.self.loops

### ** Examples

graph <- list(
  "A" = c(1, 2),
  "B" = c(2, 3),
  "C" = c(1, 3),
  "D" = c(4)
)

graph.no.self.loops <- rm.self.loops(graph)
print(graph.no.self.loops)




### Name: rm.SS.outliers
### Title: Remove outliers from a state space
### Aliases: rm.SS.outliers

### ** Examples

## Not run: 
##D # Create example data
##D set.seed(123)
##D S <- matrix(rnorm(200), ncol = 2)
##D y <- rnorm(100)
##D 
##D # Remove outliers
##D result <- rm.SS.outliers(S, y, p = 0.95)
##D 
##D # Check how many points were removed
##D sum(!result$idx)
## End(Not run)




### Name: robust.log.transform
### Title: Returns robust log transform of a non-negative vector
### Aliases: robust.log.transform

### ** Examples

# Example 1: Basic usage
x <- c(0, 1, 10, 100, 1000)
robust.log.transform(x)

# Example 2: Scale invariance property
x1 <- c(0, 1, 2, 5, 10)
x2 <- c(0, 10, 20, 50, 100)  # x1 * 10
# The relative differences are preserved
robust.log.transform(x1)
robust.log.transform(x2)



### Name: robust.transform
### Title: Returns robust transform of a non-negative vector
### Aliases: robust.transform

### ** Examples

# Example 1: Basic usage
x <- c(0, 1, 4, 16, 64)
robust.transform(x)
# Geometric mean of positive transformed values equals 1

# Example 2: Comparison with arithmetic mean normalization
x <- c(0, 1, 2, 100)  # Outlier present
# Robust (geometric mean) normalization
robust.transform(x)
# Compare to arithmetic mean normalization
x_arith <- x
x_arith[x > 0] <- x[x > 0] / mean(x[x > 0])
x_arith  # More affected by the outlier



### Name: robust.z.normalize
### Title: Perform Robust Z-Score Normalization
### Aliases: robust.z.normalize

### ** Examples

# Basic usage
data <- c(1, 2, 3, 100, 4, 5, 6)  # Note the outlier
robust.z.normalize(data)

# Handling missing values
data.with.na <- c(1, NA, 3, 100, 4, NA, 6)
robust.z.normalize(data.with.na)




### Name: robust.zscore
### Title: Robust Z-score normalization using median and MAD
### Aliases: robust.zscore

### ** Examples

# Example with random data
set.seed(123)
example.data <- matrix(rnorm(100, 5, 2), ncol=5)
# Add some outliers
example.data[1,1] <- 25
example.data[2,3] <- -15
normalized.data <- robust.zscore(example.data)



### Name: runif.sphere
### Title: Uniformly Sample Points from a Unit Sphere Surface
### Aliases: runif.sphere

### ** Examples

# Generate 10 points on a 2D circle (1-sphere)
points_circle <- runif.sphere(10, 1)

# Generate 100 points on a 3D sphere surface (2-sphere)
points_sphere <- runif.sphere(100, 2)




### Name: runif.torus
### Title: Generate Uniform Random Sample from n-dimensional Torus
### Aliases: runif.torus

### ** Examples

# Generate 100 points on a circle (1-dimensional torus)
circle_points <- runif.torus(100, dim = 1)

# Generate 100 points on a 2-dimensional torus
torus_points <- runif.torus(100, dim = 2)




### Name: S.version.graph.kernel.smoother
### Title: Graph-Based Kernel Smoother with Buffer Zone Cross-Validation
### Aliases: S.version.graph.kernel.smoother

### ** Examples

## Not run: 
##D n.pts <- 100
##D gm <- generate.1d.gaussian.mixture(
##D     n.points = n.pts,
##D     x.knot = c(0, 10),
##D     y.knot = c(10, 2.5),
##D     sd.knot = 1.5,
##D     x.offset = 3)
##D x <- sort(runif(n.pts, min = min(gm$x), max = max(gm$x)))
##D x.graph <- create.bi.kNN.chain.graph(k = 1, x = x, y = gm$y)
##D 
##D y.smooth <- approx(gm$x, gm$y, xout = x)$y
##D sigma <- 1
##D eps <- rnorm(n.pts, 0, sigma)
##D y <- y.smooth + eps
##D 
##D g <- ggraph(x.graph$adj.list, x.graph$edge.lengths)
##D plot(g, y.smooth)
##D 
##D plot(gm$x, gm$y, type = "l", las = 1, ylim = range(y), col = "red", xlab = "x", ylab = "y")
##D points(x, y)
##D legend("topright", legend = c("y.smooth", "y"), lty = c(1,NA), pch = c(NA,1),
##D col = c("red", "black"), inset = 0.1)
##D 
##D gks.res <- graph.kernel.smoother(x.graph$adj.list,
##D                                  x.graph$edge.lengths,
##D                                  y,
##D                                  min.bw.factor = 0.025,
##D                                  max.bw.factor = 0.5,
##D                                  n.bws = 20,
##D                                  log.grid = TRUE,
##D                                  vertex.hbhd.min.size = 3,
##D                                  dist.normalization.factor = 1.1,
##D                                  use.uniform.weights = FALSE,
##D                                  buffer.hops = 1,
##D                                  auto.buffer.hops = FALSE,
##D                                  kernel.type = 7L,
##D                                  n.folds = 10,
##D                                  with.bw.predictions = TRUE,
##D                                  verbose = TRUE)
##D 
##D # View results
##D print(gks.res)
##D summary(gks.res)
##D 
##D # Computing Mean Absolute and Mean Squared Errors
##D mae <- c()
##D for (i in seq(ncol(gks.res$bw_predictions))) {
##D     mae[i] <- mean(abs(y.smooth - gks.res$bw_predictions[,i]))
##D }
##D plot(mae, las = 1, type = 'b', xlab = "Bandwidth indices", ylab = "Mean Absolute Error")
##D which.min(mae)
##D gks.res$opt_bw_idx
##D 
##D plot(gks.res$bw_mean_abs_errors, las = 1, type = "b")
##D abline(v = gks.res$opt_bw_idx, lty = 2)
##D which.min(gks.res$bw_mean_abs_errors)
##D 
##D plot(gm$x, gm$y, type = "l", las = 1, ylim = range(y), col = "red")
##D points(x, y)
##D lines(x, gks.res$predictions, col = "blue")
## End(Not run)




### Name: sample.from.empirical.distribution
### Title: Sample points from empirical distribution with optional linear
###   approximation
### Aliases: sample.from.empirical.distribution

### ** Examples

# Generate bimodal test data
y <- c(rnorm(1000, 0.3, 0.1), rnorm(1000, 0.7, 0.1))

# Sample using step function (default)
samples1 <- sample.from.empirical.distribution(1000, y)

# Sample using linear approximation
samples2 <- sample.from.empirical.distribution(1000, y, use.linear.approximation=TRUE)

# Compare distributions
par(mfrow=c(1,3))
hist(y, main="Original Data")
hist(samples1, main="Step Function Sampling")
hist(samples2, main="Linear Interpolation Sampling")




### Name: scaled.log.tan
### Title: Applies a Scaled Log Transformation
### Aliases: scaled.log.tan

### ** Examples

## Not run: 
##D x <- 0.1
##D scaled.log.tan(x)
## End(Not run)



### Name: select3D.points.profiles
### Title: Select 3D Points and Show Profiles
### Aliases: select3D.points.profiles

### ** Examples

## Not run: 
##D if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
##D   set.seed(123)
##D   n <- 100
##D   X <- matrix(rnorm(n * 3), ncol = 3)
##D   Z <- matrix(rpois(n * 10, lambda = 5), ncol = 10)
##D   colnames(Z) <- paste0("Feature", 1:10)
##D 
##D   rgl::open3d()
##D   rgl::plot3d(X, col = "blue")
##D 
##D   sel <- select3D.points.profiles(X, Z, show.profiles = TRUE)
##D   rgl::points3d(X[sel$idx, , drop = FALSE], col = "red", size = 10)
##D }
## End(Not run)



### Name: select3D.points
### Title: Select Points in 3D Space
### Aliases: select3D.points

### ** Examples

## Not run: 
##D if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
##D   set.seed(1)
##D   X <- matrix(rnorm(300), ncol = 3)
##D   rgl::open3d()
##D   rgl::plot3d(X)
##D   # Return coordinates of selected points:
##D   coords <- select3D.points()
##D   # Or return (id, index) pairs:
##D   # idxmat <- select3D.points(value = FALSE)
##D }
## End(Not run)



### Name: separatrices.with.cells.plot
### Title: Plot Separatrices with Cells
### Aliases: separatrices.with.cells.plot

### ** Examples

## Not run: 
##D # Typically used with pre-computed separatrices
##D # Example with a simple two-maximum system
##D separatrices <- list(
##D   boundary_left = c(0, 0.5),
##D   boundary_right = c(1, 0.5),
##D   boundary_top = c(0.5, 1),
##D   boundary_bottom = c(0.5, 0),
##D   saddle_descending1 = matrix(c(0.5, 0.5, 0.3, 0.3), nrow = 2, byrow = TRUE),
##D   saddle_descending2 = matrix(c(0.5, 0.5, 0.7, 0.7), nrow = 2, byrow = TRUE)
##D )
##D separatrices.with.cells.plot(separatrices,
##D                             M1.pos = c(0.3, 0.7),
##D                             M2.pos = c(0.7, 0.3),
##D                             saddle_pos = c(0.5, 0.5))
## End(Not run)



### Name: set.boundary
### Title: Compute Boundary Vertices of a Subset in a Graph
### Aliases: set.boundary

### ** Examples

# Create a simple graph adjacency list
adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2), c(2, 5), c(4))

# Define a subset of vertices
U <- c(1, 2, 3)

# Compute the boundary
boundary <- set.boundary(U, adj.list)
print(boundary)  # Returns c(2) since only vertex 2 is adjacent to vertex 4 (outside U)




### Name: shifted.tan
### Title: Applies a Linear Transformation and Tangent Function
### Aliases: shifted.tan

### ** Examples

## Not run: 
##D x <- 0.001
##D shifted.tan(x)
## End(Not run)



### Name: shortest.path
### Title: Computes Shortest Path Distances for a Selected Set of Vertices
###   of a Graph
### Aliases: shortest.path

### ** Examples

graph <- list(c(2,3), c(1,3,4), c(1,2,4), c(2,3))
edge.lengths <- list(c(1,4), c(1,2,5), c(4,2,1), c(5,1))
vertices <- c(1,3,4)
result <- shortest.path(graph, edge.lengths, vertices)
print(result)




### Name: show.3d.cl
### Title: Highlight a Specific Cluster in 3D Space
### Aliases: show.3d.cl

### ** Examples

## Not run: 
##D X <- matrix(rnorm(300), ncol = 3)
##D cltr <- sample(1:5, 100, replace = TRUE)
##D 
##D # Highlight cluster 3
##D show.3d.cl(3, cltr, X)
##D 
##D # Show with labels
##D rownames(X) <- paste0("Point", 1:100)
##D show.3d.cl(3, cltr, X, show.labels = TRUE)
## End(Not run)




### Name: show.cltrs
### Title: Display Cluster Labels in 3D Space (Headless/CRAN-safe)
### Aliases: show.cltrs

### ** Examples

## Not run: 
##D # Minimal example; runs even without rgl (no plotting).
##D set.seed(1)
##D X <- matrix(rnorm(300), ncol = 3)
##D cl <- sample(1:5, nrow(X), replace = TRUE)
##D centers <- show.cltrs(X, cl, show.plot = FALSE)
##D 
##D # If rgl is available, render off-screen safely:
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
##D   # Suppose you've already drawn points via your helper:
##D   # plot3D.plain(X)
##D   show.cltrs(X, cl, cex = 1.1, adj = c(0.5, 1), show.plot = TRUE)
##D }
## End(Not run)




### Name: show.profile
### Title: Display Profile of Top Components
### Aliases: show.profile

### ** Examples

# Create example data
set.seed(123)
Z <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 10, ncol = 10)
colnames(Z) <- paste0("Feature", 1:10)

# Get profile for row 1
profile <- show.profile(1, Z, n.comp = 3)
print(profile)




### Name: skewed.gaussian
### Title: Compute Value of a Single Skewed Gaussian
### Aliases: skewed.gaussian

### ** Examples

skewed.gaussian(0, mu = 0, sigma = 1, alpha = 2)




### Name: spectral.lowess.graph.smoothing
### Title: Iterative Spectral LOWESS Graph Smoothing
### Aliases: spectral.lowess.graph.smoothing

### ** Examples

## Not run: 
##D # Create a graph and data matrix
##D graph <- create.iknn.graph(X, k = 10, pruning.thld = 0.1)
##D 
##D # Apply spectral LOWESS graph smoothing - traditional approach
##D result1 <- spectral.lowess.graph.smoothing(
##D   adj.list = graph$adj_list,
##D   weight.list = graph$weight_list,
##D   X = X,
##D   max.iterations = 10,
##D   convergence.threshold = 1e-4,
##D   convergence.type = 1,  # MAX.ABSOLUTE.DIFF
##D   k = 10,
##D   pruning.thld = 0.1,
##D   n.evectors = 8,
##D   verbose = TRUE
##D )
##D 
##D # Apply spectral LOWESS graph smoothing with boosting
##D result2 <- spectral.lowess.graph.smoothing(
##D   adj.list = graph$adj_list,
##D   weight.list = graph$weight_list,
##D   X = X,
##D   max.iterations = 10,
##D   convergence.threshold = 1e-4,
##D   convergence.type = 1,  # MAX.ABSOLUTE.DIFF
##D   k = 10,
##D   pruning.thld = 0.1,
##D   n.evectors = 8,
##D   switch.to.residuals.after = 2,  # Switch to boosting after 2 iterations
##D   verbose = TRUE
##D )
##D 
##D # Access final smoothed data matrix
##D X.smoothed <- result2$smoothed.X[[length(result2$smoothed.X)]]
##D 
##D # Plot convergence metrics for both approaches
##D plot(result1$convergence.metrics, type = "b", col = "blue",
##D      xlab = "Iteration", ylab = "Convergence Metric")
##D lines(result2$convergence.metrics, type = "b", col = "red")
##D legend("topright", legend = c("Traditional", "Boosting"), 
##D        col = c("blue", "red"), lty = 1)
## End(Not run)




### Name: standardize.string
### Title: Standardize String for File Names
### Aliases: standardize.string

### ** Examples

# Standardize a metabolite name
standardize.string("Glucose (6-phosphate)")
# Returns: "Glucose_6-phosphate"

# Standardize a complex name
standardize.string("Compound A/B + C's metabolite*")
# Returns: "Compound_A_B__Cs_metabolite"

# Use for creating file names
sample_name <- "Sample 1: Day 3 (replicate)"
file_name <- paste0(standardize.string(sample_name), ".csv")
# Results in: "Sample_1_Day_3_replicate.csv"




### Name: subdivide.path
### Title: Subdivides a path of points in R^n into a uniform grid of points
###   along the path
### Aliases: subdivide.path

### ** Examples

# Create a simple 2D path
path <- matrix(c(0,0, 1,0, 2,1, 3,1), ncol=2, byrow=TRUE)
# Subdivide into 10 evenly spaced points
subdivided <- subdivide.path(path, 10)




### Name: summarize.and.test.cst
### Title: Summarize and Test CST Patterns with Binary Outcomes
### Aliases: summarize.and.test.cst

### ** Examples

## Not run: 
##D # Simulated longitudinal CST data
##D set.seed(123)
##D n_subjects <- 50
##D n_timepoints <- 5
##D 
##D # Generate subject IDs and CSTs
##D subj_ids <- rep(1:n_subjects, each = n_timepoints)
##D csts <- sample(c("I", "II", "III", "IV"), n_subjects * n_timepoints,
##D                replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))
##D 
##D # Binary outcome (constant per subject)
##D outcomes <- rep(rbinom(n_subjects, 1, 0.3), each = n_timepoints)
##D 
##D # Analyze using most frequent CST
##D result_mode <- summarize.and.test.cst(
##D   x = csts,
##D   y = outcomes,
##D   subj.ids = subj_ids,
##D   summary.method = "most_frequent"
##D )
##D 
##D # Analyze using CST proportions
##D result_prop <- summarize.and.test.cst(
##D   x = csts,
##D   y = outcomes,
##D   subj.ids = subj_ids,
##D   summary.method = "proportions"
##D )
##D 
##D print(result_prop$summary)
## End(Not run)




### Name: summary.basin_cx
### Title: Summary Method for Basin Complex Objects
### Aliases: summary.basin_cx

### ** Examples

## Not run: 
##D # Create a basin complex
##D adj_list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2,5), c(3,4))
##D weight_list <- list(c(1,2), c(1,1,3), c(2,1,2), c(3,1), c(2,1))
##D y <- c(2.5, 1.8, 3.2, 0.7, 2.1)
##D basin_cx <- create.basin.cx(adj_list, weight_list, y)
##D 
##D # Get detailed summary
##D basin_summary <- summary(basin_cx)
##D 
##D # Access specific components
##D basin_summary$minima_info
##D basin_summary$ascending_basins_jaccard_index
## End(Not run)




### Name: summary.geodesic_stats
### Title: Summary method for geodesic_stats objects
### Aliases: summary.geodesic_stats

### ** Examples

## Not run: 
##D stats <- compute.geodesic.stats(adj.list, weight.list)
##D summary(stats)  # Default type = "rays"
##D summary(stats, type = "composite")
##D summary(stats, type = "overlap")
## End(Not run)




### Name: summary.ggflow
### Title: Summarize Gradient Flow Structure
### Aliases: summary.ggflow

### ** Examples

## Not run: 
##D flow <- construct.graph.gradient.flow(adj.list, weight.list, y, scale)
##D summary(flow)
## End(Not run)




### Name: summary.iknn_graphs
### Title: Summarize an iknn_graphs Object
### Aliases: summary.iknn_graphs

### ** Examples

# Create sample data
set.seed(123)
x <- matrix(rnorm(1000), ncol = 5)

# Generate intersection kNN graphs
iknn.res <- create.iknn.graphs(x, kmin = 3, kmax = 10)

# Summarize the geometrically pruned graphs
summary(iknn.res)

# Summarize the intersection-size pruned graphs
summary(iknn.res, use_isize_pruned = TRUE)




### Name: summary.IkNN
### Title: Summarize IkNN Graph Object
### Aliases: summary.IkNN

### ** Examples

# Generate sample data
set.seed(123)
X <- matrix(rnorm(100 * 5), ncol = 5)
result <- create.single.iknn.graph(X, k = 3)
summary(result)




### Name: summary.local_extrema
### Title: Summarize Local Extrema Detection Results
### Aliases: summary.local_extrema

### ** Examples

# Create example data
adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
y <- c(1, 3, 2, 5, 1)

# Detect and summarize maxima
maxima <- detect.local.extrema(adj.list, weight.list, y, 2, 2)
summary(maxima)




### Name: summary.mknn_graphs
### Title: Summary Method for mknn_graphs Objects
### Aliases: summary.mknn_graphs

### ** Examples

## Not run: 
##D # Create sample data
##D set.seed(123)
##D X <- matrix(rnorm(200 * 5), ncol = 5)
##D 
##D # Generate MkNN graphs
##D mknn_result <- create.mknn.graphs(X, kmin = 5, kmax = 15, compute.full = TRUE)
##D 
##D # Display summary
##D summary(mknn_result)
##D 
##D # Store summary statistics for plotting
##D stats <- summary(mknn_result)
##D plot(stats$k, stats$mean_degree, type = "b",
##D      xlab = "k", ylab = "Mean Degree",
##D      main = "Mean Vertex Degree vs k")
## End(Not run)




### Name: summary.mst_completion_graph
### Title: Summary Method for MST Completion Graph Objects
### Aliases: summary.mst_completion_graph

### ** Examples

## Not run: 
##D X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
##D graph <- create.cmst.graph(X, q.thld = 0.85)
##D summary(graph)
## End(Not run)




### Name: summary.path.graph
### Title: Summary Method for path.graph Objects
### Aliases: summary.path.graph

### ** Examples

## Not run: 
##D pg <- create.path.graph(graph, edge.lengths, h = 3)
##D summary(pg)
## End(Not run)




### Name: summary.pgmalo
### Title: Summary Method for pgmalo Objects
### Aliases: summary.pgmalo

### ** Examples

## Not run: 
##D # Create example data
##D n <- 50
##D neighbors <- vector("list", n)
##D edge_lengths <- vector("list", n)
##D for(i in 1:n) {
##D   if(i == 1) {
##D     neighbors[[i]] <- 2L
##D     edge_lengths[[i]] <- 1.0
##D   } else if(i == n) {
##D     neighbors[[i]] <- (n-1L)
##D     edge_lengths[[i]] <- 1.0
##D   } else {
##D     neighbors[[i]] <- c(i-1L, i+1L)
##D     edge_lengths[[i]] <- c(1.0, 1.0)
##D   }
##D }
##D y <- rnorm(n)
##D 
##D # Fit model and get summary
##D fit <- pgmalo(neighbors, edge_lengths, y)
##D summary(fit)
## End(Not run)




### Name: summary.upgmalo
### Title: Summary Method for upgmalo Objects
### Aliases: summary.upgmalo

### ** Examples

# See examples in upgmalo()




### Name: synthetic.1D.spline
### Title: Generates the values of a synthetic 1d spline function over a
###   uniform grid in the pre-specified range
### Aliases: synthetic.1D.spline

### ** Examples

## Not run: 
##D synth <- synthetic.1D.spline(n.lmax = 5)
##D plot(synth$x, synth$y, type = "l")
##D points(synth$x.lmax, synth$y.lmax, col = "red")
## End(Not run)



### Name: synthetic.mixture.of.gaussians
### Title: Generate a Synthetic Smooth Function in One Dimension
### Aliases: synthetic.mixture.of.gaussians

### ** Examples

## Not run: 
##D # Example usage
##D x <- seq(0, 1, length.out = 10)
##D synthetic.fn <- generate.synthetic.function(x)
##D point <- 0.5 # A point in 1D space
##D value <- synthetic.fn(point) # Evaluate the function at the given point
## End(Not run)



### Name: synthetic.xD.spline
### Title: Generate a Spline Function
### Aliases: synthetic.xD.spline

### ** Examples

## Not run: 
##D # Example usage
##D x <- seq(0, 10, length.out = 10)
##D y <- sin(x)  # Example y values
##D spline.fn <- synthetic.spline(x, y)
##D point <- 5  # An x-coordinate
##D value <- spline.fn(point)  # Evaluate the spline function at this point
## End(Not run)



### Name: thresholded.sigmoid
### Title: Sigmoidal Function with Lower and Upper Thresholds
### Aliases: thresholded.sigmoid

### ** Examples

# Define thresholds
q90 <- 0.1  # Lower threshold (e.g., 90th percentile)
q95 <- 0.15 # Upper threshold (e.g., 95th percentile)

# Apply to a sequence of values
x_seq <- seq(0, 0.3, length.out = 100)
y_values <- thresholded.sigmoid(x_seq, q90, q95)

# Plot the result
plot(x_seq, y_values, type = "l",
     xlab = "Input Value", ylab = "Transformed Value",
     main = "Thresholded Sigmoid Function")
abline(v = c(q90, q95), lty = 2, col = c("blue", "red"))




### Name: total.associations
### Title: Calculate Total Associations from Geodesic Analysis
### Aliases: total.associations

### ** Examples

## Not run: 
##D # Assuming obj is a result from E.geodesic.X()
##D # Calculate associations with distance
##D ta_results <- total.associations(obj)
##D 
##D # Calculate associations with an outcome
##D ED <- rnorm(400)  # Example outcome vector
##D ta_results_ED <- total.associations(obj, ED = ED)
## End(Not run)




### Name: triangle.plot
### Title: triangle.plot
### Aliases: triangle.plot

### ** Examples

## Not run: 
##D   library(igraph)
##D 
##D   # Define the triangle plot function
##D   triangle.plot <- function(coords, v = NULL, params) {
##D       vertex.color <- params("vertex", "color")
##D       if (length(vertex.color) != 1 && !is.null(v)) {
##D           vertex.color <- vertex.color\[v\]
##D       }
##D       vertex.frame.color <- params("vertex", "frame.color")
##D       if (length(vertex.frame.color) != 1 && !is.null(v)) {
##D           vertex.frame.color <- vertex.frame.color\[v\]
##D       }
##D       vertex.frame.width <- params("vertex", "frame.width")
##D       if (length(vertex.frame.width) != 1 && !is.null(v)) {
##D           vertex.frame.width <- vertex.frame.width\[v\]
##D       }
##D       vertex.size <- 1/200 * params("vertex", "size")
##D       if (length(vertex.size) != 1 && !is.null(v)) {
##D           vertex.size <- vertex.size\[v\]
##D       }
##D       vertex.size <- rep(vertex.size, length.out = nrow(coords))
##D 
##D       # Adjust the size for the triangle
##D       side.length <- sqrt(4 * pi / sqrt(3)) * vertex.size / 2
##D 
##D       vertex.frame.color\[vertex.frame.width <= 0\] <- NA
##D       vertex.frame.width\[vertex.frame.width <= 0\] <- 1
##D 
##D       for (i in 1:nrow(coords)) {
##D           x <- coords\[i, 1\]
##D           y <- coords\[i, 2\]
##D           size <- side.length\[i\]
##D           polygon(x + size * c(cos(pi/2), cos(7*pi/6), cos(11*pi/6)),
##D                   y + size * c(sin(pi/2), sin(7*pi/6), sin(11*pi/6)),
##D                   col = vertex.color\[i\], border = vertex.frame.color\[i\],
##D                   lwd = vertex.frame.width\[i\])
##D       }
##D   }
##D 
##D   # Add the triangle shape to igraph
##D   add_shape("triangle", clip = shapes(shape = "circle")$clip, plot = triangle.plot)
##D 
##D   # Example graph
##D   g <- make_ring(10)
##D   V(g)$shape <- rep(c("circle", "triangle"), length.out = vcount(g))
##D   plot(g, vertex.size = 15, vertex.color = "skyblue")
## End(Not run)




### Name: two.factor.analysis
### Title: Two-Factor Analysis for Contingency Tables with Proportion Tests
### Aliases: two.factor.analysis

### ** Examples

## Not run: 
##D # Create example data
##D set.seed(123)
##D treatment <- factor(rep(c("A", "B", "C"), each = 100))
##D outcome <- factor(rbinom(300, 1, rep(c(0.3, 0.5, 0.4), each = 100)),
##D                   labels = c("Failure", "Success"))
##D 
##D # Two-sample proportion tests
##D result <- two.factor.analysis(
##D   y = treatment,
##D   x = outcome,
##D   out.dir = tempdir(),
##D   label = "treatment_analysis"
##D )
##D 
##D # One-sample tests against expected proportion of 0.4
##D result2 <- two.factor.analysis(
##D   y = treatment,
##D   x = outcome,
##D   out.dir = tempdir(),
##D   label = "treatment_vs_expected",
##D   expected.prop = 0.4
##D )
## End(Not run)




### Name: uggmalo.bayesian.bootstrap.with.uncertainty
### Title: Bayesian Bootstrap with Uncertainty for UGGMALO Predictions
### Aliases: uggmalo.bayesian.bootstrap.with.uncertainty

### ** Examples

## Not run: 
##D # WARNING: This is an experimental function with known issues
##D # The following example demonstrates the intended usage but currently fails
##D # with error: "Reference vertex not found in path vertices"
##D 
##D if (requireNamespace("foreach", quietly = TRUE) &&
##D     requireNamespace("doParallel", quietly = TRUE)) {
##D 
##D   # Create example data
##D   n.vertices <- 20
##D   n.perms <- 100
##D 
##D   # Simulated permutation results
##D   perm.results <- matrix(rnorm(n.vertices * n.perms),
##D                          nrow = n.vertices, ncol = n.perms)
##D 
##D   # True predictions
##D   true.predictions <- rnorm(n.vertices)
##D 
##D   # Simple graph structure
##D   graph <- list(
##D     pruned_adj_list = lapply(1:n.vertices, function(i) {
##D       # Simple chain graph
##D       if (i == 1) c(2L)
##D       else if (i == n.vertices) c(i - 1L)
##D       else c(i - 1L, i + 1L)
##D     }),
##D     pruned_dist_list = lapply(1:n.vertices, function(i) {
##D       if (i == 1) c(1.0)
##D       else if (i == n.vertices) c(1.0)
##D       else c(1.0, 1.0)
##D     })
##D   )
##D 
##D   # Run bootstrap (with fewer iterations for example)
##D   # NOTE: This currently fails with a known bug
##D   results <- uggmalo.bayesian.bootstrap.with.uncertainty(
##D     perm.results = perm.results,
##D     true.predictions = true.predictions,
##D     graph = graph,
##D     n.bootstrap = 100,
##D     n.cores = 2
##D   )
##D 
##D   # Plot credible intervals (if successful)
##D   plot(1:n.vertices, true.predictions, pch = 19,
##D        ylim = range(c(results$ci.lower, results$ci.upper)),
##D        xlab = "Vertex", ylab = "Prediction",
##D        main = "UGGMALO Predictions with 95% Credible Intervals")
##D   arrows(1:n.vertices, results$ci.lower,
##D          1:n.vertices, results$ci.upper,
##D          length = 0.05, angle = 90, code = 3)
##D }
## End(Not run)




### Name: uggmalo
### Title: Uniform Grid Graph Model-Averaged Local Linear Regression
###   (UGGMALO)
### Aliases: uggmalo

### ** Examples

## Not run: 
##D # Create a simple chain graph
##D set.seed(123)  # For reproducibility
##D n.vertices <- 20
##D adj.list <- lapply(1:n.vertices, function(i) {
##D   # Simple chain graph
##D   if (i == 1) c(2L)
##D   else if (i == n.vertices) c(i - 1L)
##D   else c(i - 1L, i + 1L)
##D })
##D 
##D dist.list <- lapply(1:n.vertices, function(i) {
##D   if (i == 1) c(1.0)
##D   else if (i == n.vertices) c(1.0)
##D   else c(1.0, 1.0)
##D })
##D 
##D # Create parabolic signal with noise
##D x <- seq(-1, 1, length.out = n.vertices)
##D true.signal <- 3 * x^2 - 2 * x + 1  # Parabola
##D noise <- rnorm(n.vertices, mean = 0, sd = 0.3)
##D y <- true.signal + noise
##D 
##D # Run estimation with default parameters
##D result <- uggmalo(adj.list, dist.list, y)
##D 
##D # Compare true vs estimated values
##D plot(1:n.vertices, y, pch = 19, col = "gray50",
##D      xlab = "Vertex", ylab = "Value",
##D      main = "UGGMALO: True vs Estimated")
##D lines(1:n.vertices, true.signal, col = "blue", lwd = 2)
##D lines(1:n.vertices, result$predictions, col = "red", lwd = 2)
##D legend("topright", c("Observations", "True signal", "UGGMALO estimate"),
##D        col = c("gray50", "blue", "red"),
##D        pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = 2)
##D 
##D # Run with custom parameters for comparison
##D result2 <- uggmalo(adj.list, dist.list, y,
##D                    min.bw.factor = 0.1,
##D                    max.bw.factor = 0.8,
##D                    n.bws = 30,
##D                    verbose = TRUE)
## End(Not run)




### Name: uggmalog
### Title: Uniform Grid Graph Model-Averaged LOGistic regression (UGGMALOG)
### Aliases: uggmalog

### ** Examples

## Not run: 
##D # Create a simple chain graph
##D n <- 10
##D adj <- vector("list", n)
##D weights <- vector("list", n)
##D 
##D # Build chain: 1 -- 2 -- 3 -- ... -- n
##D for (i in 1:n) {
##D   adj[[i]] <- integer(0)
##D   weights[[i]] <- numeric(0)
##D 
##D   if (i > 1) {
##D     adj[[i]] <- c(adj[[i]], i - 1)
##D     weights[[i]] <- c(weights[[i]], 1.0)
##D   }
##D 
##D   if (i < n) {
##D     adj[[i]] <- c(adj[[i]], i + 1)
##D     weights[[i]] <- c(weights[[i]], 1.0)
##D   }
##D }
##D 
##D # Generate some response data
##D set.seed(123)
##D y <- sin(seq(0, pi, length.out = n)) + rnorm(n, sd = 0.1)
##D 
##D # Run UGGMALOG
##D result <- uggmalog(adj, weights, y, verbose = TRUE)
##D 
##D # Print optimal bandwidth index
##D cat("Optimal bandwidth index:", result$opt_bw_idx, "\n")
##D 
##D # More complex example: Grid graph
##D # Create a 5x5 grid graph
##D grid_size <- 5
##D n <- grid_size^2
##D adj <- vector("list", n)
##D weights <- vector("list", n)
##D 
##D # Function to convert 2D coordinates to 1D index
##D coord_to_idx <- function(i, j) (i - 1) * grid_size + j
##D 
##D # Build grid adjacency
##D for (i in 1:grid_size) {
##D   for (j in 1:grid_size) {
##D     idx <- coord_to_idx(i, j)
##D     adj[[idx]] <- integer(0)
##D     weights[[idx]] <- numeric(0)
##D 
##D     # Add horizontal edges
##D     if (j > 1) {
##D       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i, j - 1))
##D       weights[[idx]] <- c(weights[[idx]], 1.0)
##D     }
##D     if (j < grid_size) {
##D       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i, j + 1))
##D       weights[[idx]] <- c(weights[[idx]], 1.0)
##D     }
##D 
##D     # Add vertical edges
##D     if (i > 1) {
##D       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i - 1, j))
##D       weights[[idx]] <- c(weights[[idx]], 1.0)
##D     }
##D     if (i < grid_size) {
##D       adj[[idx]] <- c(adj[[idx]], coord_to_idx(i + 1, j))
##D       weights[[idx]] <- c(weights[[idx]], 1.0)
##D     }
##D   }
##D }
##D 
##D # Generate response based on distance from center
##D center <- (grid_size + 1) / 2
##D y <- numeric(n)
##D for (i in 1:grid_size) {
##D   for (j in 1:grid_size) {
##D     dist_from_center <- sqrt((i - center)^2 + (j - center)^2)
##D     y[coord_to_idx(i, j)] <- exp(-dist_from_center / 2) + rnorm(1, sd = 0.05)
##D   }
##D }
##D 
##D # Run UGGMALOG with custom parameters
##D result <- uggmalog(
##D   adj.list = adj,
##D   weight.list = weights,
##D   y = y,
##D   n.bws = 30,
##D   kernel.type = 5,  # Tricube kernel
##D   fit.quadratic = TRUE,
##D   verbose = TRUE
##D )
## End(Not run)




### Name: ulogit
### Title: Fit Univariate Logistic Regression Model
### Aliases: ulogit

### ** Examples

# Basic usage with simulated data
set.seed(123)
x <- seq(0, 1, length.out = 100)
true_prob <- 1/(1 + exp(-(2*x - 1)))
y <- rbinom(100, 1, prob = true_prob)
fit <- ulogit(x, y)

# Plot results
plot(x, y, pch = 16, col = ifelse(y == 1, "blue", "red"),
     main = "Univariate Logistic Regression")
lines(x, fit$predictions, lwd = 2)
legend("topleft", c("y = 1", "y = 0", "Fitted"),
       col = c("blue", "red", "black"),
       pch = c(16, 16, NA), lty = c(NA, NA, 1))

# Example with weights
w <- runif(100, 0.5, 1.5)
fit_weighted <- ulogit(x, y, w = w)

# Example with increased regularization
fit_regularized <- ulogit(x, y, ridge.lambda = 0.1)

# Compare LOOCV errors
cat("Standard model LOOCV error:", mean(fit$errors), "\n")
cat("Weighted model LOOCV error:", mean(fit_weighted$errors), "\n")
cat("Regularized model LOOCV error:", mean(fit_regularized$errors), "\n")




### Name: univariate.gkmm
### Title: Adaptive Neighborhood Size Graph K-Means for Univariate Data
### Aliases: univariate.gkmm

### ** Examples

## Not run: 
##D # Generate example data
##D set.seed(123)
##D x <- seq(0, 10, length.out = 100)
##D y <- sin(x) + rnorm(100, 0, 0.1)
##D 
##D # Fit model
##D result <- univariate.gkmm(x, y, h.min = 2, h.max = 20, n.bb = 100)
##D 
##D # Plot results
##D plot(result)
##D 
##D # Summary
##D summary(result)
## End(Not run)




### Name: upgmalo
### Title: Univariate Path Graph Model Averaging Local Linear Model
### Aliases: upgmalo

### ** Examples

## Not run: 
##D # Generate nonlinear data
##D n <- 50
##D x <- seq(0, 2*pi, length.out = n)
##D y_true <- sin(x)
##D y <- y_true + rnorm(n, 0, 0.2)
##D 
##D # Fit model
##D   fit <- upgmalo(x, y, y.true = y_true, h.max = 20, n.CVs = 20)
##D 
##D   # Plot results
##D   plot(fit)
## End(Not run)




### Name: validate.maximal.packing
### Title: Validate a Maximal Packing
### Aliases: validate.maximal.packing

### ** Examples

## Not run: 
##D # Create a simple path graph with 6 vertices
##D adj.list <- list(
##D   c(2),           # vertex 1 connects to 2
##D   c(1, 3),        # vertex 2 connects to 1 and 3
##D   c(2, 4),        # vertex 3 connects to 2 and 4
##D   c(3, 5),        # vertex 4 connects to 3 and 5
##D   c(4, 6),        # vertex 5 connects to 4 and 6
##D   c(5)            # vertex 6 connects to 5
##D )
##D weight.list <- list(
##D   c(1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1)
##D )
##D 
##D # Test a packing with vertices 1 and 4
##D packing <- c(1, 4)
##D radius <- 3
##D 
##D # Validate the packing
##D result <- validate.maximal.packing(adj.list, weight.list, packing, radius)
##D print(result)
## End(Not run)




### Name: vars.approx.monotonically.assoc.with.geodesic
### Title: Given TA and TAA from E.geodesic.X(), this routine identifies
###   variables monotonically associated with the distance along the
###   geodesic
### Aliases: vars.approx.monotonically.assoc.with.geodesic

### ** Examples

# Assuming X has columns TA and TAA from E.geodesic.X output
# monotonic_vars <- vars.approx.monotonically.assoc.with.geodesic(X, eps = 0.4)
# X[monotonic_vars, ] # subset to monotonic variables




### Name: verify.maximal.packing
### Title: Verify Maximal Packing Created by create.maximal.packing
### Aliases: verify.maximal.packing

### ** Examples

## Not run: 
##D # Create a simple cycle graph
##D n <- 10
##D adj.list <- lapply(1:n, function(i) {
##D   c(ifelse(i == 1, n, i - 1), ifelse(i == n, 1, i + 1))
##D })
##D weight.list <- lapply(1:n, function(i) c(1, 1))
##D 
##D # Create maximal packing
##D result <- create.maximal.packing(adj.list, weight.list, grid.size = 3)
##D 
##D # Verify the packing with detailed output
##D is_valid <- verify.maximal.packing(result, verbose = TRUE)
##D 
##D # Verify quietly
##D is_valid <- verify.maximal.packing(result, verbose = FALSE)
## End(Not run)




### Name: vert.error.bar
### Title: Add Vertical Error Bar to Plot
### Aliases: vert.error.bar

### ** Examples

plot(1:10, rnorm(10), ylim = c(-3, 3))
vert.error.bar(5, -1, 1, col = "red")




### Name: vertices.local_extrema
### Title: Extract Vertices of a Specific Local Extremum
### Aliases: vertices.local_extrema

### ** Examples

# Create example data
adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
y <- c(1, 3, 2, 5, 1)

# Detect maxima
maxima <- detect.local.extrema(adj.list, weight.list, y, 2, 2)

# Extract vertices for the first maximum (if it exists)
if (length(maxima$vertices) > 0) {
  v <- vertices(maxima, maxima$labels[1])
}




### Name: visualize.grid.function
### Title: Visualize a Function on a Grid Graph
### Aliases: visualize.grid.function

### ** Examples

## Not run: 
##D grid.size <- 20
##D x <- rep(1:grid.size, grid.size) / grid.size
##D y <- rep(1:grid.size, each = grid.size) / grid.size
##D z <- exp(-10*((x-0.3)^2 + (y-0.3)^2)) + 0.5*exp(-8*((x-0.7)^2 + (y-0.6)^2))
##D centers <- which(z > 0.9)
##D visualize.grid.function(grid.size, z, centers)
##D 
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
##D   visualize.grid.function(grid.size, z, centers)
##D }
## End(Not run)



### Name: visualize.smoothing.steps
### Title: Visualize Smoothing Steps from Gradient Flow Complex
### Aliases: visualize.smoothing.steps

### ** Examples

## Not run: 
##D if (requireNamespace("rgl", quietly = TRUE)) {
##D   # result <- create.gflow.cx(..., detailed.recording = TRUE)
##D   # plot_res <- list(layout = <n x 2 matrix>, graph = <igraph optional>)
##D   # visualize.smoothing.steps(result, plot_res)
##D }
## End(Not run)



### Name: wasserstein.distance.1D
### Title: Calculate Wasserstein Distance Between 1D Samples
### Aliases: wasserstein.distance.1D

### ** Examples

x <- rnorm(1000)
y <- rnorm(1000, mean = 1)
dist <- wasserstein.distance.1D(x, y)
print(dist)




### Name: wasserstein.distance
### Title: Compute Wasserstein Distance Between Two Datasets
### Aliases: wasserstein.distance

### ** Examples

X <- matrix(rnorm(1000), ncol = 2)
Y <- matrix(rnorm(1000, mean = 1), ncol = 2)
result <- wasserstein.distance(X, Y)
print(result)




### Name: wasserstein.divergence
### Title: Compute Wasserstein Divergence Between Two Point Sets
### Aliases: wasserstein.divergence

### ** Examples

## Not run: 
##D X <- matrix(rnorm(1000 * 3), nrow = 1000, ncol = 3)
##D Y <- matrix(rnorm(1000 * 3), nrow = 1000, ncol = 3)
##D div <- wasserstein.divergence(X, Y, k = 10)
## End(Not run)




### Name: wasserstein1d.test
### Title: Performs Permutation Test for the Wasserstein Distance Between
###   Two 1D Samples
### Aliases: wasserstein1d.test

### ** Examples

## Two samples from the same distribution
x <- rnorm(100); y <- rnorm(100)
res <- wasserstein1d.test(x, y, n.perms = 1000, n.cores = 1)
str(res)

## Two samples from different distributions
x <- rnorm(100)
y <- rnorm(100, mean = 2)
## Perform the permutation test with 1000 permutations, serially
res <- wasserstein1d.test(x, y, n.perms = 1000, n.cores = 1)
str(res)
## Perform the permutation test with 1000 permutations, using 2 cores
res <- wasserstein1d.test(x, y, n.perms = 1000, n.cores = 2)
str(res)



### Name: weighted.p.value
### Title: Weighted P-value Calculation
### Aliases: weighted.p.value

### ** Examples

# Example 1: Uncertainty represented by a normal distribution
u_sample <- rnorm(1000, mean = 1.5, sd = 0.5)

# Test against null hypothesis N(0, 1)
p_value <- weighted.p.value(u_sample, mu = 0, sigma = 1)
cat("Weighted p-value:", p_value, "\n")

# Compare with classical p-value using the mean
classical_p <- pnorm(mean(u_sample), mean = 0, sd = 1, lower.tail = FALSE)
cat("Classical p-value:", classical_p, "\n")

# Example 2: Two-sided test
u_sample2 <- rnorm(1000, mean = -0.5, sd = 0.3)
p_two_sided <- weighted.p.value(u_sample2, mu = 0, sigma = 1,
                                 alternative = "two.sided")
cat("Two-sided weighted p-value:", p_two_sided, "\n")

# Example 3: Using with bootstrap samples
# Simulate some data and bootstrap the mean
set.seed(123)
original_data <- rnorm(30, mean = 1.2, sd = 2)
boot_means <- replicate(1000, mean(sample(original_data, replace = TRUE)))

# Weighted p-value using bootstrap distribution
p_boot <- weighted.p.value(boot_means, mu = 0, sigma = 2/sqrt(30))
cat("Bootstrap-based weighted p-value:", p_boot, "\n")




### Name: weighted.p.value.summary
### Title: Weighted P-value Summary
### Aliases: weighted.p.value.summary

### ** Examples

u_sample <- rnorm(1000, mean = 1.5, sd = 0.5)
summary_results <- weighted.p.value.summary(u_sample, mu = 0, sigma = 1)
print(summary_results)




### Name: wgraph.prune.long.edges
### Title: Prune Long Edges in a Weighted Graph
### Aliases: wgraph.prune.long.edges

### ** Examples

## Not run: 
##D # Create a simple weighted graph
##D graph <- list(c(2,3), c(1,3), c(1,2))
##D edge.lengths <- list(c(1,2), c(1,3), c(2,3))
##D threshold <- 0.9
##D 
##D # Prune the graph
##D result <- wgraph.prune.long.edges(graph, edge.lengths, threshold,
##D                                   use.total.length.constraint = TRUE,
##D                                   verbose = TRUE)
##D # Examine the results
##D print(result$adj_list)
##D print(result$edge_lengths_list)
##D print(result$path_lengths)
##D print(result$edge_lengths)
## End(Not run)



### Name: winsorize
### Title: Winsorize a numeric vector
### Aliases: winsorize

### ** Examples

# Simple example
x <- 1:10
winsorize(x, p = 0.2)

# With normally distributed data
set.seed(123)
y <- rnorm(100)
y_wins <- winsorize(y, p = 0.05)

# Compare ranges
range(y)
range(y_wins)




### Name: winsorize.zscore
### Title: Winsorized Z-score normalization
### Aliases: winsorize.zscore

### ** Examples

# Example with random data
set.seed(123)
example.data <- matrix(rnorm(100, 5, 2), ncol=5)
# Add some outliers
example.data[1,1] <- 25
example.data[2,3] <- -15
normalized.data <- winsorize.zscore(example.data)



