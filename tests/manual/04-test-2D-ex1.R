## ---------------------------------------------------------------------------------------------
## gflow adaptation of ~/current_projects/msr2/R_correctness_tests/test_v4_create_basin_cx.R
## ---------------------------------------------------------------------------------------------

## =================================================================================
## 2D Example 1: Gaussian Mixture on Uniform Grid with Added Noise
## =================================================================================
## Setup: 2D uniform sample from [0,1]^2 with Gaussian mixture response + noise
## =================================================================================

## ---- Function and Grid Setup ----
f1 <- create.gaussian.mixture(A1 = 0.25, A2 = 1, sigma = 0.5)
axis.n.pts <- 100
grid <- create.grid(axis.n.pts)
f1.res <- evaluate.function.on.grid.as.vector(f1$f, grid)

## ---- Extract Full Grid Data ----
X.full <- f1.res$X
y.smooth.full <- f1.res$y.smooth
cat(sprintf("Full grid: %d points, dim = %d\n", nrow(X.full), ncol(X.full)))

## ---- Subsample Data ----
n.sample <- 250
set.seed(123)  # For reproducibility
sample.idx <- sample(nrow(X.full), size = n.sample)
X <- X.full[sample.idx, ]
y.smooth <- y.smooth.full[sample.idx]

## ---- Add Noise ----
noise.sigma <- 0.1
set.seed(124)  # Separate seed for noise
eps <- rnorm(n.sample, mean = 0, sd = noise.sigma)
y <- y.smooth + eps

## ---- Create iKNN Graphs ----
k.range <- c(min = 2, max = 10)

ex1.graphs <- create.iknn.graphs(
  X,
  kmin = k.range["min"],
  kmax = k.range["max"],
  max.path.edge.ratio.deviation.thld = 0.1,
  path.edge.ratio.percentile = 0.5,
  compute.full = TRUE,
  pca.dim = 100,
  variance.explained = 0.99,
  n.cores = NULL,
  verbose = TRUE
)

## ---- Analyze Graph Statistics ----
ex1.graphs.stats <- summary(ex1.graphs)
##   idx  k n_ccomp edges mean_degree min_degree max_degree sparsity
##     1  2      12   131        2.62          1          5  0.97354
##     2  3       1   189        3.78          1          9  0.96182
##   ...

## ---- Select Graph for Analysis ----
selected.graph.idx <- 3
ex1.graph <- ex1.graphs$geom_pruned_graphs[[selected.graph.idx]]
cat(sprintf("Selected graph: k = %d\n",
            ex1.graphs.stats$k[selected.graph.idx]))

## ---- Save Results ----
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output.dir <- "~/current_projects/gflow/tests/manual/data"
output.file <- file.path(output.dir,
                         sprintf("2D_example_1_%s.rda", timestamp))

# Ensure directory exists
if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

save(
  # Function and grid
  f1.res, grid, axis.n.pts,
  # Data
  X.full, y.smooth.full,
  X, y.smooth, y, eps,
  # Parameters
  n.sample, noise.sigma, k.range,
  # Results
  ex1.graphs, ex1.graphs.stats,
  selected.graph.idx, ex1.graph,
  file = output.file
)
cat(sprintf("Saved to: %s\n", output.file))
## Saved to: ~/current_projects/gflow/tests/manual/data/2D_example_1_20251008_183721.rda # n = 100
## Saved to: ~/current_projects/gflow/tests/manual/data/2D_example_1_20251008_184407.rda # n = 250

output.file <- "~/current_projects/gflow/tests/manual/data/2D_example_1_20251008_183721.rda"
## ---- Load Results (for continuation) ----
## load(output.file)

ex1.ggraph <- ggraph(ex1.graph$adj_list, ex1.graph$weight_list)

## ---- View (X, y.smooth) ----
plot3D.graph(ex1.ggraph,
             z = 10*y.smooth,
             layout = "kk",
             conn.points = TRUE,
             use.spheres = FALSE,
             graph.alpha = 0.7,
             z.point.size = 1,
             z.color = NULL,
             z.alpha = 1,
             edge.color = "gray70",
             edge.width = 1,
             base.plane = TRUE,
             base.vertex.size = 0.5,
             z.scale = 1,
             show.axes = TRUE,
             vertical.lines = FALSE,
             vertical.line.style = "dashed",
             dash.length = 0.05,
             gap.length = 0.05,
             vertical.line.color = "darkgray",
             vertical.line.alpha = 0.5,
             vertical.line.width = 0.5,
             basins.df = NULL,
             evertex.sphere.radius = 0.2,
             evertex.min.color = "blue",
             evertex.max.color = "red",
             evertex.cex = 1,
             evertex.adj = c(0.5, 0.5),
             evertex.label.offset = 0.3)

ex1.y.smooth.lextr.nbhds <- compute.extrema.hop.nbhds(ex1.graph$adj_list,
                                                      ex1.graph$weight_list,
                                                      y.smooth)
extr.df <- ex1.y.smooth.lextr.nbhds$extrema_df

extr.df[extr.df$hop_idx > 8,]
##    vertex hop_idx is_max label     value spurious
## 8      76      21      0    m2 0.3874571    FALSE
## 9     120      14      0    m3 0.3927188    FALSE
## 12    153     Inf      0    m1 0.3871638    FALSE
## 14    116       9      0   m11 0.8619050    FALSE
## 20    181     Inf      1    M1 1.1000452    FALSE

## ---- View (X, y) ----

plot3D.graph(ex1.ggraph,
             z = 10*y,
             layout = "kk",
             conn.points = TRUE,
             use.spheres = FALSE,
             graph.alpha = 0.7,
             z.point.size = 1,
             z.color = NULL,
             z.alpha = 1,
             edge.color = "gray70",
             edge.width = 1,
             base.plane = TRUE,
             base.vertex.size = 0.5,
             z.scale = 1,
             show.axes = TRUE,
             vertical.lines = FALSE,
             vertical.line.style = "dashed",
             dash.length = 0.05,
             gap.length = 0.05,
             vertical.line.color = "darkgray",
             vertical.line.alpha = 0.5,
             vertical.line.width = 0.5,
             basins.df = NULL,
             evertex.sphere.radius = 0.2,
             evertex.min.color = "blue",
             evertex.max.color = "red",
             evertex.cex = 1,
             evertex.adj = c(0.5, 0.5),
             evertex.label.offset = 0.3)


## ---- Estimate E[y|G] using riem_dcx ----

## # Option 1: Statistical convention
## oracle.error <- c()      # vs true smooth
## training.error <- c()    # vs observations

## # Option 2: ML convention
## true.mse <- c()          # vs y.smooth
## observed.mse <- c()      # vs y

## # Option 3: Explicit about what's being measured
## bias.squared <- c()      # (y.smooth - fitted)^2
## rss <- c()              # residual sum of squares (y - fitted)^2

## # Option 4: If thinking about decomposition
## smoothness.error <- c()  # how well you recover the smooth function
## noise.accommodation <- c()  # how much you fit the noise

gamma.grid = c(0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.45, 0.5)
oracle.error <- c()
training.error <- c()
gcv.score <- c()
fit.list <- list()
for (i in seq_along(gamma.grid)) {
    cat("\n----------------------------\n")
    cat("i:",i,"\n----------------------------\n")
    fit <- fit.knn.riem.graph.regression(X,
                                         y,
                                         k = 3,
                                         response.penalty.exp = gamma.grid[i],
                                         density.uniform.pull = 0.1,
                                         n.eigenpairs = 50,
                                         max.iterations = 10,
                                         filter.type = "cubic_spline",
                                         use.counting.measure = TRUE,
                                         density.normalization = 0,
                                         t.diffusion = 10,
                                         epsilon.y = 1e-4,
                                         epsilon.rho = 1e-4,
                                         max.ratio.threshold = 0.1,
                                         threshold.percentile = 0.5,
                                         test.stage = -1)
    fit.list[[i]] <- fit
    ## Error against the true smooth function (oracle error)
    oracle.error[i] <- sum((y.smooth - fit$fitted.values)^2)
    ## Residual sum of squares from the noisy observations
    training.error[i] <- sum((y - fit$fitted.values)^2)
    ## GCV error
    gcv.score[i] <- min(fit$gcv$gcv.optimal)
}

output.file <- "~/current_projects/gflow/tests/manual/data/rdcx_fits_of_2D_example_1_20251008_184407.rda"
save(fit.list,
     oracle.error,
     training.error,
     gcv.score,
     file = output.file
)
cat(sprintf("Saved to: %s\n", output.file))
## Saved to: ~/current_projects/gflow/tests/manual/data/rdcx_fits_of_2D_example_1_20251008_184407.rda

plot(gamma.grid, oracle.error, type = "l", main = "oracle.error", ylab = "oracle.error", xlab = "gamma.grid")
plot(gamma.grid, training.error, type = "l", main = "training.error", ylab = "training.error", xlab = "gamma.grid")
plot(gamma.grid, gcv.score, type = "l", main = "gcv.score", ylab = "gcv.score", xlab = "gamma.grid")

file <- "~/current_projects/gflow/tests/manual/pics/rdcx_fits_of_2D_example_1_errors_"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=18, height=6)
op <- par(mfrow=c(1,3), mar=c(3.5,3.5,1.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot(gamma.grid, oracle.error, type = "l", main = "oracle.error", ylab = "oracle.error", xlab = "gamma.grid")
plot(gamma.grid, training.error, type = "l", main = "training.error", ylab = "training.error", xlab = "gamma.grid")
plot(gamma.grid, gcv.score, type = "l", main = "gcv.score", ylab = "gcv.score", xlab = "gamma.grid")
par(op)
dev.off()
system(paste0("open ",file))

## ---- Selecting optimal k ----

k.oracle.error <- c()
k.training.error <- c()
k.gcv.score <- c()
k.fit.list <- list()
i <- 1
for (k in 3:20) {
    cat("\n----------------------------\n")
    cat("k",k,"\n----------------------------\n")
    fit <- fit.knn.riem.graph.regression(X,
                                         y,
                                         k = k,
                                         response.penalty.exp = 0.3,
                                         density.uniform.pull = 0.1,
                                         n.eigenpairs = 50,
                                         max.iterations = 10,
                                         filter.type = "cubic_spline",
                                         use.counting.measure = TRUE,
                                         density.normalization = 0,
                                         t.diffusion = 10,
                                         epsilon.y = 1e-4,
                                         epsilon.rho = 1e-4,
                                         max.ratio.threshold = 0.1,
                                         threshold.percentile = 0.5,
                                         test.stage = -1)
    k.fit.list[[i]] <- fit
    ## Error against the true smooth function (oracle error)
    k.oracle.error[i] <- sum((y.smooth - fit$fitted.values)^2)
    ## Residual sum of squares from the noisy observations
    k.training.error[i] <- sum((y - fit$fitted.values)^2)
    ## GCV error
    k.gcv.score[i] <- min(fit$gcv$gcv.optimal)
    i <- i + 1
}

output.file <- "~/current_projects/gflow/tests/manual/data/rdcx_fits_by_k_of_2D_example_1_20251008_184407.rda"
save(k.fit.list,
     k.oracle.error,
     k.training.error,
     k.gcv.score,
     file = output.file
)
cat(sprintf("Saved to: %s\n", output.file))
## Saved to: ~/current_projects/gflow/tests/manual/data/rdcx_fits_by_k_of_2D_example_1_20251008_184407.rda

op <- par(mfrow=c(1,3), mar=c(3.5,3.5,1.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot(3:20,k.oracle.error, type = "l", main = "oracle.error", ylab = "oracle.error", xlab = "k")
plot(3:20,k.training.error, type = "l", main = "training.error", ylab = "training.error", xlab = "k")
plot(3:20,k.gcv.score, type = "l", main = "gcv.score", ylab = "gcv.score", xlab = "k")
par(op)

file <- "~/current_projects/gflow/tests/manual/pics/rdcx_fits_by_k_of_2D_example_1_errors_"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=18, height=6)
op <- par(mfrow=c(1,3), mar=c(3.5,3.5,1.5,0.5), mgp=c(2.0,0.4,0),tcl = -0.3)
plot(3:20,k.oracle.error, type = "l", main = "oracle.error", ylab = "oracle.error", xlab = "k")
plot(3:20,k.training.error, type = "l", main = "training.error", ylab = "training.error", xlab = "k")
plot(3:20,k.gcv.score, type = "l", main = "gcv.score", ylab = "gcv.score", xlab = "k")
par(op)
dev.off()
system(paste0("open ",file))
