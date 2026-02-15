## 2D synthetic data
##
## based on the following past work -> ./R_correctness_tests/test_v4_create_basin_cx.R
##
## [~/current_projects/msr2]% find . -type f -name "*.R" -exec grep -lq 'harmonic' {} \; -print | xargs ls -lht
## find . -type f -name "*.R" -exec grep -lq 'harmonic' {} \; -print | xargs ls -lht
## -rw-r--r--  1 pgajer  staff    37K Jun 24 16:29 ./R/gflow_cx_utils.R
## -rw-r--r--  1 pgajer  staff    23K May 16 14:22 ./R/gflow_cx.R
## -rw-r--r--  1 pgajer  staff    13K May 16 14:22 ./R_correctness_tests/test_v4_create_basin_cx.R
## -rw-r--r--  1 pgajer  staff    13K May 16 12:32 ./R_correctness_tests/test_v3_create_basin_cx.R
## -rw-r--r--  1 pgajer  staff   109K May 15 14:48 ./R/basin_cx.R
## -rw-r--r--  1 pgajer  staff   9.7K May 15 14:16 ./R_correctness_tests/test_v2_create_basin_cx.R
## -rw-r--r--  1 pgajer  staff    67K May 13 12:53 ./R_correctness_tests/test_create_basin_cx.R
## -rw-r--r--  1 pgajer  staff    84K May  2 03:33 ./R/Archives/gflow_basins_M2.R
## -rw-r--r--  1 pgajer  staff    19K May  1 09:52 ./R/harmonic_smoother.R
## -rw-r--r--  1 pgajer  staff    37K Apr 30 11:14 ./R/Archives/gflow_basins_A30.R
## -rw-r--r--  1 pgajer  staff   8.4K Apr 29 09:38 ./R_correctness_tests/test_find_local_extrema.R
## -rw-r--r--  1 pgajer  staff    15K Apr 25 19:50 ./R/amagelo.R

## ---------------------------------------------------------------------------------------------
## gflow adaptation of ~/current_projects/msr2/R_correctness_tests/test_v4_create_basin_cx.R
## ---------------------------------------------------------------------------------------------

f1 <- create.gaussian.mixture(A1 = 0.25, A2 = 1, sigma = 0.5)
axis.n.pts <- 100
grid <- create.grid(axis.n.pts)
f1.res <- evaluate.function.on.grid.as.vector(f1$f, grid)
str(f1.res)
X <- f1.res$X
dim(X)
## [1] 10000     2
y.smooth <- f1.res$y.smooth
length(y.smooth)
## f1.grid <- evaluate.function.on.grid(f1$f, grid)
## image(f1.grid)
n <- 500
ii <- sample(nrow(X), size = n)
X <- X[ii,]
y.smooth <- y.smooth[ii]

sigma <- .1
eps <- rnorm(n, 0, sigma)
y <- y.smooth + eps

sigma <- .05
eps <- rnorm(n, 0, sigma)
y2 <- y.smooth + eps

pruning.thld <- 0.01
k.min <- 2
k.max <- 20
## graphs.res <- create.iknn.graphs(X,
##                                  kmin = k.min,
##                                  kmax = k.max,
##                                  pruning.thld = pruning.thld,
##                                  verbose = TRUE)
graphs.res <- create.iknn.graphs(X,
                                 kmin = k.min,
                                 kmax = k.max,
                                 max.path.edge.ratio.deviation.thld = 0.1,
                                 path.edge.ratio.percentile = 0.5,
                                 compute.full = TRUE,
                                 pca.dim = 100,
                                 variance.explained = 0.99,
                                 n.cores = NULL,
                                 verbose = TRUE)

graphs.stats <- summary(graphs.res)
## idx  k n_ccomp edges mean_degree min_degree max_degree sparsity
##   1  2      27   583        2.33          1          5  0.99533
##   2  3       2   750        3.00          1          6  0.99399
##   3  4       1   894        3.58          1          7  0.99283
##   4  5       1  1045        4.18          1          8  0.99162
##   5  6       1  1150        4.60          1          8  0.99078
##   6  7       1  1241        4.96          1          9  0.99005
##   7  8       1  1327        5.31          1         10  0.98936
##   8  9       1  1403        5.61          1         10  0.98875
##   9 10       1  1472        5.89          1         10  0.98820

## timestamp <- generate.timestamp()
## file <- "~/current_projects/msr2/data/2d_gm_example__"
## file <- paste0(file, timestamp, ".rda")
file <- "~/current_projects/msr2/data/2d_gm_example__2025_05_13_093816.rda"
## save(f1.res, X, y.smooth, y, eps, sigma, k.min, k.max, graphs.res, graph.idx, graph, file = file)
load(file)

plot(X, xlab = "", ylab = "")

graph.idx <- 4
graph <- graphs.res$geom_pruned_graphs[[graph.idx]]

g <- ggraph(graph$adj_list, graph$weight_list)
plot.res <- plot(g, vertex.size = 5, vertex.radius = 0.1, y = y.smooth, dim = 2)

##
## mst completion graph
##
cmst.graph <- create.cmst.graph(X,
                                q.thld = 0.999,
                                verbose = TRUE)
cmst.ggraph <- ggraph(cmst.graph$cmst_adj_list, cmst.graph$cmst_weight_list)
plot(cmst.ggraph, vertex.size = 5, use.saved.layout = plot.res$layout, y = y.smooth, dim = 2)

b.cx <- create.basin.cx(graph$adj_list,
                        graph$weight_list,
                        y.smooth,
                        basin.merge.overlap.thld = 0.1,
                        min.asc.desc.cell.size.thld = 10,
                        min.asc.asc.cell.size.thld = 10,
                        min.desc.desc.cell.size.thld = 20)

b.cx$initial_basin_cx$basins_df

plot3D.graph(g,
             z = 20*y.smooth,
             layout = "kk",
             conn.points = TRUE,
             use.spheres = TRUE,
             graph.alpha = 0.7,
             z.point.size = 15,
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

lextr.nbhds <- compute.extrema.hop.nbhds(graph$adj_list,
                                         graph$weight_list,
                                         y.smooth)

names(lextr.nbhds) # "lmin_hop_nbhds" "lmax_hop_nbhds" "extrema_df"

lextr.nbhds$extrema_df
##    vertex hop_idx is_max label  fn_value
## 1     490       2      0   m13 0.9330497
## 2       9       8      0   m11 0.8597970
## 3     451       1      0   m10 0.5350543
## 4     478      11      0    m3 0.4011693
## 5     117       1      0    m5 0.4405109
## 6     129      12      0    m2 0.3927716
## 7     144       2      0    m8 0.4826867
## 8     197       1      0    m9 0.4859079
## 9     343       1      0    m6 0.4505406
## 10     76       2      0    m7 0.4788771
## 11    395       1      0    m4 0.4141412
## 12    406     Inf      0    m1 0.3257621
## 13    114       1      0   m14 0.9728258
## 14    462       1      0   m12 0.9319622
## 15    471       2      1    M3 0.6092218
## 16    243       1      1    M2 0.9561888
## 17    155     Inf      1    M1 1.1011285
## 18     75       1      1    M4 0.4524172

lextr.nbhds$extrema_df[lextr.nbhds$extrema_df$hop_idx > 2,]
##    vertex hop_idx is_max label  fn_value
## 2       9       8      0   m11 0.8597970
## 4     478      11      0    m3 0.4011693
## 6     129      12      0    m2 0.3927716
## 12    406     Inf      0    m1 0.3257621
## 17    155     Inf      1    M1 1.1011285

lextr.nbhds$extrema_df[lextr.nbhds$extrema_df$hop_idx > 8,]

plot(y.smooth, y, las = 1)


## ------------------------------------------------------------------------------------------
##
## Extracting local nbhds of each local extremum of y.smooth as even in the
## smooth function case there are spurious local extrema due to the topology of
## the underlying graph
##
## ------------------------------------------------------------------------------------------

subgraphs.res <- extract.extrema.subgraphs(graph$adj_list,
                                           graph$weight_list,
                                           y = y.smooth,
                                           predictions = y.smooth,
                                           lextr.nbhds$extrema_df,
                                           hop.offset = 3,
                                           min.hop.idx = 0,
                                           max.hop.idx = Inf)

length(subgraphs.res) # 18
nrow(lextr.nbhds$extrema_df) # 18
names(subgraphs.res[[1]])

i <- 8
subgraph <- subgraphs.res[[i]]

gg <- ggraph(subgraph$adj_list, subgraph$weight_list)
sg.plot.res <- plot(gg, vertex.size = 5, y = subgraph$y)
subgraph$extremum_info
##   vertex hop_idx is_max label  fn_value
## 1     72       2      0   m13 0.9330497

subgraph$extrema_df
##    vertex hop_idx is_max label  fn_value
## 1      72       2      0   m13 0.9330497
## 2       2       8      0   m11 0.8597970
## 16     39       1      1    M2 0.9561888

extrema.df <- subgraph$extrema_df
colnames(extrema.df)[1] <- "evertex"

plot3D.graph(sg.plot.res,
             subgraph$predictions,
             base.vertex.size = 0.05,
             z.point.size = 0.05,
             basins.df = extrema.df,
             evertex.sphere.radius = 0.1,
             evertex.min.color = "blue",
             evertex.max.color = "red",
             evertex.cex = 1,
             evertex.adj = c(0.5, 0.5),
             evertex.label.offset = 0.01)

## method = c("weighted_mean", "harmonic_iterative", "harmonic_eigen", "biharmonic_harmonic", "boundary_smoothed")
sg.hext.res <- apply.harmonic.extension(subgraph,
                                        method = "weighted_mean",
                                        ##method = "harmonic_iterative",
                                        ##method = "harmonic_eigen",
                                        ##method = "biharmonic_harmonic",
                                        ## method = "boundary_smoothed",
                                        max_iterations = 100,
                                        tolerance = 1e-6,
                                        sigma = 1.0,
                                        record_iterations = TRUE,
                                        verbose = TRUE)

## names(sg.hext.res)
## [1] "original"         "smoothed"         "iterations"       "method"
## [5] "extremum_info"    "convergence_info"
## sg.hext.res$extremum_info
## plot(sg.hext.res$original, sg.hext.res$smoothed, las = 1)
plot.graph.3d(sg.plot.res,
              sg.hext.res$smoothed,
              base.vertex.size = 0.05,
              z.point.size = 0.05,
              basins.df = extrema.df,
              evertex.sphere.radius = 0.1,
              evertex.min.color = "blue",
              evertex.max.color = "red",
              evertex.cex = 1,
              evertex.adj = c(0.5, 0.5),
              evertex.label.offset = 0.01)

## -----------------------------------------------------------------------------------------------------
##
## Applying quasi-harmonic smoothing to all spurious local extrema
##
## Question 1: Does the process create new local extrema?
## Question 2: If yes, does iterative application of the qusi-harmonic smoothing fix the issue?
## Question 3: Would it help to extend the boundary set to a few hops away from the current boundary?
##
## -----------------------------------------------------------------------------------------------------

gflow.cx.res <- create.gflow.cx(graph$adj_list,
                                graph$weight_list,
                                y = y.smooth,
                                hop.idx.thld = 5,
                                smoother.type = 0,  # Default: Weighted Mean
                                max.outer.iterations = 5,
                                max.inner.iterations = 100,
                                smoothing.tolerance = 1e-6,
                                sigma = 1.0,
                                process.in.order = TRUE,
                                verbose = TRUE,
                                detailed.recording = TRUE)

names(gflow.cx.res)
## [1] "harmonic_predictions" "lmin_hop_nbhds"       "lmax_hop_nbhds"
## [4] "smoothing_history"    "extrema_df"           "extrema_df2"

gflow.cx.res$extrema_df2

he.lextr.nbhds <- compute.extrema.hop.nbhds(graph$adj_list,
                                            graph$weight_list,
                                            gflow.cx.res$harmonic_predictions)

he.lextr.nbhds$extrema_df
##   vertex hop_idx is_max label spurious     value
## 1    478      11      0  min1    FALSE 0.4011693
## 2    245       1      0  min2     TRUE 0.5974583
## 3    129      12      0  min3    FALSE 0.3927716
## 4    406     Inf      0  min4    FALSE 0.3257621
## 5    107       2      0  min5     TRUE 0.4405109
## 6     78       1      0  min6     TRUE 0.4907489
## 7      9       8      0  min7    FALSE 0.8597970
## 8    155     Inf      1  max1    FALSE 1.1011285

i2.gflow.cx.res <- create.gflow.cx(graph$adj_list,
                                graph$weight_list,
                                gflow.cx.res$harmonic_predictions,
                                hop.idx.thld = 5,
                                smoother.type = 0,  # Default: Weighted Mean
                                max.outer.iterations = 5,
                                max.inner.iterations = 100,
                                smoothing.tolerance = 1e-6,
                                sigma = 1.0,
                                process.in.order = TRUE,
                                verbose = TRUE,
                                detailed.recording = TRUE)

i2.he.lextr.nbhds <- compute.extrema.hop.nbhds(graph$adj_list,
                                               graph$weight_list,
                                               i2.gflow.cx.res$harmonic_predictions)
i2.he.lextr.nbhds$extrema_df
##   vertex hop_idx is_max label spurious     value
## 1    406     Inf      0  min1    FALSE 0.3257621
## 2    478      11      0  min2    FALSE 0.4011693
## 3    231       2      0  min3     TRUE 0.5333015
## 4    129      12      0  min4    FALSE 0.3927716
## 5      9       8      0  min5    FALSE 0.8597970
## 6    155     Inf      1  max1    FALSE 1.1011285


i3.gflow.cx.res <- create.gflow.cx(graph$adj_list,
                                graph$weight_list,
                                i2.gflow.cx.res$harmonic_predictions,
                                hop.idx.thld = 5,
                                smoother.type = 0,  # Default: Weighted Mean
                                max.outer.iterations = 5,
                                max.inner.iterations = 100,
                                smoothing.tolerance = 1e-6,
                                sigma = 1.0,
                                process.in.order = TRUE,
                                verbose = TRUE,
                                detailed.recording = TRUE)

i3.he.lextr.nbhds <- compute.extrema.hop.nbhds(graph$adj_list,
                                               graph$weight_list,
                                               i3.gflow.cx.res$harmonic_predictions)
i3.he.lextr.nbhds$extrema_df
##   vertex hop_idx is_max label spurious     value
## 1    478      11      0  min1    FALSE 0.4011693
## 2    406     Inf      0  min2    FALSE 0.3257621
## 3    129      12      0  min3    FALSE 0.3927716
## 4      9       8      0  min4    FALSE 0.8597970
## 5    155     Inf      1  max1    FALSE 1.1011285

lo.errors <- c()
for (i in 4:9) {
    k <- k.values[i]
    cat("----------------------------\n")
    cat("Fitting k =", k, "(", i, "/", n.k, ")\n")
    cat("----------------------------\n")
    ## Start timing this iteration
    iter_start <- proc.time()
    ##
    graph <- graphs$geom_pruned_graphs[[i]]
    ## graph.diffusion.smoother() was archived in legacy graph-CE cleanup (2026-02-15).
    ## Keep slot for historical comparison plots, but skip this step in current codebase.
    ds.errors[i] <- NA_real_
    cat("  gds step skipped (archived legacy API)\n\n")
    ##
    iter2_start <- proc.time()
    lo.fit <- graph.deg0.lowess.cv(graph$adj_list,
                                   graph$weight_list,
                                   y = as.double(y),
                                   min.bw.factor = 0.05,
                                   max.bw.factor = 0.5,
                                   n.bws = 20,
                                   log.grid = TRUE,
                                   dist.normalization.factor = 3,
                                   use.uniform.weights = TRUE,
                                   kernel.type = 7L,
                                   n.folds = 10,
                                   with.bw.predictions = FALSE,
                                   verbose = FALSE)
    lo.errors[i] <- lo.fit$bw_errors[lo.fit$opt_bw_idx]
    cat("  lo time:", round((proc.time() - iter2_start)["elapsed"], 2), "sec\n\n")
    ##
    iter_elapsed <- proc.time() - iter_start
    k.timing[i] <- iter_elapsed["elapsed"]
    cat("  Training error:", round(ds.errors[i], 4), "\n")
    cat("  GCV score:",      round(lo.errors[i], 4), "\n")
    cat("  Time:", round(k.timing[i], 2), "sec\n\n")
}

## out.dir = "~/current_projects/ZB/analysis_output/asv"
## (file <- paste0(out.dir, "zambia_5k_lo_errors.rda"))
## save(lo.errors, file = file)
## "~/current_projects/ZB/analysis_output/asvzambia_5k_lo_errors.rda"

plot(lo.errors, type = "b")


i <- 4
graph <- graphs$geom_pruned_graphs[[i]]

system.time(lpf.fit <- klaps.low.pass.smoother(graph$adj_list,
                                               graph$weight_list,
                                               y = as.double(y),
                                               n.evectors.to.compute = 30,
                                               min.num.eigenvectors = 2L,
                                               tau.factor = 0.1,
                                               radius.factor = 10.0,
                                               laplacian.power = 1,
                                               n.candidates = 100,
                                               log.grid = TRUE,
                                               with.k.predictions = FALSE,
                                               verbose = TRUE))

lpf.gcv.scores[i] <- lpf.fit$gcv_scores[lpf.fit$opt_k_gcv]

## ----------------------------------------------------------------------
## klaps.low.pass.smoother
## ----------------------------------------------------------------------
n.cores <- 14
registerDoParallel(n.cores)
ptm <- proc.time()
lpf.res <- foreach (i = good.idices) %dopar% {
    cat("i:",i,"\n")
    graph <- graphs$geom_pruned_graphs[[i]]
    lpf.fit <- klaps.low.pass.smoother(graph$adj_list,
                                       graph$weight_list,
                                       y = as.double(y),
                                       n.evectors.to.compute = 30,
                                       min.num.eigenvectors = 2L,
                                       tau.factor = 0.1,
                                       radius.factor = 10.0,
                                       laplacian.power = 1,
                                       n.candidates = 100,
                                       log.grid = TRUE,
                                       with.k.predictions = FALSE,
                                       verbose = TRUE)
    lpf.fit
}
elapsed.time(ptm)
stopImplicitCluster()
## DONE (1:12)

lpf.gcv.scores <- c()
for (i in seq(good.idices)) {
    lpf.fit <- lpf.res[[i]]
    lpf.gcv.scores[i] <- lpf.fit$gcv_scores[lpf.fit$opt_k_gcv]
}

## out.dir = "~/current_projects/ZB/analysis_output/asv"
## (file <- paste0(out.dir, "zambia_5k_ds_errors.rda"))
## save(lpf.fits, lpf.gcv.scores, file = file)
load("~/current_projects/ZB/analysis_output/asvzambia_5k_ds_errors.rda")

plot(lpf.gcv.scores, type = "b")

# Step 2b: OPTIONAL - Subsample for timing estimation
## ----------------------------------------------------
# Set to NULL to use full dataset, or specify sample size (e.g., 500, 600, 800)
n_subsample <- 1000  # Change to 500, 600, etc. for subsampling

if (!is.null(n_subsample) && n_subsample < nrow(data_subset$rel)) {
  cat("\n--- SUBSAMPLING FOR TIMING ESTIMATION ---\n")
  cat("Original n samples:", nrow(data_subset$rel), "\n")
  cat("Subsampling to:", n_subsample, "samples\n")

  # Set seed for reproducibility
  set.seed(123)

  # Stratified sampling to maintain sPTB/TB ratio
  sptb_idx <- which(data_subset$map$sptb.binary_int == 1)
  tb_idx <- which(data_subset$map$sptb.binary_int == 0)

  # Proportional sampling
  n_sptb_sub <- round(n_subsample * length(sptb_idx) / nrow(data_subset$rel))
  n_tb_sub <- n_subsample - n_sptb_sub

  subsample_idx <- c(
    sample(sptb_idx, min(n_sptb_sub, length(sptb_idx))),
    sample(tb_idx, min(n_tb_sub, length(tb_idx)))
  )

  # Apply subsampling
  data_subset$rel <- data_subset$rel[subsample_idx, , drop = FALSE]
  data_subset$map <- data_subset$map[subsample_idx, , drop = FALSE]
  data_subset$tag <- paste0(data_subset$tag, "_n", n_subsample)

  cat("Subsampled dataset:\n")
  cat("  n samples:", nrow(data_subset$rel), "\n")
  cat("  n sPTB:", sum(data_subset$map$sptb.binary_int == 1), "\n")
  cat("  n TB:", sum(data_subset$map$sptb.binary_int == 0), "\n")
  cat("  Proportion sPTB:", round(mean(data_subset$map$sptb.binary_int), 3), "\n")
  cat("-----------------------------------------\n\n")
}
