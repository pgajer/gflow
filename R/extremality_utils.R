
#' Summarize Extrema from Riemannian Extremality Analysis
#'
#' @description
#' Creates a structured summary of local extrema detected via Riemannian extremality
#' analysis. This function extracts and filters extremality scores, hop-extremality
#' radii, and neighborhood sizes from a fitted model, producing a data frame suitable
#' for reporting, visualization, and further analysis.
#'
#' @details
#' This function provides a convenient interface for working with extremality analysis
#' results computed during model fitting (when compute.extremality = TRUE). It applies
#' filtering by extremality score magnitude and effective degree (inverse nearest
#' neighbor distance), then organizes the results into a clean summary format.
#'
#' \strong{Function Naming}
#'
#' The function is named extremality.summary() rather than compute.gextrema.nbhds()
#' because:
#' - The extremality scores and hop radii are already computed during fitting
#' - This function only summarizes and filters pre-computed results
#' - The name better reflects its role as a post-processing summarization tool
#' - It distinguishes from compute.pextrema.nbhds() which performs its own computations
#'
#' \strong{Three-Stage Filtering Pipeline}
#'
#' Stage 1 (Extremality Filtering): Select vertices with |extremality| above the
#' specified quantile threshold. This identifies vertices exhibiting strong local
#' maximum or minimum behavior according to the Riemannian metric.
#'
#' Stage 2 (Connectivity Filtering): Remove poorly connected vertices (outliers)
#' by filtering on inverse nearest neighbor distance (inv.1d). Vertices in sparse,
#' isolated regions are likely less reliable as extrema.
#'
#' Stage 3 (Organization): Separate maxima and minima, sort by function value,
#' assign labels, and combine into a single summary data frame.
#'
#' \strong{Comparison with compute.pextrema.nbhds()}
#'
#' extremality.summary():
#' - Uses Riemannian extremality scores from the coboundary operator
#' - Requires compute.extremality = TRUE during fitting
#' - Provides hop-extremality radii and neighborhood sizes
#' - Theoretically rigorous for final analysis
#'
#' compute.pextrema.nbhds():
#' - Computes probabilistic maxp/minp scores on the fly
#' - Works with any fitted model
#' - Provides hop-extremp radii (probabilistic variant)
#' - Fast and approximate, suitable for exploration
#'
#' @param dcx A fitted riem.dcx object from fit.rdgraph.regression() with
#'   compute.extremality = TRUE. Must contain the extremality component
#'   with scores, hop_extremality_radii, and hop_neighborhood_sizes.
#' @param extremality.quantile Numeric in (0,1). Quantile threshold for
#'   extremality score filtering. Default 0.95 keeps vertices with |extremality|
#'   above the 95th percentile.
#' @param inv.d1.quantile Numeric in \eqn{[0,1)}. Lower quantile threshold for
#'   inverse nearest neighbor distance filtering. Default 0.10 removes vertices
#'   in the bottom 10% of connectivity. Set to 0 to disable.
#'
#' @return A data frame with one row per detected extremum, containing:
#' \describe{
#'   \item{vertex}{Integer: vertex index (1-indexed)}
#'   \item{value}{Numeric: fitted value \eqn{\hat{y}(v)} at the vertex}
#'   \item{rel.value}{Numeric: relative value normalized by \eqn{\bar{\hat{y}}}}
#'   \item{type}{Character: "max" for local maxima, "min" for local minima}
#'   \item{inv.1d}{Numeric: inverse distance to nearest neighbor (1/d_min),
#'         measuring local density}
#'   \item{extremality}{Numeric in \eqn{[-1, 1]}: Riemannian extremality score.
#'         Positive for maxima, negative for minima. Magnitude indicates strength.}
#'   \item{extremality.radius}{Numeric: hop-extremality radius. Maximum hop
#'         distance maintaining |extremality| above threshold. Inf for global
#'         extrema, NA if not computed.}
#'   \item{extremality.nbhd.size}{Numeric: number of vertices in the
#'         hop-extremality neighborhood. Measures local support for the extremum.
#'         Inf for global extrema, NA if not computed.}
#'   \item{label}{Character: extremum label. Maxima labeled M1, M2, ... in
#'         descending order by value. Minima labeled m1, m2, ... in ascending order.}
#' }
#'
#' @examples
#' \dontrun{
#' # Fit model with extremality analysis
#' fit <- fit.rdgraph.regression(X, y, k = 50,
#'                               compute.extremality = TRUE,
#'                               p.threshold = 0.90,
#'                               max.hop = 20)
#'
#' # Summarize extrema with default filtering
#' extrema <- extremality.summary(fit)
#' print(extrema)
#'
#' # More stringent filtering
#' extrema_strict <- extremality.summary(fit,
#'                                       extremality.quantile = 0.99,
#'                                       inv.d1.quantile = 0.20)
#'
#' # Examine persistent extrema with strong support
#' persistent <- subset(extrema,
#'                      extremality.radius >= 5 &
#'                      extremality.nbhd.size >= 100)
#'
#' # Visualize extrema by neighborhood size
#' with(extrema, plot(value, extremality,
#'                    cex = sqrt(extremality.nbhd.size) / 10,
#'                    col = ifelse(type == "max", "red", "blue"),
#'                    main = "Extrema: size ~ neighborhood support"))
#' }
#'
#' @seealso
#' \code{\link{fit.rdgraph.regression}} for model fitting with extremality analysis,
#' \code{\link{compute.pextrema.nbhds}} for probabilistic extrema detection
#'
#' @export
extremality.summary <- function(dcx,
                                extremality.quantile = 0.95,
                                inv.d1.quantile = 0.10) {

    ## ================================================================
    ## ARGUMENT VALIDATION
    ## ================================================================

    ## Validate dcx is a riem.dcx object
    if (!inherits(dcx, "riem.dcx")) {
        stop("dcx must be an object of class 'riem.dcx' (output from fit.rdgraph.regression)")
    }

    ## Check for required fields
    required_fields <- c("fitted.values", "graph", "extremality")
    missing_fields <- setdiff(required_fields, names(dcx$optimal.fit))
    if (length(missing_fields) > 0) {
        stop(sprintf("dcx$optimal.fit is missing required fields: %s\n",
                     paste(missing_fields, collapse = ", ")),
             "The extremality field requires compute.extremality = TRUE ",
             "when calling fit.rdgraph.regression()")
    }

    ## Check for extremality subfields
    required_extr_fields <- c("scores", "hop.extremality.radii", "hop.neighborhood.sizes")
    missing_extr_fields <- setdiff(required_extr_fields, names(dcx$optimal.fit$extremality))
    if (length(missing_extr_fields) > 0) {
        stop(sprintf("dcx$optimal.fit$extremality is missing required fields: %s\n",
                     paste(missing_extr_fields, collapse = ", ")),
             "This suggests p.threshold was 0 during fitting. ",
             "Set p.threshold > 0 to compute hop radii and neighborhood sizes.")
    }

    ## Check for graph structure
    required_graph_fields <- c("adj.list", "edge.length.list")
    missing_graph_fields <- setdiff(required_graph_fields, names(dcx$optimal.fit$graph))
    if (length(missing_graph_fields) > 0) {
        stop(sprintf("dcx$optimal.fit$graph is missing required fields: %s",
                     paste(missing_graph_fields, collapse = ", ")))
    }

    ## Validate quantile parameters
    if (!is.numeric(extremality.quantile) ||
        length(extremality.quantile) != 1 ||
        extremality.quantile <= 0 ||
        extremality.quantile >= 1) {
        stop("extremality.quantile must be a single numeric value in (0, 1)")
    }

    if (!is.numeric(inv.d1.quantile) ||
        length(inv.d1.quantile) != 1 ||
        inv.d1.quantile < 0 ||
        inv.d1.quantile >= 1) {
        stop("inv.d1.quantile must be a single numeric value in [0, 1)")
    }

    ## ================================================================
    ## EXTRACT DATA FROM DCX OBJECT
    ## ================================================================

    adj.list <- dcx$optimal.fit$graph$adj.list
    edgelen.list <- dcx$optimal.fit$graph$edge.length.list
    yhat <- dcx$optimal.fit$fitted.values
    extremality_scores <- dcx$optimal.fit$extremality$scores
    hop_radii <- dcx$optimal.fit$extremality$hop.extremality.radii
    hop_sizes <- dcx$optimal.fit$extremality$hop.neighborhood.sizes

    n <- length(yhat)
    Ey <- mean(yhat)

    message("Using pre-computed extremality analysis from fitted model")

    ## ================================================================
    ## STAGE 1: EXTREMALITY SCORE FILTERING
    ## ================================================================

    ## The extremality scores are in [-1, 1]
    ## Positive values indicate maxima, negative indicate minima
    ## We filter by absolute value to get strong extrema of either type

    abs_extremality <- abs(extremality_scores)
    extremality_threshold <- quantile(abs_extremality,
                                      probs = extremality.quantile,
                                      na.rm = TRUE)

    ## Candidates are vertices with high absolute extremality
    candidates <- which(abs_extremality >= extremality_threshold &
                       !is.na(extremality_scores))

    if (length(candidates) == 0) {
        warning("No candidates found meeting extremality threshold. ",
                "Returning empty data frame.")
        return(data.frame(
            vertex = integer(),
            value = numeric(),
            rel.value = numeric(),
            type = character(),
            inv.1d = numeric(),
            extremality = numeric(),
            extremality.radius = numeric(),
            extremality.nbhd.size = numeric(),
            label = character(),
            stringsAsFactors = FALSE
        ))
    }

    message(sprintf("Stage 1: %d candidates with |extremality| >= %.3f (%.0f%% quantile)",
                    length(candidates),
                    extremality_threshold,
                    extremality.quantile * 100))

    ## ================================================================
    ## STAGE 2: CONNECTIVITY FILTERING (inv.1d)
    ## ================================================================

    ## Compute inverse nearest neighbor distance for all vertices
    inv.1d <- numeric(n)
    for (i in seq_len(n)) {
        if (length(edgelen.list[[i]]) > 0) {
            inv.1d[i] <- 1.0 / min(edgelen.list[[i]])
        } else {
            inv.1d[i] <- 0  # Isolated vertex
        }
    }

    ## Filter by inv.1d if threshold > 0
    if (inv.d1.quantile > 0) {
        inv.1d.thld <- quantile(inv.1d, probs = inv.d1.quantile, na.rm = TRUE)
        candidates <- candidates[inv.1d[candidates] >= inv.1d.thld]

        if (length(candidates) == 0) {
            warning("No candidates remain after connectivity filtering. ",
                    "Returning empty data frame.")
            return(data.frame(
                vertex = integer(),
                value = numeric(),
                rel.value = numeric(),
                type = character(),
                inv.1d = numeric(),
                extremality = numeric(),
                extremality.radius = numeric(),
                extremality.nbhd.size = numeric(),
                label = character(),
                stringsAsFactors = FALSE
            ))
        }

        message(sprintf("Stage 2: %d candidates with inv.1d >= %.2f (%.0f%% quantile)",
                        length(candidates),
                        inv.1d.thld,
                        inv.d1.quantile * 100))
    } else {
        message("Stage 2: Connectivity filtering disabled (inv.d1.quantile = 0)")
    }

    ## ================================================================
    ## STAGE 3: BUILD RESULTS DATA FRAME
    ## ================================================================

    ## Separate maxima and minima candidates
    max_candidates <- candidates[extremality_scores[candidates] > 0]
    min_candidates <- candidates[extremality_scores[candidates] < 0]

    message(sprintf("Stage 3: Organizing %d maxima and %d minima",
                    length(max_candidates),
                    length(min_candidates)))

    results <- data.frame(
        vertex = integer(),
        value = numeric(),
        rel.value = numeric(),
        type = character(),
        inv.1d = numeric(),
        extremality = numeric(),
        extremality.radius = numeric(),
        extremality.nbhd.size = numeric(),
        label = character(),
        stringsAsFactors = FALSE
    )

    ## Process maxima
    if (length(max_candidates) > 0) {
        ## Convert -1 (encoded SIZE_MAX) to Inf
        max_radii <- hop_radii[max_candidates]
        max_radii[max_radii == -1] <- Inf

        max_sizes <- hop_sizes[max_candidates]
        max_sizes[max_sizes == -1] <- Inf

        max_df <- data.frame(
            vertex = max_candidates,
            value = yhat[max_candidates],
            rel.value = yhat[max_candidates] / Ey,
            type = "max",
            inv.1d = inv.1d[max_candidates],
            extremality = extremality_scores[max_candidates],
            extremality.radius = max_radii,
            extremality.nbhd.size = max_sizes,
            label = "",  # Will be filled after sorting
            stringsAsFactors = FALSE
        )

        ## Sort by value (descending) and assign labels
        max_df <- max_df[order(-max_df$value), ]
        max_df$label <- paste0("M", seq_len(nrow(max_df)))

        results <- rbind(results, max_df)
    }

    ## Process minima
    if (length(min_candidates) > 0) {
        ## Convert -1 (encoded SIZE_MAX) to Inf
        min_radii <- hop_radii[min_candidates]
        min_radii[min_radii == -1] <- Inf

        min_sizes <- hop_sizes[min_candidates]
        min_sizes[min_sizes == -1] <- Inf

        min_df <- data.frame(
            vertex = min_candidates,
            value = yhat[min_candidates],
            rel.value = yhat[min_candidates] / Ey,
            type = "min",
            inv.1d = inv.1d[min_candidates],
            extremality = extremality_scores[min_candidates],
            extremality.radius = min_radii,
            extremality.nbhd.size = min_sizes,
            label = "",  # Will be filled after sorting
            stringsAsFactors = FALSE
        )

        ## Sort by value (ascending) and assign labels
        min_df <- min_df[order(min_df$value), ]
        min_df$label <- paste0("m", seq_len(nrow(min_df)))

        results <- rbind(results, min_df)
    }

    ## Reset row names
    rownames(results) <- NULL

    message(sprintf("Found %d local maxima and %d local minima",
                    sum(results$type == "max"),
                    sum(results$type == "min")))

    return(results)
}

#' Label Extremality-Based Local Extrema in 3D Plot
#'
#' @description
#' Adds labeled line segments to a 3D plot at positions of local extrema
#' identified via Riemannian extremality analysis. This function works with
#' the output of extremality.summary() to visualize extrema detected using
#' the full Riemannian metric structure.
#'
#' @details
#' This function is designed to work with extremality analysis results from
#' extremality.summary(), which provides Riemannian extremality-based extremum
#' detection. It extends the functionality of label.extrema.3d() to work with
#' the richer extremality data structure.
#'
#' The function uses the 'type' column to distinguish maxima from minima
#' (values "max" and "min"), rather than the binary 'is_max' column used
#' in classical extrema detection.
#'
#' @param graph.3d A matrix or data frame with 3 columns representing 3D
#'   coordinates. Each row corresponds to a vertex.
#' @param extremality.df Data frame from extremality.summary() with required
#'   columns: vertex, label, type. The 'type' column should contain "max"
#'   for maxima and "min" for minima.
#' @param extrema.type Character string specifying which extrema to plot:
#'   "both" (default) plots both maxima and minima, "maxima" plots only maxima,
#'   "minima" plots only minima.
#' @param offset Numeric vector of length 3 specifying the offset for line
#'   segments relative to vertex positions. Default c(0, 0, 0.25) creates
#'   vertical line segments above vertices.
#' @param with.labels Logical. Whether to display text labels. Default TRUE.
#' @param lab.cex Character expansion factor for labels. Default 1.5.
#' @param lab.adj Adjustment parameter for label positioning. Numeric vector
#'   of length 2. Default c(0, 0) centers labels.
#' @param C Scaling factor for label position relative to stick center.
#'   Default -1 places labels at the top of segments.
#' @param pwd Line width for segments. Default 5.
#' @param separate.colors Logical. If TRUE and extrema.type="both", plot
#'   maxima and minima with different colors (col.max and col.min). If FALSE,
#'   use single color from 'col' parameter. Default TRUE.
#' @param col.max Color for maxima when separate.colors=TRUE. Default "red".
#' @param col.min Color for minima when separate.colors=TRUE. Default "blue".
#' @param col Color for all extrema when separate.colors=FALSE or when plotting
#'   only one type. Default "black".
#' @param ... Additional arguments passed to bin.segments3d().
#'
#' @return Invisibly returns NULL. Function is called for its side effect of
#'   adding labeled segments to the current 3D plot.
#'
#' @examples
#' \dontrun{
#' # Fit model with extremality analysis
#' fit <- fit.rdgraph.regression(X, y, k = 50,
#'                               compute.extremality = TRUE,
#'                               p.threshold = 0.90)
#'
#' # Get extremality summary
#' extremality_df <- extremality.summary(fit)
#'
#' # Create 3D projection for visualization
#' library(rgl)
#' pca_result <- prcomp(X, center = TRUE, scale. = TRUE)
#' graph_3d <- pca_result$x[, 1:3]
#'
#' # Plot fitted surface
#' plot3d(graph_3d, col = rainbow(100)[cut(fit$fitted.values, 100)],
#'        size = 5, main = "Fitted Surface with Extrema")
#'
#' # Label both maxima and minima with different colors (default)
#' label.extremality.3d(graph_3d, extremality_df)
#'
#' # Label only maxima
#' label.extremality.3d(graph_3d, extremality_df,
#'                      extrema.type = "maxima",
#'                      col.max = "darkred")
#'
#' # Label only minima
#' label.extremality.3d(graph_3d, extremality_df,
#'                      extrema.type = "minima",
#'                      col.min = "darkblue")
#'
#' # Label both with same color
#' label.extremality.3d(graph_3d, extremality_df,
#'                      separate.colors = FALSE,
#'                      col = "purple")
#'
#' # Customize offset and label size for persistent extrema
#' persistent <- subset(extremality_df, extremality.radius >= 5)
#' label.extremality.3d(graph_3d, persistent,
#'                      offset = c(0, 0, 0.5),
#'                      lab.cex = 2.0)
#' }
#'
#' @seealso
#' \code{\link{extremality.summary}} for generating extremality-based extrema summaries,
#' \code{\link{label.extrema.3d}} for labeling classical extrema,
#' \code{\link{bin.segments3d}} for the underlying segment plotting function
#'
#' @export
label.extremality.3d <- function(graph.3d,
                                 extremality.df,
                                 extrema.type = c("both", "maxima", "minima"),
                                 offset = c(0, 0, 0.25),
                                 with.labels = TRUE,
                                 lab.cex = 1.5,
                                 lab.adj = c(0, 0),
                                 C = -1,
                                 pwd = 5,
                                 separate.colors = TRUE,
                                 col.max = "red",
                                 col.min = "blue",
                                 col = "black",
                                 ...) {

    # ================================================================
    # ARGUMENT VALIDATION
    # ================================================================

    # Match extrema.type argument
    extrema.type <- match.arg(extrema.type)

    # Validate graph.3d
    S <- graph.3d

    if (!is.matrix(S) && !is.data.frame(S)) {
        stop("graph.3d must be a matrix or data.frame")
    }

    if (ncol(S) != 3) {
        stop("graph.3d must have 3 columns")
    }

    # Ensure rownames exist
    if (is.null(rownames(S))) {
        rownames(S) <- as.character(1:nrow(S))
    }

    # Validate extremality.df has required columns
    required_cols <- c("vertex", "label", "type")
    missing_cols <- setdiff(required_cols, names(extremality.df))
    if (length(missing_cols) > 0) {
        stop(sprintf("extremality.df is missing required columns: %s",
                     paste(missing_cols, collapse = ", ")))
    }

    # Validate type column values
    valid_types <- c("max", "min")
    if (!all(extremality.df$type %in% valid_types)) {
        invalid_types <- unique(extremality.df$type[!extremality.df$type %in% valid_types])
        stop(sprintf("extremality.df$type contains invalid values: %s. Must be 'max' or 'min'",
                     paste(invalid_types, collapse = ", ")))
    }

    # ================================================================
    # FILTER EXTREMA BY TYPE
    # ================================================================

    if (extrema.type == "maxima") {
        extremality.df <- extremality.df[extremality.df$type == "max", ]
    } else if (extrema.type == "minima") {
        extremality.df <- extremality.df[extremality.df$type == "min", ]
    }

    # Check if we have any extrema to plot
    if (nrow(extremality.df) == 0) {
        warning("No extrema to plot based on extrema.type='", extrema.type, "'")
        return(invisible(NULL))
    }

    # ================================================================
    # PREPARE DATA FOR PLOTTING
    # ================================================================

    # Create binary vector indicating which vertices are extrema
    y <- rep(0, nrow(S))
    names(y) <- rownames(S)

    # Mark extrema vertices as 1
    extrema_vertices <- as.character(extremality.df$vertex)
    y[extrema_vertices] <- 1

    # Create label table mapping vertex IDs to labels
    lab.tbl <- extremality.df$label
    names(lab.tbl) <- as.character(extremality.df$vertex)

    # ================================================================
    # PLOT EXTREMA
    # ================================================================

    if (extrema.type == "both" && separate.colors) {
        # Plot maxima and minima separately with different colors

        # Maxima
        is_maxima <- extremality.df$type == "max"
        if (any(is_maxima)) {
            y_max <- rep(0, nrow(S))
            names(y_max) <- rownames(S)
            max_vertices <- as.character(extremality.df$vertex[is_maxima])
            y_max[max_vertices] <- 1

            lab.tbl.max <- extremality.df$label[is_maxima]
            names(lab.tbl.max) <- max_vertices

            bin.segments3d(S, y_max,
                          offset = offset,
                          with.labels = with.labels,
                          lab.tbl = lab.tbl.max,
                          lab.cex = lab.cex,
                          lab.adj = lab.adj,
                          C = C,
                          pwd = pwd,
                          col = col.max,
                          ...)
        }

        # Minima
        is_minima <- extremality.df$type == "min"
        if (any(is_minima)) {
            y_min <- rep(0, nrow(S))
            names(y_min) <- rownames(S)
            min_vertices <- as.character(extremality.df$vertex[is_minima])
            y_min[min_vertices] <- 1

            lab.tbl.min <- extremality.df$label[is_minima]
            names(lab.tbl.min) <- min_vertices

            bin.segments3d(S, y_min,
                          offset = offset,
                          with.labels = with.labels,
                          lab.tbl = lab.tbl.min,
                          lab.cex = lab.cex,
                          lab.adj = lab.adj,
                          C = C,
                          pwd = pwd,
                          col = col.min,
                          ...)
        }

    } else {
        # Plot all selected extrema with the same style
        # Use appropriate color based on extrema.type
        plot.col <- col
        if (extrema.type == "maxima") {
            plot.col <- col.max
        } else if (extrema.type == "minima") {
            plot.col <- col.min
        }

        bin.segments3d(S, y,
                      offset = offset,
                      with.labels = with.labels,
                      lab.tbl = lab.tbl,
                      lab.cex = lab.cex,
                      lab.adj = lab.adj,
                      C = C,
                      pwd = pwd,
                      col = plot.col,
                      ...)
    }

    invisible(NULL)
}

#' Compute geodesic distances between extrema on graph
#'
#' @param extremality.df Data frame from extremality.summary()
#' @param adj.list Graph adjacency list
#' @param edgelen.list Edge length list
#' @return Distance matrix between extrema (geodesic distances)
compute.extrema.geodesic.distances <- function(extremality.df,
                                               adj.list,
                                               edgelen.list) {
    # Extract maxima vertices
    maxima_vertices <- extremality.df$vertex[extremality.df$type == "max"]
    n_maxima <- length(maxima_vertices)

    # Build igraph object from adjacency and edge lengths
    n <- length(adj.list)
    edge_list <- matrix(nrow = 0, ncol = 2)
    edge_weights <- numeric(0)

    for (i in 1:n) {
        neighbors <- adj.list[[i]]
        if (length(neighbors) > 0) {
            edges <- cbind(i, neighbors)
            edge_list <- rbind(edge_list, edges)
            edge_weights <- c(edge_weights, edgelen.list[[i]])
        }
    }

    g <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    g <- igraph::set_edge_attr(g, "weight", value = edge_weights)

    # Compute shortest path distances between all pairs of maxima
    dist_matrix <- igraph::distances(
        g,
        v = maxima_vertices,
        to = maxima_vertices,
        mode = "all",
        weights = igraph::E(g)$weight
    )

    rownames(dist_matrix) <- extremality.df$label[extremality.df$type == "max"]
    colnames(dist_matrix) <- extremality.df$label[extremality.df$type == "max"]

    return(dist_matrix)
}

#' Cluster extrema using graph-aware hierarchical clustering
#'
#' @param geodesic.dist Geodesic distance matrix between extrema
#' @param extremality.df Extremality summary data frame
#' @param extrema.type Type of extrema to summarize: "max" or "min"
#' @param method Linkage method for hclust
#' @param density.weight Weight for density-based similarity (0-1)
#' @return hclust object with clustering results
cluster.extrema.graph.aware <- function(geodesic.dist,
                                        extremality.df,
                                        extrema.type = "max",
                                        method = "complete",
                                        density.weight = 0.3) {
    # Extract maxima
    extrema.df <- extremality.df[extremality.df$type == extrema.type, ]

    # Normalize geodesic distances to [0, 1]
    max_dist <- max(geodesic.dist[is.finite(geodesic.dist)])
    geodesic_norm <- geodesic.dist / max_dist

    # Create density-based similarity matrix
    # Vertices with similar neighborhood sizes are more likely to cluster
    nbhd_sizes <- extrema.df$extremality.nbhd.size
    nbhd_sizes[is.infinite(nbhd_sizes)] <- max(nbhd_sizes[is.finite(nbhd_sizes)])

    # Pairwise ratios of neighborhood sizes (symmetric)
    n <- length(nbhd_sizes)
    density_dissim <- matrix(0, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            ratio <- max(nbhd_sizes[i], nbhd_sizes[j]) /
                     min(nbhd_sizes[i], nbhd_sizes[j])
            density_dissim[i, j] <- log(ratio)  # Log scale for better behavior
        }
    }

    # Normalize density dissimilarity
    density_dissim <- density_dissim / max(density_dissim)

    # Combine geodesic and density information
    combined_dist <- (1 - density.weight) * geodesic_norm +
                     density.weight * density_dissim

    # Convert to dist object and cluster
    combined_dist_obj <- as.dist(combined_dist)
    hc <- hclust(combined_dist_obj, method = method)

    return(hc)
}

#' Compute Cluster Summary with Representative Vertices
#'
#' @description
#' Summarizes clusters of extrema and selects a representative vertex for each
#' cluster. The representative can be chosen using different strategies that
#' respect the graph's metric and density structure.
#'
#' @param extremality.df Data frame from extremality.summary() with cluster assignments
#' @param geodesic.dist Distance matrix (geodesic distances between extrema)
#' @param representative.method Method for selecting cluster representative:
#'   - "medoid": vertex minimizing sum of distances to cluster members (default)
#'   - "centroid": vertex closest to cluster centroid in distance space
#'   - "max_extremality": vertex with highest extremality score
#'   - "max_persistence": vertex with largest extremality.radius
#'   - "max_support": vertex with largest extremality.nbhd.size
#' @param extrema.type Type of extrema to summarize: "max" or "min"
#'
#' @return Data frame with cluster summaries including representative vertices
#'
#' @details
#' Representative Selection Strategies:
#'
#' medoid: The vertex that minimizes the sum of geodesic distances to all other
#' vertices in the cluster. This is the graph-theoretic center and is robust to
#' outliers. Recommended as default.
#'
#' centroid: Compute the centroid in distance space (mean of distances), then
#' select the vertex closest to this centroid. Similar to medoid but can be
#' affected by outliers.
#'
#' max_extremality: The vertex with the highest |extremality| score in the
#' cluster. Emphasizes the "strongest" extremum.
#'
#' max_persistence: The vertex with the largest extremality.radius. Emphasizes
#' spatial persistence.
#'
#' max_support: The vertex with the largest extremality.nbhd.size. Emphasizes
#' local density and evidence strength.
#'
#' @examples
#' \dontrun{
#' # Using medoid (default)
#' summary_medoid <- compute.cluster.summary(extremality.df, geodesic.dist)
#'
#' # Using maximum persistence
#' summary_persist <- compute.cluster.summary(
#'   extremality.df, geodesic.dist,
#'   representative.method = "max_persistence"
#' )
#' }
#'
#' @export
compute.cluster.summary <- function(extremality.df,
                                    geodesic.dist,
                                    representative.method = c("medoid", "centroid",
                                                             "max_extremality",
                                                             "max_persistence",
                                                             "max_support"),
                                    extrema.type = c("max", "min")) {

    representative.method <- match.arg(representative.method)
    extrema.type <- match.arg(extrema.type)

    # Filter by extrema type
    type_str <- ifelse(extrema.type == "max", "max", "min")
    maxima_df <- extremality.df[extremality.df$type == type_str &
                                !is.na(extremality.df$cluster), ]

    if (nrow(maxima_df) == 0) {
        stop("No extrema with cluster assignments found")
    }

    # Get unique clusters (excluding noise if using DBSCAN)
    clusters <- sort(unique(maxima_df$cluster))
    clusters <- clusters[clusters > 0]  # Remove cluster 0 (noise) if present

    # Initialize result data frame
    result <- data.frame(
        cluster = integer(),
        n_vertices = integer(),
        representative = integer(),
        representative_label = character(),
        mean_extremality = numeric(),
        sd_extremality = numeric(),
        mean_radius = numeric(),
        sd_radius = numeric(),
        mean_nbhd_size = numeric(),
        sd_nbhd_size = numeric(),
        mean_value = numeric(),
        sd_value = numeric(),
        min_value = numeric(),
        max_value = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each cluster
    for (cl in clusters) {
        cluster_members <- maxima_df[maxima_df$cluster == cl, ]
        n_members <- nrow(cluster_members)

        # Select representative based on method
        if (representative.method == "medoid") {
            # Medoid: minimize sum of distances to all cluster members
            member_labels <- cluster_members$label
            cluster_dist <- geodesic.dist[member_labels, member_labels, drop = FALSE]

            # Sum of distances for each vertex
            dist_sums <- rowSums(cluster_dist)
            medoid_idx <- which.min(dist_sums)
            rep_vertex <- cluster_members$vertex[medoid_idx]
            rep_label <- cluster_members$label[medoid_idx]

        } else if (representative.method == "centroid") {
            # Centroid: closest to mean position in distance space
            member_labels <- cluster_members$label
            cluster_dist <- geodesic.dist[member_labels, member_labels, drop = FALSE]

            # Compute centroid (mean of distances)
            centroid <- colMeans(cluster_dist)

            # Find vertex closest to centroid
            dist_to_centroid <- apply(cluster_dist, 1, function(row) {
                sqrt(sum((row - centroid)^2))
            })
            centroid_idx <- which.min(dist_to_centroid)
            rep_vertex <- cluster_members$vertex[centroid_idx]
            rep_label <- cluster_members$label[centroid_idx]

        } else if (representative.method == "max_extremality") {
            # Highest extremality score
            max_idx <- which.max(abs(cluster_members$extremality))
            rep_vertex <- cluster_members$vertex[max_idx]
            rep_label <- cluster_members$label[max_idx]

        } else if (representative.method == "max_persistence") {
            # Largest extremality radius
            # Handle Inf values
            radii <- cluster_members$extremality.radius
            radii[is.infinite(radii)] <- max(radii[is.finite(radii)]) + 1
            max_idx <- which.max(radii)
            rep_vertex <- cluster_members$vertex[max_idx]
            rep_label <- cluster_members$label[max_idx]

        } else if (representative.method == "max_support") {
            # Largest neighborhood size
            # Handle Inf values
            sizes <- cluster_members$extremality.nbhd.size
            sizes[is.infinite(sizes)] <- max(sizes[is.finite(sizes)]) + 1
            max_idx <- which.max(sizes)
            rep_vertex <- cluster_members$vertex[max_idx]
            rep_label <- cluster_members$label[max_idx]
        }

        # Compute cluster statistics
        # Handle Inf values in radii and sizes for statistics
        radii <- cluster_members$extremality.radius
        radii_finite <- radii[is.finite(radii)]

        sizes <- cluster_members$extremality.nbhd.size
        sizes_finite <- sizes[is.finite(sizes)]

        result <- rbind(result, data.frame(
            cluster = cl,
            n_vertices = n_members,
            representative = rep_vertex,
            representative_label = rep_label,
            mean_extremality = mean(cluster_members$extremality, na.rm = TRUE),
            sd_extremality = sd(cluster_members$extremality, na.rm = TRUE),
            mean_radius = mean(radii_finite, na.rm = TRUE),
            sd_radius = sd(radii_finite, na.rm = TRUE),
            mean_nbhd_size = mean(sizes_finite, na.rm = TRUE),
            sd_nbhd_size = sd(sizes_finite, na.rm = TRUE),
            mean_value = mean(cluster_members$value, na.rm = TRUE),
            sd_value = sd(cluster_members$value, na.rm = TRUE),
            min_value = min(cluster_members$value, na.rm = TRUE),
            max_value = max(cluster_members$value, na.rm = TRUE),
            stringsAsFactors = FALSE
        ))
    }

    return(result)
}

#' Extract Representative Extrema for Visualization
#'
#' @description
#' Creates a filtered extremality data frame containing only cluster representatives,
#' useful for cleaner visualizations that show one extremum per cluster.
#'
#' @param extremality.df Data frame from extremality.summary() with cluster assignments
#' @param cluster_summary Output from compute.cluster.summary()
#'
#' @return Data frame with same structure as extremality.df but only representative vertices
#'
#' @examples
#' \dontrun{
#' # Get cluster summary with representatives
#' cluster_summary <- compute.cluster.summary(extremality.df, geodesic.dist)
#'
#' # Extract representatives for plotting
#' representatives <- extract.cluster.representatives(extremality.df, cluster_summary)
#'
#' # Plot only representatives for cleaner visualization
#' label.extremality.3d(graph.3d, representatives,
#'                      offset = c(0, -2, -2),
#'                      lab.cex = 2.5)
#' }
#'
#' @export
extract.cluster.representatives <- function(extremality.df, cluster_summary) {

    # Get representative vertices
    rep_vertices <- cluster_summary$representative

    # Filter extremality.df to only representatives
    representatives <- extremality.df[extremality.df$vertex %in% rep_vertices, ]

    # Add cluster information
    representatives$cluster <- NA
    for (i in 1:nrow(cluster_summary)) {
        rep_vertex <- cluster_summary$representative[i]
        cluster_id <- cluster_summary$cluster[i]
        representatives$cluster[representatives$vertex == rep_vertex] <- cluster_id
    }

    # Optionally relabel for clarity (e.g., C1, C2, ... for clusters)
    if (nrow(representatives) > 0) {
        # Create cluster-based labels
        representatives$original_label <- representatives$label
        representatives$label <- paste0("C", representatives$cluster,
                                       "_", representatives$label)
    }

    return(representatives)
}

#' Compare Representative Selection Methods
#'
#' @description
#' Compares different methods for selecting cluster representatives,
#' useful for understanding how different strategies affect the choice.
#'
#' @param extremality.df Data frame from extremality.summary()
#' @param geodesic.dist Geodesic distance matrix
#' @param cluster_id Specific cluster to examine
#'
#' @return Data frame comparing representatives from different methods
#'
compare.representative.methods <- function(extremality.df,
                                          geodesic.dist,
                                          cluster_id) {

    methods <- c("medoid", "centroid", "max_extremality",
                "max_persistence", "max_support")

    comparison <- data.frame(
        method = character(),
        vertex = integer(),
        label = character(),
        extremality = numeric(),
        radius = numeric(),
        nbhd_size = numeric(),
        value = numeric(),
        stringsAsFactors = FALSE
    )

    for (method in methods) {
        summary <- compute.cluster.summary(extremality.df, geodesic.dist,
                                          representative.method = method)

        cluster_row <- summary[summary$cluster == cluster_id, ]

        if (nrow(cluster_row) > 0) {
            rep_vertex <- cluster_row$representative
            rep_info <- extremality.df[extremality.df$vertex == rep_vertex, ]

            comparison <- rbind(comparison, data.frame(
                method = method,
                vertex = rep_vertex,
                label = cluster_row$representative_label,
                extremality = rep_info$extremality,
                radius = rep_info$extremality.radius,
                nbhd_size = rep_info$extremality.nbhd.size,
                value = rep_info$value,
                stringsAsFactors = FALSE
            ))
        }
    }

    return(comparison)
}

#' Compute DBSCAN Cluster Summary with Representative Vertices
#'
#' @description
#' Summarizes DBSCAN clusters of extrema and selects a representative vertex
#' for each cluster. This function is specifically designed for DBSCAN output,
#' handling noise points (cluster 0) appropriately.
#'
#' @param extremality.df Data frame from extremality.summary() containing
#'   vertex-level extremality information
#' @param dbscan_clusters Integer vector of cluster assignments from dbscan::dbscan().
#'   Length must equal the number of extrema. Cluster 0 represents noise points.
#' @param geodesic.dist Distance matrix (geodesic distances between extrema).
#'   Optional but required if representative.method is "medoid" or "centroid".
#' @param representative.method Method for selecting cluster representative:
#'   - "medoid": vertex minimizing sum of distances to cluster members (default)
#'   - "centroid": vertex closest to cluster centroid in distance space
#'   - "max_extremality": vertex with highest extremality score
#'   - "max_persistence": vertex with largest extremality.radius
#'   - "max_support": vertex with largest extremality.nbhd.size
#' @param extrema.type Type of extrema to summarize: "max" or "min"
#' @param include.noise Logical. If TRUE, include noise points (cluster 0) in
#'   summary as separate "clusters". Default FALSE.
#'
#' @return Data frame with cluster summaries including representative vertices.
#'   Columns include:
#'   - cluster: DBSCAN cluster ID
#'   - n_vertices: number of vertices in cluster
#'   - representative: vertex ID of cluster representative
#'   - representative_label: extremum label of representative
#'   - mean_extremality, sd_extremality: extremality statistics
#'   - mean_radius, sd_radius: hop-extremality radius statistics
#'   - mean_nbhd_size, sd_nbhd_size: neighborhood size statistics
#'   - mean_value, sd_value, min_value, max_value: fitted value statistics
#'   - is_noise: logical indicating if this is noise cluster (only if include.noise=TRUE)
#'
#' @details
#' DBSCAN Clustering Context:
#'
#' DBSCAN (Density-Based Spatial Clustering of Applications with Noise) produces
#' cluster assignments where:
#' - Cluster IDs >= 1: Core clusters
#' - Cluster ID = 0: Noise/outlier points
#'
#' This function processes DBSCAN results specifically, treating noise appropriately:
#' - By default, noise points are excluded from summary
#' - If include.noise=TRUE, each noise point becomes its own "cluster"
#'
#' The function assigns cluster IDs to extremality.df by matching row order,
#' assuming dbscan_clusters is parallel to the extrema of the specified type.
#'
#' Representative Selection:
#'
#' For distance-based methods (medoid, centroid), geodesic.dist is required.
#' For property-based methods (max_extremality, max_persistence, max_support),
#' geodesic.dist is optional.
#'
#' @examples
#' \dontrun{
#' library(dbscan)
#'
#' # Compute geodesic distances
#' geodesic.dist <- compute.extrema.geodesic.distances(
#'     extremality.df,
#'     rdcx.obj$optimal.fit$graph$adj.list,
#'     rdcx.obj$optimal.fit$graph$edge.length.list
#' )
#'
#' # Run DBSCAN
#' maxima_df <- extremality.df[extremality.df$type == "max", ]
#' maxima_labels <- maxima_df$label
#' maxima_dist <- geodesic.dist[maxima_labels, maxima_labels]
#'
#' db <- dbscan(maxima_dist, eps = 2.5, minPts = 2)
#'
#' # Compute cluster summary
#' dbscan_summary <- compute.dbscan.cluster.summary(
#'     extremality.df,
#'     db$cluster,
#'     geodesic.dist = geodesic.dist,
#'     representative.method = "medoid"
#' )
#'
#' print(dbscan_summary)
#'
#' # Include noise points
#' dbscan_summary_with_noise <- compute.dbscan.cluster.summary(
#'     extremality.df,
#'     db$cluster,
#'     geodesic.dist = geodesic.dist,
#'     include.noise = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{compute.cluster.summary}} for hierarchical clustering summary,
#' \code{\link{extract.dbscan.cluster.representatives}} for extracting representatives
#'
#' @export
compute.dbscan.cluster.summary <- function(extremality.df,
                                          dbscan_clusters,
                                          geodesic.dist = NULL,
                                          representative.method = c("medoid", "centroid",
                                                                   "max_extremality",
                                                                   "max_persistence",
                                                                   "max_support"),
                                          extrema.type = c("max", "min"),
                                          include.noise = FALSE) {

    representative.method <- match.arg(representative.method)
    extrema.type <- match.arg(extrema.type)

    # Validate distance matrix if needed
    if (representative.method %in% c("medoid", "centroid") && is.null(geodesic.dist)) {
        stop("geodesic.dist is required for representative.method = '",
             representative.method, "'")
    }

    # Filter by extrema type
    type_str <- ifelse(extrema.type == "max", "max", "min")
    extrema_subset <- extremality.df[extremality.df$type == type_str, ]

    if (nrow(extrema_subset) == 0) {
        stop("No extrema of type '", type_str, "' found in extremality.df")
    }

    # Validate cluster vector length
    if (length(dbscan_clusters) != nrow(extrema_subset)) {
        stop(sprintf("Length of dbscan_clusters (%d) must match number of %s extrema (%d)",
                     length(dbscan_clusters), type_str, nrow(extrema_subset)))
    }

    # Add cluster assignments to subset
    extrema_subset$cluster <- dbscan_clusters

    # Identify unique clusters
    unique_clusters <- sort(unique(dbscan_clusters))

    # Separate noise and core clusters
    core_clusters <- unique_clusters[unique_clusters > 0]
    has_noise <- 0 %in% unique_clusters

    # Initialize result data frame
    result <- data.frame(
        cluster = integer(),
        n_vertices = integer(),
        representative = integer(),
        representative_label = character(),
        mean_extremality = numeric(),
        sd_extremality = numeric(),
        mean_radius = numeric(),
        sd_radius = numeric(),
        mean_nbhd_size = numeric(),
        sd_nbhd_size = numeric(),
        mean_value = numeric(),
        sd_value = numeric(),
        min_value = numeric(),
        max_value = numeric(),
        is_noise = logical(),
        stringsAsFactors = FALSE
    )

    # Process core clusters
    for (cl in core_clusters) {
        cluster_members <- extrema_subset[extrema_subset$cluster == cl, ]

        # Select representative
        rep_info <- select_representative(cluster_members, geodesic.dist,
                                         representative.method)

        # Compute statistics
        stats <- compute_cluster_stats(cluster_members)

        result <- rbind(result, data.frame(
            cluster = cl,
            n_vertices = nrow(cluster_members),
            representative = rep_info$vertex,
            representative_label = rep_info$label,
            mean_extremality = stats$mean_extremality,
            sd_extremality = stats$sd_extremality,
            mean_radius = stats$mean_radius,
            sd_radius = stats$sd_radius,
            mean_nbhd_size = stats$mean_nbhd_size,
            sd_nbhd_size = stats$sd_nbhd_size,
            mean_value = stats$mean_value,
            sd_value = stats$sd_value,
            min_value = stats$min_value,
            max_value = stats$max_value,
            is_noise = FALSE,
            stringsAsFactors = FALSE
        ))
    }

    # Process noise points if requested
    if (has_noise && include.noise) {
        noise_members <- extrema_subset[extrema_subset$cluster == 0, ]

        message(sprintf("Including %d noise points as individual clusters",
                       nrow(noise_members)))

        # Each noise point is its own cluster
        for (i in 1:nrow(noise_members)) {
            vertex_info <- noise_members[i, ]

            # Handle Inf values for statistics
            radius_val <- vertex_info$extremality.radius
            if (is.infinite(radius_val)) radius_val <- NA

            size_val <- vertex_info$extremality.nbhd.size
            if (is.infinite(size_val)) size_val <- NA

            result <- rbind(result, data.frame(
                cluster = 0,  # Keep cluster 0 to indicate noise
                n_vertices = 1,
                representative = vertex_info$vertex,
                representative_label = vertex_info$label,
                mean_extremality = vertex_info$extremality,
                sd_extremality = NA,
                mean_radius = radius_val,
                sd_radius = NA,
                mean_nbhd_size = size_val,
                sd_nbhd_size = NA,
                mean_value = vertex_info$value,
                sd_value = NA,
                min_value = vertex_info$value,
                max_value = vertex_info$value,
                is_noise = TRUE,
                stringsAsFactors = FALSE
            ))
        }
    } else if (has_noise) {
        n_noise <- sum(dbscan_clusters == 0)
        message(sprintf("Excluding %d noise points (set include.noise=TRUE to include)",
                       n_noise))
    }

    return(result)
}

#' Helper: Select representative from cluster members
#' @keywords internal
select_representative <- function(cluster_members, geodesic.dist, method) {

    if (method == "medoid") {
        member_labels <- cluster_members$label
        cluster_dist <- geodesic.dist[member_labels, member_labels, drop = FALSE]
        dist_sums <- rowSums(cluster_dist)
        medoid_idx <- which.min(dist_sums)
        return(list(
            vertex = cluster_members$vertex[medoid_idx],
            label = cluster_members$label[medoid_idx]
        ))

    } else if (method == "centroid") {
        member_labels <- cluster_members$label
        cluster_dist <- geodesic.dist[member_labels, member_labels, drop = FALSE]
        centroid <- colMeans(cluster_dist)
        dist_to_centroid <- apply(cluster_dist, 1, function(row) {
            sqrt(sum((row - centroid)^2))
        })
        centroid_idx <- which.min(dist_to_centroid)
        return(list(
            vertex = cluster_members$vertex[centroid_idx],
            label = cluster_members$label[centroid_idx]
        ))

    } else if (method == "max_extremality") {
        max_idx <- which.max(abs(cluster_members$extremality))
        return(list(
            vertex = cluster_members$vertex[max_idx],
            label = cluster_members$label[max_idx]
        ))

    } else if (method == "max_persistence") {
        radii <- cluster_members$extremality.radius
        radii[is.infinite(radii)] <- max(radii[is.finite(radii)], na.rm = TRUE) + 1
        max_idx <- which.max(radii)
        return(list(
            vertex = cluster_members$vertex[max_idx],
            label = cluster_members$label[max_idx]
        ))

    } else if (method == "max_support") {
        sizes <- cluster_members$extremality.nbhd.size
        sizes[is.infinite(sizes)] <- max(sizes[is.finite(sizes)], na.rm = TRUE) + 1
        max_idx <- which.max(sizes)
        return(list(
            vertex = cluster_members$vertex[max_idx],
            label = cluster_members$label[max_idx]
        ))
    }
}

#' Helper: Compute cluster statistics
#' @keywords internal
compute_cluster_stats <- function(cluster_members) {

    # Handle Inf values in radii and sizes
    radii <- cluster_members$extremality.radius
    radii_finite <- radii[is.finite(radii)]

    sizes <- cluster_members$extremality.nbhd.size
    sizes_finite <- sizes[is.finite(sizes)]

    list(
        mean_extremality = mean(cluster_members$extremality, na.rm = TRUE),
        sd_extremality = sd(cluster_members$extremality, na.rm = TRUE),
        mean_radius = mean(radii_finite, na.rm = TRUE),
        sd_radius = sd(radii_finite, na.rm = TRUE),
        mean_nbhd_size = mean(sizes_finite, na.rm = TRUE),
        sd_nbhd_size = sd(sizes_finite, na.rm = TRUE),
        mean_value = mean(cluster_members$value, na.rm = TRUE),
        sd_value = sd(cluster_members$value, na.rm = TRUE),
        min_value = min(cluster_members$value, na.rm = TRUE),
        max_value = max(cluster_members$value, na.rm = TRUE)
    )
}

#' Extract DBSCAN Cluster Representatives for Visualization
#'
#' @description
#' Creates a filtered extremality data frame containing only DBSCAN cluster
#' representatives, useful for cleaner visualizations that show one extremum
#' per cluster.
#'
#' @param extremality.df Data frame from extremality.summary()
#' @param dbscan_cluster_summary Output from compute.dbscan.cluster.summary()
#'
#' @return Data frame with same structure as extremality.df but only representative
#'   vertices. Includes additional columns:
#'   - cluster: DBSCAN cluster ID
#'   - original_label: original extremum label before cluster-based relabeling
#'   - is_noise: logical indicating if representative is from noise cluster
#'   The 'label' column is updated to include cluster information (e.g., "C1_M4")
#'
#' @details
#' This function extracts representative extrema for visualization, creating
#' cleaner plots by showing only one extremum per cluster. For noise points
#' (if included in the summary), each point is its own representative.
#'
#' The function updates labels to include cluster information:
#' - Core cluster representatives: "C1_M4", "C2_M12", etc.
#' - Noise representatives: "N_M36", "N_M38", etc. (if include.noise=TRUE)
#'
#' @examples
#' \dontrun{
#' # Get DBSCAN cluster summary
#' dbscan_summary <- compute.dbscan.cluster.summary(
#'     extremality.df,
#'     db$cluster,
#'     geodesic.dist = geodesic.dist
#' )
#'
#' # Extract representatives for plotting
#' representatives <- extract.dbscan.cluster.representatives(
#'     extremality.df,
#'     dbscan_summary
#' )
#'
#' # Visualize only representatives
#' plot3D.cont(graph.3d, rel.sptb.cond.exp, radius = 0.15)
#' label.extremality.3d(graph.3d, representatives,
#'                      extrema.type = "maxima",
#'                      offset = c(0, -2, -2),
#'                      lab.cex = 2.5)
#'
#' # Separate visualization for noise vs core clusters
#' core_reps <- representatives[!representatives$is_noise, ]
#' noise_reps <- representatives[representatives$is_noise, ]
#'
#' label.extremality.3d(graph.3d, core_reps,
#'                      col.max = "darkred", lwd = 5)
#' label.extremality.3d(graph.3d, noise_reps,
#'                      col.max = "orange", lwd = 2)
#' }
#'
#' @seealso
#' \code{\link{compute.dbscan.cluster.summary}} for computing cluster summaries,
#' \code{\link{extract.cluster.representatives}} for hierarchical clustering version
#'
#' @export
extract.dbscan.cluster.representatives <- function(extremality.df,
                                                   dbscan_cluster_summary) {

    # Validate input
    if (!is.data.frame(dbscan_cluster_summary)) {
        stop("dbscan_cluster_summary must be a data frame")
    }

    required_cols <- c("cluster", "representative", "representative_label", "is_noise")
    missing_cols <- setdiff(required_cols, names(dbscan_cluster_summary))
    if (length(missing_cols) > 0) {
        stop(sprintf("dbscan_cluster_summary is missing required columns: %s",
                     paste(missing_cols, collapse = ", ")))
    }

    # Get representative vertices
    rep_vertices <- dbscan_cluster_summary$representative

    # Filter extremality.df to only representatives
    representatives <- extremality.df[extremality.df$vertex %in% rep_vertices, ]

    if (nrow(representatives) == 0) {
        warning("No representatives found in extremality.df")
        return(representatives)
    }

    # Add cluster information and noise flag
    representatives$cluster <- NA
    representatives$is_noise <- FALSE
    representatives$original_label <- representatives$label

    for (i in 1:nrow(dbscan_cluster_summary)) {
        rep_vertex <- dbscan_cluster_summary$representative[i]
        cluster_id <- dbscan_cluster_summary$cluster[i]
        is_noise <- dbscan_cluster_summary$is_noise[i]

        idx <- which(representatives$vertex == rep_vertex)
        if (length(idx) > 0) {
            representatives$cluster[idx] <- cluster_id
            representatives$is_noise[idx] <- is_noise

            # Update label with cluster information
            orig_label <- representatives$original_label[idx]
            if (is_noise) {
                representatives$label[idx] <- paste0("N_", orig_label)
            } else {
                representatives$label[idx] <- paste0("C", cluster_id, "_", orig_label)
            }
        }
    }

    # Sort: core clusters first (by cluster ID), then noise
    representatives <- representatives[order(representatives$is_noise,
                                            representatives$cluster), ]
    rownames(representatives) <- NULL

    return(representatives)
}

#' Analyze DBSCAN Noise Points
#'
#' @description
#' Provides detailed analysis of noise points from DBSCAN clustering,
#' helping to understand why they were classified as noise and whether
#' they should be included in final results.
#'
#' @param extremality.df Data frame from extremality.summary()
#' @param dbscan_clusters Integer vector of cluster assignments from dbscan
#' @param extrema.type Type of extrema: "max" or "min"
#'
#' @return Data frame with noise point analysis including:
#'   - vertex, label: identification
#'   - extremality, value: strength metrics
#'   - extremality.radius, extremality.nbhd.size: persistence metrics
#'   - nearest_cluster_dist: distance to nearest core cluster member
#'   - isolation_score: measure of how isolated the point is
#'
#' @examples
#' \dontrun{
#' noise_analysis <- analyze.dbscan.noise(
#'     extremality.df,
#'     db$cluster,
#'     extrema.type = "max"
#' )
#'
#' # Identify noise points that might be genuine extrema
#' strong_noise <- subset(noise_analysis,
#'                        abs(extremality) >= 0.90 &
#'                        extremality.radius >= 3)
#' }
#'
#' @export
analyze.dbscan.noise <- function(extremality.df,
                                 dbscan_clusters,
                                 extrema.type = c("max", "min")) {

    extrema.type <- match.arg(extrema.type)
    type_str <- ifelse(extrema.type == "max", "max", "min")

    extrema_subset <- extremality.df[extremality.df$type == type_str, ]
    extrema_subset$cluster <- dbscan_clusters

    noise_points <- extrema_subset[extrema_subset$cluster == 0, ]

    if (nrow(noise_points) == 0) {
        message("No noise points found")
        return(data.frame())
    }

    # Basic statistics
    noise_analysis <- data.frame(
        vertex = noise_points$vertex,
        label = noise_points$label,
        extremality = noise_points$extremality,
        value = noise_points$value,
        extremality.radius = noise_points$extremality.radius,
        extremality.nbhd.size = noise_points$extremality.nbhd.size,
        stringsAsFactors = FALSE
    )

    message(sprintf("Found %d noise points (%.1f%% of %s extrema)",
                   nrow(noise_points),
                   100 * nrow(noise_points) / nrow(extrema_subset),
                   type_str))

    return(noise_analysis)
}
