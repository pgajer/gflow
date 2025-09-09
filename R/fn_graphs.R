#' Construct a Function-Aware Graph
#'
#' Creates a new graph with edge weights modified based on function values
#'
#' @param adj.list A list where each element contains the indices of vertices adjacent to vertex i.
#'                The list should be 1-indexed (as is standard in R).
#' @param weight.list A list of the same structure as adj.list, where each element contains the
#'                   weights of edges connecting vertex i to its adjacent vertices.
#' @param function.values A numeric vector containing function values at each vertex of the graph.
#' @param weight.type Type of weight modification:
#'   0: Inverse relationship: w_new = w_old / (|f(i) - f(j)| + epsilon)
#'   1: Direct relationship: w_new = w_old * |f(i) - f(j)|
#'   2: Exponential decay: w_new = w_old * exp(-lambda * |f(i) - f(j)|)
#'   3: Power law: w_new = w_old * |f(i) - f(j)|^(-alpha)
#'   4: Sigmoid: w_new = w_old * 1/(1 + exp(beta * (|f(i) - f(j)| - tau)))
#'   5: L_p embedding: w_new = (w_old^p + alpha*|f(i) - f(j)|^q)^(1/r)
#' @param epsilon Small constant to avoid division by zero (for weight.type 0)
#' @param lambda Decay rate parameter (for weight.type 2)
#' @param alpha Power law exponent (for weight.type 3) or scaling factor (for weight.type 5)
#' @param beta Sigmoid steepness parameter (for weight.type 4)
#' @param tau Sigmoid threshold parameter (for weight.type 4)
#' @param p Power for feature distance term (for weight.type 5)
#' @param q Power for function difference term (for weight.type 5)
#' @param r Power for the overall normalization (for weight.type 5)
#' @param normalize Whether to normalize weights after modification
#' @param weight.thld Threshold for pruning edges; edges with weight > weight.thld will be pruned.
#'                   Set to a negative value to disable pruning.
#'
#' @return A new graph object with modified weights
#' @export
construct.function.aware.graph <- function(adj.list,
                                           weight.list,
                                           function.values,
                                           weight.type = 0,
                                           epsilon = 1e-6,
                                           lambda = 1.0,
                                           alpha = 1.0,
                                           beta = 5.0,
                                           tau = 0.1,
                                           p = 2.0,
                                           q = 2.0,
                                           r = 2.0,
                                           normalize = FALSE,
                                           weight.thld = -1.0) {
  ## Validate inputs
  if (length(function.values) != length(adj.list)) {
    stop("Length of function.values must match the number of vertices in the graph")
  }

  adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

  ## Call C++ function
  result <- .Call("S_construct_function_aware_graph",
                  adj.list.0based,
                  weight.list,
                  function.values,
                  as.integer(weight.type),
                  as.numeric(epsilon),
                  as.numeric(lambda),
                  as.numeric(alpha),
                  as.numeric(beta),
                  as.numeric(tau),
                  as.numeric(p),
                  as.numeric(q),
                  as.numeric(r),
                  as.logical(normalize),
                  as.numeric(weight.thld))

  return(result)
}

#' Analyze Function-Aware Weight Modifications
#'
#' Analyzes how different weighting schemes would modify edge weights based on function values.
#'
#' @param adj.list A list where each element contains the indices of vertices adjacent to vertex i.
#'                The list should be 1-indexed (as is standard in R).
#' @param weight.list A list of the same structure as adj.list, where each element contains the
#'                   weights of edges connecting vertex i to its adjacent vertices.
#' @param function.values A numeric vector containing function values at each vertex of the graph.
#' @param weight.types A vector of weight.type values to analyze (see construct.function.aware.graph
#'                    for details on weight types)
#' @param epsilon Small constant to avoid division by zero (for weight.type 0)
#' @param lambda Decay rate parameter (for weight.type 2)
#' @param alpha Power law exponent (for weight.type 3) or scaling factor (for weight.type 5)
#' @param beta Sigmoid steepness parameter (for weight.type 4)
#' @param tau Sigmoid threshold parameter (for weight.type 4)
#' @param p Power for feature distance term (for weight.type 5)
#' @param q Power for function difference term (for weight.type 5)
#' @param r Power for the overall normalization (for weight.type 5)
#'
#' @return A list where each element corresponds to a weight.type and contains a vector of
#'         modified weights for all edges in the graph
#' @export
analyze.function.aware.weights <- function(adj.list,
                                          weight.list,
                                          function.values,
                                          weight.types = 0:5,
                                          epsilon = 1e-6,
                                          lambda = 1.0,
                                          alpha = 1.0,
                                          beta = 5.0,
                                          tau = 0.1,
                                          p = 2.0,
                                          q = 2.0,
                                          r = 2.0) {
  ## Validate inputs
  if (length(function.values) != length(adj.list)) {
    stop("Length of function.values must match the number of vertices in the graph")
  }

  adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

  ## Call C++ function
  result <- .Call("S_analyze_function_aware_weights",
                  adj.list.0based,
                  weight.list,
                  function.values,
                  as.integer(weight.types),
                  as.numeric(epsilon),
                  as.numeric(lambda),
                  as.numeric(alpha),
                  as.numeric(beta),
                  as.numeric(tau),
                  as.numeric(p),
                  as.numeric(q),
                  as.numeric(r))

  return(result)
}

#' Plot Weight Modification Comparison
#'
#' Creates a visualization comparing how different weighting schemes
#' modify edge weights based on function value differences.
#'
#' @param adj.list A list where each element contains the indices of vertices adjacent to vertex i.
#' @param weight.list A list of the same structure as adj.list, with edge weights.
#' @param function.values Vector of function values at each vertex
#' @param weight.types Vector of weight types to analyze (default: 0:5)
#' @param epsilon Small constant to avoid division by zero (for weight.type 0)
#' @param lambda Decay rate parameter (for weight.type 2)
#' @param alpha Power law exponent (for weight.type 3) or scaling factor (for weight.type 5)
#' @param beta Sigmoid steepness parameter (for weight.type 4)
#' @param tau Sigmoid threshold parameter (for weight.type 4)
#' @param p Power for feature distance term (for weight.type 5)
#' @param q Power for function difference term (for weight.type 5)
#' @param r Power for the overall normalization (for weight.type 5)
#' @param log.scale Whether to use log scale for the y-axis
#'
#' @return Invisibly returns the data frame used for plotting
#' @export
function.aware.weights.plot <- function(adj.list,
                                        weight.list,
                                        function.values,
                                        weight.types = 0:5,
                                        epsilon = 1e-6,
                                        lambda = 1.0,
                                        alpha = 1.0,
                                        beta = 5.0,
                                        tau = 0.1,
                                        p = 2.0,
                                        q = 2.0,
                                        r = 2.0,
                                        log.scale = TRUE) {

    ## Get weight distributions
    weight_distributions <- analyze.function.aware.weights(
        adj.list,
        weight.list,
        function.values,
        weight.types,
        epsilon,
        lambda,
        alpha,
        beta,
        tau,
        p,
        q,
        r
    )

    ## Collect edge function differences
    edge_diffs <- c()
    edge_weights <- c()

    ## Process each vertex and its edges to collect function differences
    for (i in 1:length(adj.list)) {
        for (j_idx in seq_along(adj.list[[i]])) {
            j <- adj.list[[i]][j_idx]
            if (i < j) {  ## Process each edge only once
                edge_diffs <- c(edge_diffs, abs(function.values[i] - function.values[j]))
                edge_weights <- c(edge_weights, weight.list[[i]][j_idx])
            }
        }
    }

    ## Create data frame for plotting
    plot_data <- data.frame(
        f_diff = rep(edge_diffs, length(weight.types)),
        original_weight = rep(edge_weights, length(weight.types)),
        weight_type = factor(rep(names(weight_distributions), each = length(edge_diffs))),
        modified_weight = unlist(weight_distributions)
    )

    ## Create the plot using base R graphics
    
    ## Set up plot parameters
    oldpar <- par(mar = c(5, 4, 4, 8), xpd = TRUE)
    on.exit(par(oldpar))
    
    ## Define colors for different weight types
    weight_type_names <- unique(as.character(plot_data$weight_type))
    n_types <- length(weight_type_names)
    colors <- rainbow(n_types, alpha = 0.5)
    color_map <- setNames(colors, weight_type_names)
    
    ## Set up the plot
    if (log.scale) {
        ## Filter out non-positive values for log scale
        valid_data <- plot_data[plot_data$modified_weight > 0, ]
        
        plot(valid_data$f_diff, valid_data$modified_weight,
             log = "y",
             type = "n",
             xlab = "Function Value Difference |f(i) - f(j)|",
             ylab = "Modified Edge Weight (log scale)",
             main = "Edge Weight Modification by Function Difference")
    } else {
        plot(plot_data$f_diff, plot_data$modified_weight,
             type = "n",
             xlab = "Function Value Difference |f(i) - f(j)|",
             ylab = "Modified Edge Weight",
             main = "Edge Weight Modification by Function Difference")
    }
    
    ## Plot points for each weight type
    for (wt in weight_type_names) {
        subset_data <- plot_data[plot_data$weight_type == wt, ]
        if (log.scale) {
            subset_data <- subset_data[subset_data$modified_weight > 0, ]
        }
        points(subset_data$f_diff, subset_data$modified_weight, 
               col = color_map[wt], pch = 19, cex = 0.5)
        
        ## Add smoothed trend line using lowess
        if (nrow(subset_data) > 3) {
            smooth_fit <- lowess(subset_data$f_diff, subset_data$modified_weight)
            lines(smooth_fit, col = substr(color_map[wt], 1, 7), lwd = 2)
        }
    }
    
    ## Add legend
    legend("topright", 
           inset = c(-0.3, 0),
           legend = weight_type_names,
           col = substr(color_map, 1, 7),
           lty = 1,
           lwd = 2,
           pch = 19,
           pt.cex = 0.8,
           title = "Weight Type",
           bty = "n")
    
    invisible(plot_data)
}
