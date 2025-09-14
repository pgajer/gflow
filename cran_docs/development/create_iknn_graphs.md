# Create ikNN graph notes

Great start! A few correctness nits + one or two logic mismatches to fix. Here’s what I’d change and why, followed by a cleaned-up wrapper you can paste in.

## Key fixes

1. **`.Call()` target**
   Use a *string* for the native symbol. Your current call `.Call(S_create_iknn_graphs, …)` will fail unless you passed an external pointer. With routine registration, the standard pattern is:

```r
.Call("S_create_iknn_graphs", ..., PACKAGE = "gflow")
```

(If you’re using Rcpp-generated wrappers, the symbol is usually `"_gflow_S_create_iknn_graphs"`.)

2. **Wrong variable name**
   You pass `n_cores` but your R arg is `n.cores`. That’s a bug.

```r
if (is.null(n.cores)) NULL else as.integer(n.cores)
```

3. **Don’t shift the geometric threshold by +1.0**
   You currently do:

```r
as.double(max.path.edge.ratio.deviation.thld + 1.0)
```

Your C++ expects the *deviation* threshold (e.g., 0.1), not a ratio (e.g., 1.1). Passing `+1.0` will effectively disable/alter pruning behavior.

4. **K + 1 policy needs to be consistent**
   You comment that ANN includes self, so you send `kmin+1`, `kmax+1`. That’s fine **only if**:

* the C++ labeling (e.g., `edge_stats_list` names, `k_statistics` “k” column) reflects the *original* k (without +1), or
* you document that internal “k” is `k+1`.

Right now your docs claim a `k` column, but your C++ `k_statistics` (as you pasted earlier) has **7 columns without k**. Either:

* add `k` in C++ as the first column, or
* remove “k” from the documented columns and the R-side `colnames` fill-in.

5. **Input docs vs behavior**
   Docs say “X cannot be a data frame”, but code `as.matrix(X)` happily accepts a numeric data frame. Reword to “numeric matrix or object coercible to a numeric matrix”.

6. **Ensure double storage**
   C++ requires `REALSXP`. Since you error in C++ if not double, coerce in R for a better UX:

```r
storage.mode(X) <- "double"
```

7. **Ranges & messages**

* You already enforce `0 ≤ percentile ≤ 1` and `0 ≤ deviation < 0.2`. Good.
* You check `nrow(X) > kmax`. Good. Consider also bounding `kmax` against `n-1` in the help text.

8. **PCA attributes**
   Storing `pca_projected_X` inside `attr(result,"pca")` can be large; if memory is a concern, omit it or make it optional. (Not required—just a note.)

## Drop-in revised wrapper

```r
#' Create intersection k-nearest neighbor graphs with dual pruning
#'
#' @description
#' For each \eqn{k \in [k_{\mathrm{min}},\,k_{\mathrm{max}}]}, builds an
#' intersection-weighted k-NN graph and applies two pruning schemes:
#' (1) geometric (path-to-edge ratio) and (2) intersection-size.
#' Optionally performs PCA before graph construction.
#'
#' @param X A numeric matrix (or object coercible to a numeric matrix) with
#'   rows = observations and columns = features.
#' @param kmin Integer \eqn{\ge 1}, the minimum k.
#' @param kmax Integer \eqn{> k_{\mathrm{min}}}, the maximum k.
#' @param max.path.edge.ratio.deviation.thld Numeric in \eqn{[0, 0.2)}. Geometric
#'   pruning removes an edge \eqn{(i,j)} when there exists an alternative path
#'   between \eqn{i} and \eqn{j} whose path/edge length ratio minus 1.0 is
#'   \emph{less than} this threshold.
#' @param path.edge.ratio.percentile Numeric in \eqn{[0,1]}. Only edges with
#'   length above this percentile are considered for geometric pruning.
#' @param compute.full Logical. If `TRUE`, return the pruned graphs; if `FALSE`,
#'   return only edge statistics.
#' @param pca.dim Positive integer or `NULL`. If not `NULL` and `ncol(X) > pca.dim`,
#'   PCA is used to reduce to at most `pca.dim` components.
#' @param variance.explained Numeric in \eqn{(0,1]} or `NULL`. If not `NULL`,
#'   choose the smallest number of PCs whose cumulative variance explained
#'   exceeds this threshold, capped by `pca.dim`.
#' @param n.cores Integer or `NULL`. Number of CPU cores. `NULL` uses the
#'   maximum available (OpenMP build only).
#' @param verbose Logical; print progress and timing.
#'
#' @return A list of class `"iknn_graphs"` with entries:
#' \describe{
#'   \item{k_statistics}{Matrix of per-\eqn{k} edge counts and reductions.
#'     (If the C++ side supplies column names, they’re preserved.
#'     Otherwise we add names consistent with what the C++ returns.)}
#'   \item{geom_pruned_graphs}{If `compute.full=TRUE`, list of geometrically
#'     pruned graphs (adjacency + weights); otherwise `NULL`.}
#'   \item{isize_pruned_graphs}{If `compute.full=TRUE`, list of intersection-size
#'     pruned graphs; otherwise `NULL`.}
#'   \item{edge_pruning_stats}{List (per \eqn{k}) of matrices with edge-level
#'     statistics (lengths, path/edge ratios, etc.).}
#' }
#'
#' @details
#' Geometric pruning uses the deviation threshold
#' `max.path.edge.ratio.deviation.thld` and the filtering percentile
#' `path.edge.ratio.percentile`. Intersection-size pruning currently uses
#' a maximum alternative path length of 2.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(200), 100, 2)
#' res <- create.iknn.graphs(X, kmin = 3, kmax = 10, compute.full = FALSE)
#' head(res$k_statistics)
#' }
#' @export
create.iknn.graphs <- function(X,
                               kmin,
                               kmax,
                               max.path.edge.ratio.deviation.thld = 0.1,
                               path.edge.ratio.percentile = 0.5,
                               compute.full = TRUE,
                               pca.dim = 100,
                               variance.explained = 0.99,
                               n.cores = NULL,
                               verbose = FALSE) {

  ## Coerce & basic checks
  if (!is.matrix(X)) {
    X <- try(as.matrix(X), silent = TRUE)
    if (inherits(X, "try-error"))
      stop("X must be a matrix or coercible to a numeric matrix.")
  }
  if (!is.numeric(X)) stop("X must be numeric.")
  if (any(!is.finite(X))) stop("X cannot contain NA/NaN/Inf.")
  storage.mode(X) <- "double"

  n <- nrow(X)
  if (n < 2) stop("X must have at least 2 rows (observations).")

  if (!is.numeric(kmin) || length(kmin) != 1 || kmin < 1 || kmin != floor(kmin))
    stop("kmin must be a positive integer.")
  if (!is.numeric(kmax) || length(kmax) != 1 || kmax < kmin || kmax != floor(kmax))
    stop("kmax must be an integer not smaller than kmin.")
  if (n <= kmax)
    stop("Number of observations (nrow(X)) must be greater than kmax.")

  if (!is.numeric(max.path.edge.ratio.deviation.thld) || length(max.path.edge.ratio.deviation.thld) != 1)
    stop("max.path.edge.ratio.deviation.thld must be numeric.")
  if (max.path.edge.ratio.deviation.thld < 0 || max.path.edge.ratio.deviation.thld >= 0.2)
    stop("max.path.edge.ratio.deviation.thld must be in [0, 0.2).")

  if (!is.numeric(path.edge.ratio.percentile) || length(path.edge.ratio.percentile) != 1 ||
      path.edge.ratio.percentile < 0 || path.edge.ratio.percentile > 1)
    stop("path.edge.ratio.percentile must be in [0, 1].")

  if (!is.logical(compute.full) || length(compute.full) != 1)
    stop("compute.full must be TRUE/FALSE.")
  if (!is.logical(verbose) || length(verbose) != 1)
    stop("verbose must be TRUE/FALSE.")

  if (!is.null(pca.dim)) {
    if (!is.numeric(pca.dim) || length(pca.dim) != 1 || pca.dim < 1 || pca.dim != floor(pca.dim))
      stop("pca.dim must be a positive integer or NULL.")
  }
  if (!is.null(variance.explained)) {
    if (!is.numeric(variance.explained) || length(variance.explained) != 1 ||
        variance.explained <= 0 || variance.explained > 1)
      stop("variance.explained must be in (0, 1], or NULL.")
  }

  ## PCA (optional)
  pca_info <- NULL
  if (!is.null(pca.dim) && ncol(X) > pca.dim) {
    if (verbose) message("High-dimensional data detected. Performing PCA.")
    original_dim <- ncol(X)
    if (!is.null(variance.explained)) {
      pca_analysis <- pca.optimal.components(
        X, variance.threshold = variance.explained, max.components = pca.dim
      )
      n_components <- pca_analysis$n.components
      if (verbose) {
        message(sprintf("Using %d PCs (explains %.2f%% variance)",
                        n_components, 100 * pca_analysis$variance.explained))
      }
      X <- pca.project(X, pca_analysis$pca.result, n_components)
      pca_info <- list(
        original_dim = original_dim,
        n_components = n_components,
        variance_explained = pca_analysis$variance.explained,
        cumulative_variance = pca_analysis$cumulative.variance
      )
    } else {
      if (verbose) message(sprintf("Projecting to first %d PCs", pca.dim))
      pca_result <- prcomp(X)
      X <- pca.project(X, pca_result, pca.dim)
      variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)
      pca_info <- list(
        original_dim = original_dim,
        n_components = pca.dim,
        variance_explained = variance_explained
      )
    }
  }

  ## Note on k: if your C++ assumes ANN returns self in its kNN sets,
  ## and you want to request that extra neighbor, pass k+1 here.
  ## Ensure the C++ labels/columns reflect the *original* k.
  result <- .Call("S_create_iknn_graphs",
                  X,
                  as.integer(kmin + 1L),
                  as.integer(kmax + 1L),
                  as.double(max.path.edge.ratio.deviation.thld),  # <-- no +1.0
                  as.double(path.edge.ratio.percentile),
                  as.logical(compute.full),
                  if (is.null(n.cores)) NULL else as.integer(n.cores),
                  as.logical(verbose),
                  PACKAGE = "gflow")

  ## Normalize optional matrices (placeholders may be empty)
  if (!is.null(result$birth_death_matrix)) {
    if (!is.matrix(result$birth_death_matrix) || nrow(result$birth_death_matrix) == 0) {
      result$birth_death_matrix <- matrix(
        numeric(0), nrow = 0, ncol = 4,
        dimnames = list(NULL, c("start","end","birth_time","death_time"))
      )
    } else if (is.null(colnames(result$birth_death_matrix))) {
      colnames(result$birth_death_matrix) <- c("start","end","birth_time","death_time")
    }
  }
  if (!is.null(result$double_birth_death_matrix)) {
    if (!is.matrix(result$double_birth_death_matrix) || nrow(result$double_birth_death_matrix) == 0) {
      result$double_birth_death_matrix <- matrix(
        numeric(0), nrow = 0, ncol = 4,
        dimnames = list(NULL, c("start","end","birth_time","death_time"))
      )
    }
  }

  ## Add column names to k_statistics only if missing, and match column count
  if (!is.null(result$k_statistics) && is.matrix(result$k_statistics) &&
      is.null(colnames(result$k_statistics))) {
    nc <- ncol(result$k_statistics)
    # Common layouts: with k (8 cols) or without k (7 cols)
    if (nc == 8L) {
      colnames(result$k_statistics) <- c("k",
        "n_edges",
        "n_edges_in_geom_pruned_graph",
        "n_geom_removed_edges",
        "geom_edge_reduction_ratio",
        "n_edges_in_isize_pruned_graph",
        "n_isize_removed_edges",
        "isize_edge_reduction_ratio"
      )
    } else if (nc == 7L) {
      colnames(result$k_statistics) <- c(
        "n_edges",
        "n_edges_in_geom_pruned_graph",
        "n_geom_removed_edges",
        "geom_edge_reduction_ratio",
        "n_edges_in_isize_pruned_graph",
        "n_isize_removed_edges",
        "isize_edge_reduction_ratio"
      )
    }
  }

  attr(result, "kmin") <- kmin
  attr(result, "kmax") <- kmax
  attr(result, "max_path_edge_ratio_deviation_thld") <- max.path.edge.ratio.deviation.thld
  attr(result, "path_edge_ratio_percentile") <- path.edge.ratio.percentile
  if (!is.null(pca_info)) attr(result, "pca") <- pca_info
  class(result) <- "iknn_graphs"
  result
}
```

### Final sanity checklist

* [x] `.Call` uses a string symbol and correct `PACKAGE`.
* [x] `n.cores` is passed (not the non-existent `n_cores`).
* [x] No spurious `+1.0` to the deviation threshold.
* [x] `kmin/kmax` `+1` policy is explicit; labels handled downstream.
* [x] `k_statistics` column names guarded by actual column count.
* [x] Docs match behavior (matrix or coercible; percentile and threshold ranges).




--------- original fn

#' Create Intersection k-Nearest Neighbor Graphs with Dual Pruning Methods
#'
#' @description
#' Computes a sequence of intersection-weighted k-nearest neighbor graphs for k in \eqn{[\text{kmin, \text{kmax}]}}
#' with two different edge pruning methods: geometric pruning and intersection-size pruning.
#' The function can optionally perform dimensionality reduction via PCA before graph construction.
#'
#' @param X numeric matrix where rows represent observations and columns represent features.
#'        Cannot be a data frame.
#' @param kmin integer, minimum number of nearest neighbors (>= 1)
#' @param kmax integer, maximum number of nearest neighbors (> kmin)
#' @param max.path.edge.ratio.deviation.thld numeric, threshold for geometric pruning based on path-to-edge ratio.
#'         If > 0, removes edges where the ratio of alternative path length to direct edge length
#'         minus 1.0 is less than this value. If <= 0, geometric pruning uses a default method.
#'         Must be in the interval [0, 0.2).
#' @param path.edge.ratio.percentile numeric in \eqn{[0,1]}, percentile threshold for edge lengths
#'         considered in geometric pruning. Only edges with length greater than this
#'         percentile are evaluated for path-ratio pruning.
#' @param compute.full logical, if TRUE returns all pruned graphs, if FALSE returns only
#'        edge statistics
#' @param pca.dim Maximum number of principal components to use if dimensionality reduction
#'        is applied (default: 100). Set to NULL to skip dimensionality reduction.
#' @param variance.explained Percentage of variance to be explained by the principal components
#'        (default: 0.99). If this threshold can be met with fewer components than pca.dim,
#'        the smaller number will be used. Set to NULL to use exactly pca.dim components.
#' @param n.cores The number of cores to use. Set to NULL to use maximal number of cores.
#' @param verbose Logical. If TRUE, print progress messages and timing information.
#'        Default is FALSE.
#'
#' @return A list of class "iknn_graphs" containing:
#' \describe{
#'   \item{k_statistics}{Matrix with columns: k, number of edges in original graph,
#'         number of edges after geometric pruning, number of removed edges,
#'         edge reduction ratio, number of edges after intersection-size pruning,
#'         additional edges removed in intersection-size pruning,
#'         intersection-size edge reduction ratio}
#'   \item{geom_pruned_graphs}{If compute_full=TRUE, list of geometrically pruned graphs for each k.
#'         Each graph contains adjacency lists and edge weights.
#'         If compute_full=FALSE, NULL}
#'   \item{isize_pruned_graphs}{If compute_full=TRUE, list of intersection-size pruned graphs for each k.
#'         Each graph contains adjacency lists and edge weights.
#'         If compute_full=FALSE, NULL}
#'   \item{edge_pruning_stats}{List of matrices, one for each k value, containing edge pruning statistics
#'         including edge lengths and path-to-edge length ratios}
#' }
#'
#' @details
#' The function applies two different pruning methods to construct efficient graph representations:
#'
#' 1. Geometric pruning: Based on path-to-edge ratios
#'    - Removes edges where alternative paths exist with similar or better geometric properties
#'    - Controlled by max.path.edge.ratio.deviation.thld and path.edge.ratio.percentile parameters
#'
#' 2. Intersection-size pruning: Based on intersection sizes of k-NN sets
#'    - Removes edges where alternative paths exist with all edges having larger intersection sizes
#'    - Uses a fixed maximum alternative path length of 2
#'
#' Note: The current implementation does not compute edge birth/death times. The birth_death_matrix
#' and double_birth_death_matrix fields may be present but will be empty.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(100 * 5), 100, 5)
#'
#' # Basic usage
#' result <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10,
#'   compute.full = FALSE
#' )
#'
#' # With custom pruning parameters
#' result <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10,
#'   max.path.edge.ratio.deviation.thld = 0.1,
#'   path.edge.ratio.percentile = 0.5,
#'   compute.full = TRUE,
#'   n.cores = 1,
#'   verbose = TRUE
#' )
#'
#' # View statistics for each k
#' print(result$k_statistics)
#' }
#'
#' @export
create.iknn.graphs <- function(X,
                               kmin,
                               kmax,
                               ## pruning parameters
                               max.path.edge.ratio.deviation.thld = 0.1,
                               path.edge.ratio.percentile = 0.5,
                               ## other
                               compute.full = TRUE,
                               pca.dim = 100,
                               variance.explained = 0.99,
                               n.cores = NULL,
                               verbose = FALSE) {
    ## Input validation
    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }
    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }
    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }
    n <- nrow(X)
    if (n < 2) {
        stop("X must contain at least 2 data points")
    }

    ## Check kmin
    if (!is.numeric(kmin) || length(kmin) != 1 || kmin < 1 || kmin != floor(kmin)) {
        stop("kmin must be a positive integer")
    }
    ## Check kmax
    if (!is.numeric(kmax) || length(kmax) != 1 || kmax < kmin || kmax != floor(kmax)) {
        stop("kmax must be an integer not smaller than kmin")
    }
    ## Check max.path.edge.ratio.deviation.thld
    if (!is.numeric(max.path.edge.ratio.deviation.thld) || length(max.path.edge.ratio.deviation.thld) != 1) {
        stop("max.path.edge.ratio.deviation.thld must be a numeric value")
    }

    if (max.path.edge.ratio.deviation.thld < 0 || max.path.edge.ratio.deviation.thld >= 0.2) {
        stop("max.path.edge.ratio.deviation.thld must be in the interval [0, 0.2)")
    }

    ## Check path.edge.ratio.percentile
    if (!is.numeric(path.edge.ratio.percentile) || length(path.edge.ratio.percentile) != 1 ||
        path.edge.ratio.percentile < 0.0 || path.edge.ratio.percentile > 1.0) {
        stop("path.edge.ratio.percentile must be a numeric value between 0.0 and 1.0")
    }
    ## Check compute.full
    if (!is.logical(compute.full) || length(compute.full) != 1) {
        stop("compute.full must be a logical value (TRUE/FALSE)")
    }
    ## Check verbose
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a logical value (TRUE/FALSE)")
    }
    ## Check pca.dim
    if (!is.null(pca.dim)) {
        if (!is.numeric(pca.dim) || length(pca.dim) != 1 || pca.dim < 1 || pca.dim != floor(pca.dim)) {
            stop("pca.dim must be a positive integer or NULL")
        }
    }
    ## Check variance.explained
    if (!is.null(variance.explained)) {
        if (!is.numeric(variance.explained) || length(variance.explained) != 1 ||
            variance.explained <= 0 || variance.explained > 1) {
            stop("variance.explained must be a numeric value between 0 and 1, or NULL")
        }
    }
    ## Check if number of observations is sufficient
    if (nrow(X) <= kmax) {
        stop("Number of observations must be greater than kmax")
    }

    ## PCA dimensionality reduction if needed
    pca_info <- NULL
    if (!is.null(pca.dim) && ncol(X) > pca.dim) {
        if (verbose) {
            message("High-dimensional data detected. Performing dimensionality reduction.")
        }

        # Store original dimensions for reporting
        original_dim <- ncol(X)

        # If variance.explained is specified, analyze variance
        if (!is.null(variance.explained)) {
            pca_analysis <- pca.optimal.components(X,
                                                   variance.threshold = variance.explained,
                                                   max.components = pca.dim)

            # Number of components to use (based on variance or pca.dim)
            n_components <- pca_analysis$n.components

            if (verbose) {
                message(sprintf("Using %d principal components (explains %.2f%% of variance)",
                                n_components, pca_analysis$variance.explained * 100))
            }

            # Project data onto selected components
            X <- pca.project(X, pca_analysis$pca.result, n_components)

            # Store PCA information
            pca_info <- list(
                original_dim = original_dim,
                n_components = n_components,
                variance_explained = pca_analysis$variance.explained,
                cumulative_variance = pca_analysis$cumulative.variance,
                pca_projected_X = X
            )
        } else {
            # Use fixed number of components (pca.dim)
            if (verbose) {
                message(sprintf("Projecting data onto first %d principal components", pca.dim))
            }

            # Perform PCA and projection
            pca_result <- prcomp(X)
            X <- pca.project(X, pca_result, pca.dim)

            # Calculate variance explained by pca.dim components
            variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)

            # Store PCA information
            pca_info <- list(
                original_dim = original_dim,
                n_components = pca.dim,
                variance_explained = variance_explained,
                pca_projected_X = X
            )
        }
    }

    ## Call the C++ function with updated signature
    ## Note: The C++ function expects k values incremented by 1 to account for
    ## ANN library including the reference vertex in its kNN set
    result <- .Call(S_create_iknn_graphs,
                    X,
                    as.integer(kmin + 1),
                    as.integer(kmax + 1),
                    as.double(max.path.edge.ratio.deviation.thld + 1.0),
                    as.double(path.edge.ratio.percentile),
                    as.logical(compute.full),
                    if (is.null(n.cores)) NULL else as.integer(n.cores),
                    as.logical(verbose),
                    PACKAGE = "gflow")

    ## The C++ function may return placeholder birth_death matrices
    ## We'll keep them for backward compatibility but note they may be empty
    if (!is.null(result$birth_death_matrix)) {
        if (!is.matrix(result$birth_death_matrix) || nrow(result$birth_death_matrix) == 0) {
            ## Create an empty matrix with correct structure
            result$birth_death_matrix <- matrix(numeric(0),
                                                nrow = 0,
                                                ncol = 4,
                                                dimnames = list(NULL,
                                                                c("start", "end",
                                                                  "birth_time", "death_time")))
        } else if (is.null(colnames(result$birth_death_matrix))) {
            ## Set column names if they weren't set in C++
            colnames(result$birth_death_matrix) <- c("start", "end",
                                                     "birth_time", "death_time")
        }
    }

    ## Similar handling for double_birth_death_matrix if it exists
    if (!is.null(result$double_birth_death_matrix)) {
        if (!is.matrix(result$double_birth_death_matrix) || nrow(result$double_birth_death_matrix) == 0) {
            result$double_birth_death_matrix <- matrix(numeric(0),
                                                       nrow = 0,
                                                       ncol = 4,
                                                       dimnames = list(NULL,
                                                                       c("start", "end",
                                                                         "birth_time", "death_time")))
        }
    }

    ## Add names to k_statistics columns if they weren't set in C++
    if (!is.null(result$k_statistics) && is.null(colnames(result$k_statistics))) {
        colnames(result$k_statistics) <- c("k",
                                           "n_edges",
                                           "n_edges_in_geom_pruned_graph",
                                           "n_geom_removed_edges",
                                           "geom_edge_reduction_ratio",
                                           "n_edges_in_isize_pruned_graph",
                                           "n_isize_removed_edges",
                                           "isize_edge_reduction_ratio")
    }

    ## Parameter attributes
    attr(result, "kmin") <- kmin
    attr(result, "kmax") <- kmax
    attr(result, "max_path_edge_ratio_deviation_thld") <- max.path.edge.ratio.deviation.thld
    attr(result, "path_edge_ratio_percentile") <- path.edge.ratio.percentile

    # Add PCA-related attributes if PCA was performed
    if (!is.null(pca_info)) {
        attr(result, "pca") <- pca_info
    }

    class(result) <- "iknn_graphs"

    return(result)
}
