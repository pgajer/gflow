#' KLaPS (Kernel Laplacian Power Smoothing) Low-Pass Smoother
#'
#' @description
#' Applies a family of low‑pass filters by projecting the vertex‑valued signal `y`
#' onto the subspace spanned by the first \(k\) Laplacian eigenvectors, for a grid of
#' candidate \(k\) values.  Three selection criteria (eigengap, GCV, spectral‑energy)
#' are computed to choose the optimal \(k\).
#'
#' @param adj.list      List of integer vectors; adjacency list (1‑based indices).
#' @param weight.list   List of numeric vectors; edge weights matching `adj.list`.
#' @param y             Numeric vector of length \(\vert V\vert\); signal on each vertex.
#' @param n.evectors.to.compute Integer; number of Laplacian eigenvectors to precompute.
#' @param min.num.eigenvectors  Integer; minimum \(k\) to test (≥1).
#' @param max.num.eigenvectors  Integer; maximum \(k\) to test (≤ `n.evectors.to.compute`).
#' @param tau_factor
#'     Positive real scaling factor that determines the kernel bandwidth \(\tau\) as a fraction of the graph diameter.
#'     The kernel bandwidth is computed as:
#'         \[
#'         \tau = \text{tau\_factor} \times \text{graph diameter}
#'         \]
#'     Smaller `tau_factor` values result in more localized smoothing; larger values make smoothing more global.
#'
#' @param radius_factor A real scaling factor of tau radius that is not less than 1. Default 10.
#'
#' @param laplacian_power
#'     Positive odd integer specifying the power to which (I - L) is raised.
#'     Higher values apply stronger smoothing by repeatedly reinforcing low-pass filtering.
#'     Typically, larger `laplacian_power` requires smaller `tau_factor` to maintain the smoothing scale.
#' @param n.candidates Integer; number of candidate \(k\) values (grid size).
#' @param log.grid     Logical; if `TRUE`, log‑spaced grid of \(k\), else linear.
#' @param energy.threshold Numeric; fraction of total spectral energy for the energy criterion (e.g. 0.9).
#' @param with.k.predictions Logical; if `TRUE`, return full matrix of reconstructions for each \(k\).
#' @param verbose      Logical; if `TRUE`, print progress messages.
#'
#' @return
#' An object of class `"klaps_low_pass_smoother"`—a list with components:
#' \describe{
#'   \item{evalues}{Numeric vector of Laplacian eigenvalues.}
#'   \item{evectors}{Matrix of Laplacian eigenvectors (columns).}
#'   \item{candidate.ks}{Integer vector of tested \(k\) values.}
#'   \item{eigengaps}{Numeric vector \(\lambda_{i+1}-\lambda_i\).}
#'   \item{gcv.scores}{Numeric vector of GCV scores per candidate \(k\).}
#'   \item{spectral.energy}{Numeric vector of cumulative energy per \(k\).}
#'   \item{opt.k.eigengap}{Index in `candidate.ks` of largest eigengap.}
#'   \item{opt.k.gcv}{Index of minimal GCV score.}
#'   \item{opt.k.spectral.energy}{Index of first \(k\) meeting `energy.threshold`.}
#'   \item{used.method}{String code of the method actually used (e.g. `"GCV"`).}
#'   \item{predictions}{Numeric vector—final low‑pass output at optimal \(k\).}
#'   \item{k.predictions}{Optional matrix of full reconstructions (if `with.k.predictions=TRUE`).}
#' }
#' @export
klaps.low.pass.smoother <- function(
                                    adj.list,
                                    weight.list,
                                    y,
                                    n.evectors.to.compute = min(10, length(adj.list)),
                                    min.num.eigenvectors = 1L,
                                    max.num.eigenvectors = n.evectors.to.compute,
                                    tau.factor = 1 / n.evectors.to.compute,
                                    radius.factor = 10.0,
                                    laplacian.power = 3,
                                    n.candidates = 20L,
                                    log.grid = FALSE,
                                    energy.threshold = 0.9,
                                    with.k.predictions = TRUE,
                                    verbose = TRUE
                                    ) {
    ## --- Input validation ---
    if (!is.list(adj.list) || !is.list(weight.list))
        stop("'adj.list' and 'weight.list' must both be lists")

    if (length(adj.list) != length(weight.list))
        stop("Must have same number of vertices in 'adj.list' and 'weight.list'")

    n.vertices <- length(adj.list)

    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]]))
            stop(sprintf("Vertex %d: lengths of adj.list and weight.list differ", i))
    }

    if (!is.numeric(y) || length(y) != n.vertices)
        stop("'y' must be numeric of length ", n.vertices)

    ## Check numeric/integer parameters
    if (!is.numeric(n.evectors.to.compute) || length(n.evectors.to.compute)!=1
        || n.evectors.to.compute < 1 || n.evectors.to.compute > n.vertices)
        stop("'n.evectors.to.compute' must be integer in [1, n.vertices]")

    if (!is.numeric(min.num.eigenvectors)||length(min.num.eigenvectors)!=1
        || min.num.eigenvectors < 1)
        stop("'min.num.eigenvectors' must be integer ≥ 1")

    if (!is.numeric(max.num.eigenvectors)||length(max.num.eigenvectors)!=1
        || max.num.eigenvectors < min.num.eigenvectors
        || max.num.eigenvectors > n.evectors.to.compute)
        stop("'max.num.eigenvectors' must satisfy min ≤ max ≤ n.evectors.to.compute")

    if (!is.numeric(tau.factor)||length(tau.factor)!=1
        || tau.factor <= 0 || tau.factor > 1)
        stop("'tau.factor' must be in (0,1]")

    if (!is.numeric(radius.factor)||length(radius.factor)!=1
        || radius.factor < 1)
        stop("'radius.factor' must be in greater or equal to 1")

    if (!is.numeric(laplacian.power) || length(laplacian.power)!=1 # || !(laplacian.power %% 2 ==1
        || laplacian.power < 1 || laplacian.power > 10)
        stop("'laplacian.power' must be integer in [1, 10]")

    if (!is.numeric(n.candidates)||length(n.candidates)!=1||n.candidates < 1)
        stop("'n.candidates' must be integer ≥ 1")

    if (!is.logical(log.grid)||length(log.grid)!=1)
        stop("'log.grid' must be a single logical")

    if (!is.numeric(energy.threshold)||length(energy.threshold)!=1
        || energy.threshold <= 0 || energy.threshold > 1)
        stop("'energy.threshold' must be in (0,1]")

    if (!is.logical(with.k.predictions)||length(with.k.predictions)!=1)
        stop("'with.k.predictions' must be a single logical")

    if (!is.logical(verbose)||length(verbose)!=1)
        stop("'verbose' must be a single logical")

    ## --- Prepare inputs for C++ (0-based indexing) ---
    adj.list.0based   <- lapply(adj.list, function(v) as.integer(v - 1L))

    ## --- Call into C++ ---
    result <- .Call("S_klaps_low_pass_smoother",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.integer(n.evectors.to.compute),
                    as.integer(min.num.eigenvectors),
                    as.integer(max.num.eigenvectors),
                    as.numeric(tau.factor),
                    as.numeric(radius.factor),
                    as.integer(laplacian.power),
                    as.integer(n.candidates),
                    as.logical(log.grid),
                    as.numeric(energy.threshold),
                    as.logical(with.k.predictions),
                    as.logical(verbose)
                    )

    class(result) <- "klaps_low_pass_smoother"
    result
}

#' Print Method for Graph Low-Pass Smoother
#'
#' @param x   Object of class \code{klaps_low_pass_smoother}
#' @param ... Ignored
#' @export
print.klaps_low_pass_smoother <- function(x, ...) {
    cat("Graph Low‑Pass Smoother\n")
    cat("========================\n")
    cat(sprintf("Number of vertices:           %d\n", length(x$predictions)))
    cat(sprintf("Number of candidates tested:  %d\n", length(x$candidate.ks)))
    cat(sprintf("Method used:                  %s\n", x$used.method))
    cat("Call summary() for detailed diagnostics.\n")
    invisible(x)
}

#' Summary Method for Graph Low-Pass Smoother
#'
#' @param object   Object of class \code{klaps_low_pass_smoother}
#' @param ...      Ignored
#' @return A summary object of class \code{summary.klaps_low_pass_smoother}
#' @export
summary.klaps_low_pass_smoother <- function(object, ...) {
    out <- list(
        n_vertices        = length(object$predictions),
        n_candidates      = length(object$candidate.ks),
        opt_k_eigengap    = object$opt_k_eigengap,
        opt_k_gcv         = object$opt_k_gcv,
        opt_k_spectral    = object$opt_k_spectral_energy,
        used_method       = object$used.method
    )
    class(out) <- "summary.klaps_low_pass_smoother"
    out
}

#' Print Method for Summary of Graph Low-Pass Smoother
#'
#' @param x   Summary object from \code{summary.klaps_low_pass_smoother}
#' @param ... Ignored
#' @export
print.summary.klaps_low_pass_smoother <- function(x, ...) {
    cat("Summary of Graph Low‑Pass Smoother\n")
    cat("----------------------------------\n")
    cat(sprintf("Vertices:                 %d\n", x$n_vertices))
    cat(sprintf("Candidates tested:        %d\n", x$n_candidates))
    cat(sprintf("Optimal k (eigengap):     %d\n", x$opt_k_eigengap))
    cat(sprintf("Optimal k (GCV):          %d\n", x$opt_k_gcv))
    cat(sprintf("Optimal k (energy):       %d\n", x$opt_k_spectral))
    cat(sprintf("Method chosen:            %s\n", x$used_method))
    invisible(x)
}
