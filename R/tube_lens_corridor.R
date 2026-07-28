#' Compute a tube-lens corridor between two vertices in a weighted graph
#'
#' @description
#' Computes the weighted shortest path \eqn{\gamma_{uv}} between two vertices,
#' forms a graph-metric tube around that path with radius equal to a fraction of
#' the path length, and intersects that tube with the endpoint lens
#' \eqn{d(u,w) \le |\gamma_{uv}|} and \eqn{d(v,w) \le |\gamma_{uv}|}.
#'
#' Optionally, the function also computes the stricter excess-filtered variant
#' that keeps only vertices satisfying
#' \eqn{d(u,w) + d(v,w) - |\gamma_{uv}| \le \tau}.
#'
#' @param adj.list A list of integer vectors describing graph adjacency with
#'   1-based vertex ids.
#' @param weight.list A list of numeric vectors of edge lengths aligned with
#'   `adj.list`.
#' @param start.vertex Integer scalar giving the first endpoint (1-based).
#' @param end.vertex Integer scalar giving the second endpoint (1-based).
#' @param path.relative.radius Non-negative scalar. The tube radius is computed
#'   as `path.relative.radius * path.length`.
#' @param excess.tol Optional non-negative scalar giving the excess tolerance
#'   \eqn{\tau}. If omitted for `mode = "excess"` or `"both"`, it defaults to
#'   the computed tube radius.
#' @param mode Character string; one of `"base"`, `"excess"`, or `"both"`.
#'
#' @return A list with components:
#'   \item{path.vertices}{Weighted shortest path between the endpoints.}
#'   \item{path.arc.length}{Normalized arc-length coordinate of each path vertex in
#'     `path.vertices`, ranging from 0 at `start.vertex` to 1 at `end.vertex`.}
#'   \item{path.length}{Length of that path.}
#'   \item{tube.radius}{Absolute tube radius used for the path-centered tube.}
#'   \item{excess.tolerance}{Excess tolerance used when the excess variant is requested.}
#'   \item{tube.vertices}{Vertices in the path-centered tube.}
#'   \item{tube.geodesic.distances}{Geodesic distance from each tube vertex to the path.}
#'   \item{corridor.vertices}{Vertices in the tube-lens corridor.}
#'   \item{t.balance}{Balanced longitudinal coordinate for each vertex in
#'     `corridor.vertices`, normalized to `[0,1]`.}
#'   \item{harmonic.t}{Harmonic longitudinal coordinate for each vertex in
#'     `corridor.vertices`, normalized to `[0,1]`.}
#'   \item{distance.to.path}{Geodesic distance from each corridor vertex to the
#'     shortest path centerline.}
#'   \item{excess}{Geodesic excess `d(u,w) + d(v,w) - path.length` for each
#'     vertex in `corridor.vertices`.}
#'   \item{excess.vertices}{Vertices in the excess-filtered corridor, or `NULL`
#'     when `mode = "base"`.}
#'   \item{selected.vertices}{The primary selected set for the requested mode.}
#'   \item{mode}{The requested mode.}
#'
#' @export
compute.tube.lens.corridor <- function(adj.list,
                                       weight.list,
                                       start.vertex,
                                       end.vertex,
                                       path.relative.radius = 0.1,
                                       excess.tol = NA_real_,
                                       mode = c("base", "excess", "both")) {
    mode <- match.arg(mode)

    if (!is.list(adj.list) || !all(vapply(adj.list, is.numeric, logical(1)))) {
        stop("adj.list must be a list of numeric vectors")
    }
    if (!is.list(weight.list) || !all(vapply(weight.list, is.numeric, logical(1)))) {
        stop("weight.list must be a list of numeric vectors")
    }
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }
    for (i in seq_along(adj.list)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Mismatch between adjacency and weight lists at vertex %d", i))
        }
    }

    n.vertices <- length(adj.list)
    if (!is.numeric(start.vertex) || length(start.vertex) != 1L ||
        start.vertex < 1L || start.vertex > n.vertices ||
        start.vertex != as.integer(start.vertex)) {
        stop(sprintf("start.vertex must be an integer between 1 and %d", n.vertices))
    }
    if (!is.numeric(end.vertex) || length(end.vertex) != 1L ||
        end.vertex < 1L || end.vertex > n.vertices ||
        end.vertex != as.integer(end.vertex)) {
        stop(sprintf("end.vertex must be an integer between 1 and %d", n.vertices))
    }
    if (!is.numeric(path.relative.radius) || length(path.relative.radius) != 1L ||
        !is.finite(path.relative.radius) || path.relative.radius < 0) {
        stop("path.relative.radius must be a finite non-negative number")
    }
    if (!is.numeric(excess.tol) || length(excess.tol) != 1L) {
        stop("excess.tol must be a single numeric value or NA")
    }
    if (is.finite(excess.tol) && excess.tol < 0) {
        stop("excess.tol must be non-negative when supplied")
    }

    adj.list.0based <- lapply(adj.list, function(neighbors) as.integer(neighbors - 1L))

    result <- .Call(
        "S_compute_tube_lens_corridor",
        adj.list.0based,
        weight.list,
        as.integer(start.vertex - 1L),
        as.integer(end.vertex - 1L),
        as.numeric(path.relative.radius),
        as.numeric(excess.tol),
        mode != "base"
    )

    result$selected.vertices <- switch(
        mode,
        base = result$corridor.vertices,
        excess = result$excess.vertices,
        both = result$corridor.vertices
    )
    result$mode <- mode
    result
}
