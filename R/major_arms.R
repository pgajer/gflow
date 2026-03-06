#' Detect Major Arms in a Weighted Graph
#'
#' Identifies major arms using a geometry-only strategy:
#' core-by-eccentricity, distance-to-core, peak persistence, and branch mass.
#'
#' @param adj.list Graph adjacency list (1-based vertex indices).
#' @param weight.list Edge-length list aligned with `adj.list`.
#' @param core.quantile Numeric in (0,1) defining the low-eccentricity core.
#' @param use.approx.eccentricity Logical; if `TRUE`, uses landmark approximation.
#' @param n.landmarks Integer >= 1; number of landmarks for approximation.
#' @param min.arm.size Integer >= 1; minimum number of vertices assigned to a tip.
#' @param min.persistence.quantile Numeric in [0,1]; quantile threshold over tip persistence.
#' @param min.length.quantile Numeric in [0,1]; quantile threshold over tip length.
#' @param max.arms Optional positive integer cap on selected arms.
#' @param seed Integer seed for landmark initialization.
#' @param verbose Logical; print C++ progress.
#'
#' @return A list of class `major_arms` with selected arm tips and diagnostics.
#' Major components include:
#' \describe{
#'   \item{major.arms}{Selected arm tips (1-based).}
#'   \item{summary}{Per-vertex diagnostics table (1-based indices).}
#'   \item{arm.table}{Per-tip diagnostics table and selection flags.}
#' }
#'
#' @export
detect.major.arms <- function(adj.list,
                              weight.list,
                              core.quantile = 0.10,
                              use.approx.eccentricity = TRUE,
                              n.landmarks = 64L,
                              min.arm.size = 50L,
                              min.persistence.quantile = 0.80,
                              min.length.quantile = 0.80,
                              max.arms = NULL,
                              seed = 1L,
                              verbose = FALSE) {
    if (!is.list(adj.list)) stop("'adj.list' must be a list.")
    if (!is.list(weight.list)) stop("'weight.list' must be a list.")
    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length.")
    }
    if (!is.numeric(core.quantile) || length(core.quantile) != 1L ||
        !is.finite(core.quantile) || core.quantile <= 0 || core.quantile >= 1) {
        stop("'core.quantile' must be a finite scalar in (0, 1).")
    }
    if (!is.numeric(n.landmarks) || length(n.landmarks) != 1L || !is.finite(n.landmarks) || n.landmarks < 1) {
        stop("'n.landmarks' must be a finite scalar >= 1.")
    }
    if (!is.numeric(min.arm.size) || length(min.arm.size) != 1L || !is.finite(min.arm.size) || min.arm.size < 1) {
        stop("'min.arm.size' must be a finite scalar >= 1.")
    }
    if (!is.numeric(min.persistence.quantile) || length(min.persistence.quantile) != 1L ||
        !is.finite(min.persistence.quantile) || min.persistence.quantile < 0 || min.persistence.quantile > 1) {
        stop("'min.persistence.quantile' must be a finite scalar in [0, 1].")
    }
    if (!is.numeric(min.length.quantile) || length(min.length.quantile) != 1L ||
        !is.finite(min.length.quantile) || min.length.quantile < 0 || min.length.quantile > 1) {
        stop("'min.length.quantile' must be a finite scalar in [0, 1].")
    }
    if (!is.null(max.arms)) {
        if (!is.numeric(max.arms) || length(max.arms) != 1L || !is.finite(max.arms) || max.arms < 1) {
            stop("'max.arms' must be NULL or a finite scalar >= 1.")
        }
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
        stop("'seed' must be a finite scalar.")
    }
    if (!is.logical(use.approx.eccentricity) || length(use.approx.eccentricity) != 1L) {
        stop("'use.approx.eccentricity' must be a scalar logical.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("'verbose' must be a scalar logical.")
    }

    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1L))
    max.arms.int <- if (is.null(max.arms)) 0L else as.integer(max.arms)

    res <- .Call(
        "S_detect_major_arms",
        adj.list.0,
        weight.list,
        as.double(core.quantile),
        as.logical(use.approx.eccentricity),
        as.integer(n.landmarks),
        as.integer(min.arm.size),
        as.double(min.persistence.quantile),
        as.double(min.length.quantile),
        as.integer(max.arms.int),
        as.integer(seed),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    idx.fields <- c("major_arms", "local_maxima", "core_vertices")
    for (nm in idx.fields) {
        if (!is.null(res[[nm]])) res[[nm]] <- as.integer(res[[nm]]) + 1L
    }
    if (!is.null(res$diagnostics$landmarks)) {
        res$diagnostics$landmarks <- as.integer(res$diagnostics$landmarks) + 1L
    }
    if (!is.null(res$diagnostics$candidate_tips)) {
        res$diagnostics$candidate_tips <- as.integer(res$diagnostics$candidate_tips) + 1L
    }
    if (!is.null(res$assigned_tip)) {
        aa <- as.integer(res$assigned_tip)
        aa[aa >= 0] <- aa[aa >= 0] + 1L
        aa[aa < 0] <- NA_integer_
        res$assigned_tip <- aa
    }
    if (!is.null(res$summary) && is.data.frame(res$summary)) {
        if ("vertex" %in% names(res$summary)) res$summary$vertex <- as.integer(res$summary$vertex) + 1L
        if ("assigned_tip" %in% names(res$summary)) {
            aa <- as.integer(res$summary$assigned_tip)
            aa[aa >= 0] <- aa[aa >= 0] + 1L
            aa[aa < 0] <- NA_integer_
            res$summary$assigned_tip <- aa
        }
    }

    tip.vertex <- seq_along(res$tip_basin_size)
    arm.table <- data.frame(
        vertex = as.integer(tip.vertex),
        is.local.max = as.logical(tip.vertex %in% as.integer(res$local_maxima)),
        is.major.arm.tip = as.logical(res$is_tip_selected),
        arm.rank = as.integer(res$arm_rank_of_tip),
        basin.size = as.integer(res$tip_basin_size),
        persistence = as.numeric(res$tip_persistence),
        length = as.numeric(res$tip_length),
        score = as.numeric(res$tip_score),
        stringsAsFactors = FALSE
    )
    arm.table <- arm.table[arm.table$is.local.max %in% TRUE, , drop = FALSE]
    arm.table <- arm.table[order(is.na(arm.table$arm.rank), arm.table$arm.rank, -arm.table$score), , drop = FALSE]

    names(res)[names(res) == "major_arms"] <- "major.arms"
    names(res)[names(res) == "core_vertices"] <- "core.vertices"
    names(res)[names(res) == "distance_to_core"] <- "distance.to.core"
    names(res)[names(res) == "is_core"] <- "is.core"
    names(res)[names(res) == "core_threshold"] <- "core.threshold"
    names(res)[names(res) == "persistence_threshold"] <- "persistence.threshold"
    names(res)[names(res) == "length_threshold"] <- "length.threshold"

    res$arm.table <- arm.table
    class(res) <- c("major_arms", class(res))
    res
}
