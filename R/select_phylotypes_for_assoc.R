#' Select Phylotypes for Association Analysis Over a Vertex Set
#'
#' Filters phylotypes from a relative abundance matrix using prevalence,
#' outcome-stratified detection counts, and optional spatial coverage across a
#' distance axis (binned into quantile-based bins).
#'
#' @param phi A numeric matrix of relative abundances with rows as samples and
#'   columns as phylotypes.
#' @param outcome An integer/logical vector of outcomes of length \code{nrow(phi)}.
#'   By default, \code{1} indicates sPTB and \code{0} indicates term.
#' @param vertices Integer vector of row indices (1-based) specifying the subset
#'   of samples (vertices) used for selection.
#' @param dist Optional numeric vector of length \code{nrow(phi)} giving a
#'   distance (or pseudotime) coordinate for spatial coverage filtering. If
#'   \code{NULL}, the spatial criterion is skipped.
#' @param n.bins Integer; number of bins used to assess spatial coverage.
#'   Default is 5.
#' @param min.prevalence Minimum overall prevalence proportion required.
#'   Default is 0.20.
#' @param min.sptb.detections Minimum number of detected sPTB samples required.
#'   Default is 3.
#' @param min.term.detections Minimum number of detected term samples required.
#'   Default is 10.
#' @param min.bins Minimum number of bins with at least one detection required.
#'   Only used if \code{dist} is not \code{NULL}. Default is 3.
#' @param sort.by Column name used to sort the returned summary data frame.
#'   Default is \code{"n.prevalence"}.
#' @param decreasing Logical; if TRUE, sort descending by \code{sort.by}.
#'   Default is TRUE.
#' @param verbose Logical; if TRUE, print a short summary of sample counts and
#'   filter pass counts. Default is TRUE.
#'
#' @return An object of class \code{"phylotype_selection"}: a list with elements
#' \itemize{
#'   \item \code{selected.phylotypes}: character vector of selected phylotype names
#'   \item \code{phi.filtered}: filtered abundance matrix (subset of \code{phi})
#'   \item \code{selected.phylotypes.df}: ranked data.frame of per-phylotype metrics
#'   \item \code{filters}: logical vectors for each criterion and the combined filter
#'   \item \code{metrics}: named vectors of computed metrics (prevalence, counts, bins)
#'   \item \code{thresholds}: list of threshold parameters used
#'   \item \code{counts}: list with \code{n.samples}, \code{n.sptb}, \code{n.term}
#' }
#'
select.phylotypes.for.assoc <- function(phi,
                                        outcome,
                                        vertices,
                                        dist = NULL,
                                        n.bins = 5L,
                                        min.prevalence = 0.20,
                                        min.sptb.detections = 3L,
                                        min.term.detections = 10L,
                                        min.bins = 3L,
                                        sort.by = "n.prevalence",
                                        decreasing = TRUE,
                                        verbose = TRUE) {

    ## ------------------------------------------------------------------------
    ## Input validation
    ## ------------------------------------------------------------------------

    if (is.null(phi) || !(is.matrix(phi) || is.data.frame(phi))) {
        stop("phi must be a matrix or data.frame.")
    }
    phi <- as.matrix(phi)

    if (!is.numeric(phi)) {
        stop("phi must be numeric.")
    }

    if (is.null(outcome) || length(outcome) != nrow(phi)) {
        stop("outcome must have length nrow(phi).")
    }
    outcome <- as.integer(outcome)

    if (anyNA(outcome)) stop("outcome contains NA values.")
    if (!all(outcome %in% c(0L, 1L))) {
        stop("outcome must be binary with values 0 (term) and 1 (sPTB).")
    }

    if (is.null(vertices) || length(vertices) == 0L) {
        stop("vertices must be a non-empty integer vector of row indices.")
    }
    vertices <- as.integer(vertices)
    if (anyNA(vertices)) stop("vertices contains NA values.")
    if (any(vertices < 1L | vertices > nrow(phi))) stop("vertices out of range for phi rows.")

    if (!is.null(dist)) {
        if (!is.numeric(dist) || length(dist) != nrow(phi)) {
            stop("dist must be a numeric vector of length nrow(phi) (or NULL).")
        }
        if (anyNA(dist)) stop("dist contains NA values.")
    }

    n.bins <- as.integer(n.bins)
    if (!is.finite(n.bins) || n.bins < 2L) stop("n.bins must be an integer >= 2.")

    min.sptb.detections <- as.integer(min.sptb.detections)
    min.term.detections <- as.integer(min.term.detections)
    min.bins <- as.integer(min.bins)

    if (!is.finite(min.prevalence) || min.prevalence < 0 || min.prevalence > 1) {
        stop("min.prevalence must be in [0,1].")
    }

    ## ------------------------------------------------------------------------
    ## Subset to the vertex set
    ## ------------------------------------------------------------------------

    agg.outcome <- outcome[vertices]
    agg.phi <- phi[vertices, , drop = FALSE]

    n.samples <- nrow(agg.phi)
    n.sptb <- sum(agg.outcome == 1L)
    n.term <- sum(agg.outcome == 0L)

    if (isTRUE(verbose)) {
        cat(sprintf("Samples: %d total (%d sPTB, %d term)\n", n.samples, n.sptb, n.term))
    }

    ## ------------------------------------------------------------------------
    ## Criterion 1: Overall prevalence
    ## ------------------------------------------------------------------------

    phylo.prevalence <- colSums(agg.phi > 0)
    phylo.prevalence.prop <- phylo.prevalence / n.samples

    ## ------------------------------------------------------------------------
    ## Criterion 2: Detection in both outcome groups
    ## ------------------------------------------------------------------------

    if (n.sptb > 0L) {
        phylo.n.sptb <- colSums(agg.phi[agg.outcome == 1L, , drop = FALSE] > 0)
    } else {
        phylo.n.sptb <- rep(0L, ncol(agg.phi))
        names(phylo.n.sptb) <- colnames(agg.phi)
    }

    if (n.term > 0L) {
        phylo.n.term <- colSums(agg.phi[agg.outcome == 0L, , drop = FALSE] > 0)
    } else {
        phylo.n.term <- rep(0L, ncol(agg.phi))
        names(phylo.n.term) <- colnames(agg.phi)
    }

    ## ------------------------------------------------------------------------
    ## Criterion 3: Spatial coverage (optional)
    ## ------------------------------------------------------------------------

    use.spatial <- !is.null(dist)

    if (use.spatial) {

        agg.dist <- dist[vertices]

        ## Quantile breaks may contain duplicates if agg.dist has ties; handle robustly.
        dist.breaks <- as.numeric(stats::quantile(agg.dist,
                                                  probs = seq(0, 1, length.out = n.bins + 1L),
                                                  names = FALSE,
                                                  type = 7))
        dist.breaks <- unique(dist.breaks)

        if (length(dist.breaks) < 2L) {
            stop("dist values are too degenerate to form bins (all equal).")
        }

        ## If too few unique breaks, reduce bins to feasible number
        if (length(dist.breaks) < (n.bins + 1L)) {
            if (isTRUE(verbose)) {
                cat(sprintf("Note: reducing n.bins from %d to %d due to tied quantiles.\n",
                            n.bins, length(dist.breaks) - 1L))
            }
            n.bins <- length(dist.breaks) - 1L
        }

        dist.bins <- cut(agg.dist,
                         breaks = dist.breaks,
                         include.lowest = TRUE,
                         labels = FALSE)

        phylo.n.bins <- apply(agg.phi > 0, 2, function(detected) {
            if (!any(detected)) return(0L)
            length(unique(dist.bins[detected]))
        })

    } else {
        phylo.n.bins <- rep(NA_integer_, ncol(agg.phi))
        names(phylo.n.bins) <- colnames(agg.phi)
    }

    ## ------------------------------------------------------------------------
    ## Combined filtering
    ## ------------------------------------------------------------------------

    filter.prevalence <- phylo.prevalence.prop >= min.prevalence
    filter.sptb <- phylo.n.sptb >= min.sptb.detections
    filter.term <- phylo.n.term >= min.term.detections

    if (use.spatial) {
        filter.spatial <- phylo.n.bins >= min.bins
    } else {
        filter.spatial <- rep(TRUE, ncol(agg.phi))
        names(filter.spatial) <- colnames(agg.phi)
    }

    filter.all <- filter.prevalence & filter.sptb & filter.term & filter.spatial

    if (isTRUE(verbose)) {
        cat("\nPhylotype filtering summary:\n")
        cat(sprintf("  Total phylotypes: %d\n", ncol(agg.phi)))
        cat(sprintf("  Pass prevalence (>= %.0f%%): %d\n", 100 * min.prevalence, sum(filter.prevalence)))
        cat(sprintf("  Pass sPTB detection (>= %d): %d\n", min.sptb.detections, sum(filter.sptb)))
        cat(sprintf("  Pass term detection (>= %d): %d\n", min.term.detections, sum(filter.term)))
        if (use.spatial) {
            cat(sprintf("  Pass spatial coverage (>= %d bins): %d\n", min.bins, sum(filter.spatial)))
        } else {
            cat("  Spatial coverage: skipped (dist is NULL)\n")
        }
        cat(sprintf("  Pass all criteria: %d\n", sum(filter.all)))
    }

    selected.phylotypes <- colnames(agg.phi)[filter.all]
    phi.filtered <- agg.phi[, filter.all, drop = FALSE]

    if (isTRUE(verbose)) {
        cat(sprintf("\nSelected %d phylotypes for association analysis\n", length(selected.phylotypes)))
    }

    ## ------------------------------------------------------------------------
    ## Summary table (selected.phylotypes.df)
    ## ------------------------------------------------------------------------

    base.rate <- if (n.samples > 0L) (n.sptb / n.samples) else NA_real_

    ## Core columns always present
    selected.phylotypes.df <- data.frame(
        n.prevalence = phylo.prevalence[selected.phylotypes],
        p.prevalence = phylo.prevalence.prop[selected.phylotypes],
        n.sptb = phylo.n.sptb[selected.phylotypes],
        n.term = phylo.n.term[selected.phylotypes],
        stringsAsFactors = FALSE
    )

    ## Optional spatial column
    if (use.spatial) {
        selected.phylotypes.df$n.bins <- phylo.n.bins[selected.phylotypes]
    }

    ## Add ratio and rel.ratio
    ## ratio: per-phylotype sPTB/term detection ratio (avoid division by zero)
    ratio <- rep(NA_real_, nrow(selected.phylotypes.df))
    idx.term.pos <- selected.phylotypes.df$n.term > 0
    ratio[idx.term.pos] <- selected.phylotypes.df$n.sptb[idx.term.pos] / selected.phylotypes.df$n.term[idx.term.pos]
    selected.phylotypes.df$ratio <- ratio

    ## rel.ratio: ratio normalized by the baseline sPTB rate among selected vertices
    ## rel.ratio = (n.sptb/n.term) / (n.sptb.total/n.samples.total)
    if (is.finite(base.rate) && base.rate > 0) {
        selected.phylotypes.df$rel.ratio <- selected.phylotypes.df$ratio / base.rate
    } else {
        selected.phylotypes.df$rel.ratio <- NA_real_
    }

    ## If spatial column exists but is entirely NA, drop it
    if ("n.bins" %in% colnames(selected.phylotypes.df)) {
        if (all(is.na(selected.phylotypes.df$n.bins))) {
            selected.phylotypes.df$n.bins <- NULL
        }
    }

    ## Row names
    rownames(selected.phylotypes.df) <- selected.phylotypes

    ## Sort
    if (!sort.by %in% colnames(selected.phylotypes.df)) {
        stop(sprintf("sort.by='%s' is not a column of selected.phylotypes.df.", sort.by))
    }
    o <- order(selected.phylotypes.df[[sort.by]], decreasing = isTRUE(decreasing), na.last = TRUE)
    selected.phylotypes.df <- selected.phylotypes.df[o, , drop = FALSE]

    out <- list(
        selected.phylotypes = selected.phylotypes,
        phi.filtered = phi.filtered,
        selected.phylotypes.df = selected.phylotypes.df,
        filters = list(
            prevalence = filter.prevalence,
            sptb = filter.sptb,
            term = filter.term,
            spatial = filter.spatial,
            all = filter.all
        ),
        metrics = list(
            phylo.prevalence = phylo.prevalence,
            phylo.prevalence.prop = phylo.prevalence.prop,
            phylo.n.sptb = phylo.n.sptb,
            phylo.n.term = phylo.n.term,
            phylo.n.bins = phylo.n.bins
        ),
        thresholds = list(
            min.prevalence = min.prevalence,
            min.sptb.detections = min.sptb.detections,
            min.term.detections = min.term.detections,
            min.bins = if (use.spatial) min.bins else NA_integer_,
            n.bins = if (use.spatial) n.bins else NA_integer_
        ),
        counts = list(
            n.samples = n.samples,
            n.sptb = n.sptb,
            n.term = n.term
        ),
        call = match.call()
    )

    class(out) <- c("phylotype_selection", "list")
    out
}

#' Print Method for phylotype_selection Objects
#'
#' Prints a concise summary of a phylotype selection result produced by
#' \code{select.phylotypes.for.assoc()}.
#'
#' @param x An object of class \code{"phylotype_selection"}.
#' @param ... Additional arguments (ignored).
#' @param top.n Integer; number of top phylotypes (by \code{n.prevalence}) to print.
#'   Default is 10.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.phylotype_selection <- function(x, ..., top.n = 10L) {

    cat("Phylotype selection result\n")
    cat("==========================\n")

    if (!is.null(x$call)) {
        cat("Call:\n  ")
        dput(x$call)
    }

    if (!is.null(x$counts)) {
        cat(sprintf("Samples: %d total (%d sPTB, %d term)\n",
                    x$counts$n.samples, x$counts$n.sptb, x$counts$n.term))
    }

    if (!is.null(x$thresholds)) {
        th <- x$thresholds
        cat("Thresholds:\n")
        cat(sprintf("  min.prevalence: %.3f\n", th$min.prevalence))
        cat(sprintf("  min.sptb.detections: %d\n", th$min.sptb.detections))
        cat(sprintf("  min.term.detections: %d\n", th$min.term.detections))
        if (is.finite(th$n.bins) && is.finite(th$min.bins)) {
            cat(sprintf("  spatial coverage: min.bins=%d of n.bins=%d\n", th$min.bins, th$n.bins))
        } else {
            cat("  spatial coverage: skipped\n")
        }
    }

    n.sel <- length(x$selected.phylotypes %||% character(0))
    cat(sprintf("Selected phylotypes: %d\n", n.sel))

    if (!is.null(x$selected.phylotypes.df) && nrow(x$selected.phylotypes.df) > 0L) {
        top.n <- as.integer(top.n)
        if (!is.finite(top.n) || top.n < 1L) top.n <- 10L
        cat("\nTop selected phylotypes:\n")
        print(utils::head(x$selected.phylotypes.df, top.n))
    }

    invisible(x)
}
