#' Extract Sample and Feature Membership from ISA Biclusters
#'
#' @description
#' Converts ISA biclustering results into interpretable sample clusters (cells)
#' and feature clusters (modules) by thresholding membership scores.
#'
#' @param isa.result Result object from isa::isa()
#' @param score.threshold Minimum absolute score for membership. Default is 0.
#' @param min.samples Minimum number of samples for a valid bicluster. Default is 5.
#' @param min.features Minimum number of features for a valid bicluster. Default is 2.
#'
#' @return List with components:
#' \describe{
#'   \item{biclusters}{List where each element contains sample.ids, feature.ids,
#'     sample.scores, feature.scores for one bicluster}
#'   \item{n.biclusters}{Number of valid biclusters}
#'   \item{bicluster.sizes}{Data frame with rows (samples) and cols (features) per bicluster}
#' }
#'
#' @export
extract.isa.biclusters <- function(isa.result,
                                   score.threshold = 0,
                                   min.samples = 5,
                                   min.features = 2) {

    n.biclusters <- ncol(isa.result$rows)
    biclusters <- vector("list", n.biclusters)

    valid.idx <- integer()

    for (i in seq_len(n.biclusters)) {
        ## Extract membership scores
        row.scores <- isa.result$rows[, i]
        col.scores <- isa.result$columns[, i]

        ## Identify members (above threshold)
        sample.members <- which(abs(row.scores) > score.threshold)
        feature.members <- which(abs(col.scores) > score.threshold)

        ## Check minimum size requirements
        if (length(sample.members) >= min.samples &&
            length(feature.members) >= min.features) {

            biclusters[[i]] <- list(
                sample.ids = sample.members,
                feature.ids = feature.members,
                sample.scores = row.scores[sample.members],
                feature.scores = col.scores[feature.members],
                n.samples = length(sample.members),
                n.features = length(feature.members)
            )

            valid.idx <- c(valid.idx, i)
        }
    }

    ## Keep only valid biclusters
    biclusters <- biclusters[valid.idx]

    ## Create size summary
    bicluster.sizes <- data.frame(
        bicluster = valid.idx,
        n.samples = sapply(biclusters, function(x) x$n.samples),
        n.features = sapply(biclusters, function(x) x$n.features)
    )

    return(list(
        biclusters = biclusters,
        n.biclusters = length(biclusters),
        bicluster.sizes = bicluster.sizes
    ))
}


#' Filter Redundant ISA Biclusters
#'
#' @description
#' ISA often returns many overlapping biclusters. This function identifies
#' and removes highly similar biclusters, keeping the most robust ones.
#'
#' @param isa.result Result object from isa::isa()
#' @param biclusters.extracted Output from extract.isa.biclusters()
#' @param overlap.threshold Maximum Jaccard similarity for considering biclusters
#'   as distinct. Default is 0.8 (biclusters with >80% overlap are merged).
#' @param prefer.robust If TRUE, prefer biclusters with higher robustness scores.
#'   Default is TRUE.
#'
#' @return List similar to extract.isa.biclusters() but with redundant biclusters removed
#'
#' @export
filter.redundant.biclusters <- function(isa.result,
                                       biclusters.extracted,
                                       overlap.threshold = 0.8,
                                       prefer.robust = TRUE) {

    biclusters <- biclusters.extracted$biclusters
    n <- length(biclusters)

    if (n <= 1) {
        return(biclusters.extracted)
    }

    ## Compute pairwise Jaccard similarity for samples
    jaccard.matrix <- matrix(0, n, n)

    for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
            samples.i <- biclusters[[i]]$sample.ids
            samples.j <- biclusters[[j]]$sample.ids

            intersection <- length(intersect(samples.i, samples.j))
            union <- length(union(samples.i, samples.j))

            jaccard <- intersection / union
            jaccard.matrix[i, j] <- jaccard
            jaccard.matrix[j, i] <- jaccard
        }
    }

    ## Identify redundant biclusters
    keep <- rep(TRUE, n)

    for (i in seq_len(n)) {
        if (!keep[i]) next

        ## Find highly overlapping biclusters
        redundant <- which(jaccard.matrix[i, ] > overlap.threshold & seq_len(n) != i)

        if (length(redundant) > 0) {
            ## Decide which to keep based on robustness or size
            if (prefer.robust && !is.null(isa.result$seeddata$rob)) {
                ## Get original indices
                orig.idx <- biclusters.extracted$bicluster.sizes$bicluster

                scores <- isa.result$seeddata$rob[orig.idx]
                keep.this <- which.max(scores[c(i, redundant)])

                if (keep.this == 1) {
                    ## Keep i, remove others
                    keep[redundant] <- FALSE
                } else {
                    ## Keep one of the redundant, remove i and others
                    keep.idx <- redundant[keep.this - 1]
                    keep[c(i, setdiff(redundant, keep.idx))] <- FALSE
                }
            } else {
                ## Keep largest bicluster
                sizes <- sapply(biclusters[c(i, redundant)], function(x) x$n.samples)
                keep.this <- which.max(sizes)

                if (keep.this == 1) {
                    keep[redundant] <- FALSE
                } else {
                    keep.idx <- redundant[keep.this - 1]
                    keep[c(i, setdiff(redundant, keep.idx))] <- FALSE
                }
            }
        }
    }

    ## Filter biclusters
    filtered.biclusters <- biclusters[keep]
    filtered.sizes <- biclusters.extracted$bicluster.sizes[keep, ]
    filtered.sizes$bicluster <- seq_len(sum(keep))

    return(list(
        biclusters = filtered.biclusters,
        n.biclusters = sum(keep),
        bicluster.sizes = filtered.sizes
    ))
}


#' Create Hard Sample Assignment from ISA Biclusters
#'
#' @description
#' Assigns each sample to its best-matching bicluster based on membership scores.
#' This converts overlapping biclusters into non-overlapping clusters similar
#' to Louvain output.
#'
#' @param biclusters.extracted Output from extract.isa.biclusters() or
#'   filter.redundant.biclusters()
#' @param n.samples Total number of samples in original data
#'
#' @return Integer vector of length n.samples with cluster assignments (1-based).
#'   Samples not assigned to any bicluster get value 0.
#'
#' @export
isa.to.hard.clusters <- function(biclusters.extracted, n.samples) {

    biclusters <- biclusters.extracted$biclusters
    n.biclusters <- length(biclusters)

    ## Initialize assignment vector
    assignment <- integer(n.samples)
    max.scores <- numeric(n.samples)

    ## For each sample, find bicluster with highest absolute score
    for (i in seq_along(biclusters)) {
        bc <- biclusters[[i]]

        for (j in seq_along(bc$sample.ids)) {
            sample.id <- bc$sample.ids[j]
            score <- abs(bc$sample.scores[j])

            if (score > max.scores[sample.id]) {
                max.scores[sample.id] <- score
                assignment[sample.id] <- i
            }
        }
    }

    return(assignment)
}


#' Summarize ISA Bicluster Quality
#'
#' @description
#' Computes quality metrics for each bicluster including coherence,
#' enrichment relative to background, and association with outcome.
#'
#' @param isa.result Result object from isa::isa()
#' @param biclusters.extracted Output from extract.isa.biclusters()
#' @param data.matrix Original data matrix (samples × features)
#' @param outcome Optional outcome vector for computing associations
#'
#' @return Data frame with quality metrics per bicluster
#'
#' @export
summarize.bicluster.quality <- function(isa.result,
                                        biclusters.extracted,
                                        data.matrix,
                                        outcome = NULL) {

    biclusters <- biclusters.extracted$biclusters
    n.biclusters <- length(biclusters)

    quality <- data.frame(
        bicluster = seq_len(n.biclusters),
        n.samples = integer(n.biclusters),
        n.features = integer(n.biclusters),
        mean.row.score = numeric(n.biclusters),
        mean.col.score = numeric(n.biclusters),
        coherence = numeric(n.biclusters),
        robustness = numeric(n.biclusters)
    )

    ## Get robustness from seeddata if available
    orig.indices <- biclusters.extracted$bicluster.sizes$bicluster
    if (!is.null(isa.result$seeddata$rob)) {
        quality$robustness <- isa.result$seeddata$rob[orig.indices]
    }

    for (i in seq_len(n.biclusters)) {
        bc <- biclusters[[i]]

        quality$n.samples[i] <- bc$n.samples
        quality$n.features[i] <- bc$n.features
        quality$mean.row.score[i] <- mean(abs(bc$sample.scores))
        quality$mean.col.score[i] <- mean(abs(bc$feature.scores))

        ## Compute coherence (mean correlation within submatrix)
        submatrix <- data.matrix[bc$sample.ids, bc$feature.ids, drop = FALSE]

        if (nrow(submatrix) > 1 && ncol(submatrix) > 1) {
            ## Feature-feature correlation
            if (ncol(submatrix) > 2) {
                feature.cors <- cor(submatrix)
                quality$coherence[i] <- mean(feature.cors[upper.tri(feature.cors)])
            } else {
                quality$coherence[i] <- NA
            }
        } else {
            quality$coherence[i] <- NA
        }
    }

    ## Add outcome association if provided
    if (!is.null(outcome)) {
        quality$outcome.mean <- sapply(biclusters, function(bc) {
            mean(outcome[bc$sample.ids])
        })

        quality$outcome.sd <- sapply(biclusters, function(bc) {
            sd(outcome[bc$sample.ids])
        })
    }

    return(quality)
}

#' Select ISA Biclusters by Robustness
#'
#' @description
#' Uses robustness scores from ISA to select statistically significant biclusters.
#' This is the most principled approach as robustness measures stability under
#' random perturbations.
#'
#' @param isa.result ISA result with seeddata$rob component
#' @param rob.threshold Robustness threshold (higher = more stringent)
#' @param plot.distribution If TRUE, plots robustness distribution
#'
#' @export
select.biclusters.by.robustness <- function(isa.result,
                                           rob.threshold = NULL,
                                           plot.distribution = TRUE) {

    rob.scores <- isa.result$seeddata$rob

    if (plot.distribution) {
        par(mfrow = c(1, 2))

        ## Histogram
        hist(rob.scores, breaks = 50,
             main = "Robustness Score Distribution",
             xlab = "Robustness",
             col = "lightblue")

        if (!is.null(rob.threshold)) {
            abline(v = rob.threshold, col = "red", lwd = 2, lty = 2)
        }

        ## Sorted robustness (elbow plot)
        sorted.rob <- sort(rob.scores, decreasing = TRUE)
        plot(sorted.rob, type = "b", pch = 19,
             xlab = "Bicluster Rank",
             ylab = "Robustness Score",
             main = "Robustness by Rank (Look for Elbow)")

        if (!is.null(rob.threshold)) {
            abline(h = rob.threshold, col = "red", lwd = 2, lty = 2)
        }
    }

    ## Auto-select threshold if not provided
    if (is.null(rob.threshold)) {
        ## Use elbow detection or top quantile
        sorted.rob <- sort(rob.scores, decreasing = TRUE)

        ## Method 1: Find elbow using second derivative
        if (length(sorted.rob) > 10) {
            diffs <- diff(sorted.rob)
            second.diffs <- diff(diffs)
            elbow.idx <- which.max(abs(second.diffs)) + 1
            rob.threshold <- sorted.rob[elbow.idx]

            cat(sprintf("Auto-selected threshold: %.3f (elbow at rank %d)\n",
                       rob.threshold, elbow.idx))
        } else {
            ## Fallback: use 75th percentile
            rob.threshold <- quantile(rob.scores, 0.75)
            cat(sprintf("Auto-selected threshold: %.3f (75th percentile)\n",
                       rob.threshold))
        }
    }

    selected.idx <- which(rob.scores >= rob.threshold)

    cat(sprintf("\nSelected %d out of %d biclusters (robustness >= %.3f)\n",
               length(selected.idx), length(rob.scores), rob.threshold))

    return(list(
        selected.idx = selected.idx,
        threshold = rob.threshold,
        robustness.scores = rob.scores[selected.idx]
    ))
}

#' Select Bicluster Size Thresholds Using Gap Statistics
#'
#' @description
#' Identifies natural breaks in bicluster size distributions to separate
#' signal from noise.
#'
#' @export
select.thresholds.by.gaps <- function(n.samples.per.bc,
                                      n.features.per.bc,
                                      plot.diagnostics = TRUE) {

    ## Find gaps in sample distribution
    sorted.samples <- sort(n.samples.per.bc, decreasing = TRUE)
    sample.gaps <- diff(sorted.samples)

    ## Find largest gap in top half (exclude tiny biclusters)
    top.half.idx <- seq_len(ceiling(length(sorted.samples) / 2))
    largest.gap.idx <- which.max(abs(sample.gaps[top.half.idx]))
    sample.threshold <- sorted.samples[largest.gap.idx + 1]

    ## Find gaps in feature distribution
    sorted.features <- sort(n.features.per.bc, decreasing = TRUE)
    feature.gaps <- diff(sorted.features)

    top.half.idx <- seq_len(ceiling(length(sorted.features) / 2))
    largest.gap.idx <- which.max(abs(feature.gaps[top.half.idx]))
    feature.threshold <- sorted.features[largest.gap.idx + 1]

    if (plot.diagnostics) {
        par(mfrow = c(2, 2))

        ## Sample distribution
        plot(sorted.samples, type = "b", pch = 19,
             main = "Samples per Bicluster (Sorted)",
             xlab = "Rank", ylab = "Number of Samples")
        abline(h = sample.threshold, col = "red", lwd = 2, lty = 2)

        ## Sample gaps
        plot(abs(sample.gaps), type = "h",
             main = "Gaps in Sample Distribution",
             xlab = "Position", ylab = "Gap Size")

        ## Feature distribution
        plot(sorted.features, type = "b", pch = 19,
             main = "Features per Bicluster (Sorted)",
             xlab = "Rank", ylab = "Number of Features")
        abline(h = feature.threshold, col = "red", lwd = 2, lty = 2)

        ## Feature gaps
        plot(abs(feature.gaps), type = "h",
             main = "Gaps in Feature Distribution",
             xlab = "Position", ylab = "Gap Size")
    }

    cat(sprintf("Gap-based sample threshold: %d\n", sample.threshold))
    cat(sprintf("Gap-based feature threshold: %d\n", feature.threshold))

    meets.criteria <- n.samples.per.bc >= sample.threshold &
                      n.features.per.bc >= feature.threshold

    cat(sprintf("\n%d biclusters meet gap-based criteria\n", sum(meets.criteria)))

    return(list(
        sample.threshold = sample.threshold,
        feature.threshold = feature.threshold,
        selected.idx = which(meets.criteria)
    ))
}
