.validate.s3.contract <- function(x,
                                  class.name,
                                  fields = character(),
                                  attributes = character(),
                                  storage = "any") {
    if (!inherits(x, class.name)) {
        stop("Expected an object inheriting from '", class.name, "'.")
    }

    storage.ok <- switch(
        storage,
        any = TRUE,
        list = is.list(x),
        matrix = is.matrix(x),
        array = is.array(x),
        matrix_or_list = is.matrix(x) || is.list(x),
        array_or_list = is.array(x) || is.list(x),
        stop("Unknown S3 storage contract: ", storage)
    )
    if (!storage.ok) {
        stop("Invalid storage for class '", class.name, "': expected ", storage, ".")
    }

    missing.fields <- setdiff(fields, names(x))
    if (length(missing.fields)) {
        stop(
            "Invalid '", class.name, "' object; missing fields: ",
            paste(missing.fields, collapse = ", "), "."
        )
    }

    missing.attributes <- attributes[
        vapply(attributes, function(name) is.null(attr(x, name)), logical(1L))
    ]
    if (length(missing.attributes)) {
        stop(
            "Invalid '", class.name, "' object; missing attributes: ",
            paste(missing.attributes, collapse = ", "), "."
        )
    }

    invisible(x)
}

.validate.gfc.modulation <- function(x) {
    .validate.s3.contract(
        x, "gfc_modulation",
        fields = c("modulated.weight.list", "modulation", "n.vertices", "n.edges"),
        storage = "list"
    )
}

.validate.distance.quantile.bins <- function(x) {
    .validate.s3.contract(
        x, "distance.quantile.bins",
        fields = c("bins", "tests", "bin.assignment", "keep"),
        storage = "list"
    )
    if (!is.data.frame(x$bins)) {
        stop("Invalid 'distance.quantile.bins' object; 'bins' must be a data frame.")
    }
    invisible(x)
}

.validate.gfcor <- function(x) {
    .validate.s3.contract(
        x, "gfcor",
        fields = c("global", "basin_character", "n_vertices", "polarity_scale", "epsilon"),
        storage = "list"
    )
}

.validate.summary.gfcor <- function(x) {
    .validate.s3.contract(
        x, "summary.gfcor",
        fields = c("global", "vertex_distribution"),
        storage = "list"
    )
}

.validate.gflow.graph <- function(x) {
    .validate.s3.contract(
        x, "gflow_graph",
        fields = c(
            "adjacency.list", "weight.list", "intersection.matrix",
            "basin.metadata", "n.ascending", "n.descending",
            "edge.type", "min.intersection"
        ),
        storage = "list"
    )
}

.validate.summary.gflow.graph <- function(x) {
    .validate.s3.contract(
        x, "summary.gflow_graph",
        fields = c("n.vertices", "n.edges", "degree.stats", "basin.size.stats"),
        storage = "list"
    )
}

.validate.harmonic.extension <- function(x) {
    .validate.s3.contract(
        x, "harmonic_extension",
        fields = c(
            "trajectory", "trajectory.length", "tubular.vertices",
            "extended.coords", "hop.distances", "geodesic.distances",
            "nearest.traj.idx"
        ),
        storage = "list"
    )
}

.validate.lcor.posterior <- function(x) {
    .validate.s3.contract(
        x, "lcor.posterior",
        fields = c(
            "mean", "sd", "lower", "upper", "n.samples",
            "n.features", "n.vertices", "credible.level", "mode"
        ),
        storage = "list"
    )
}

.validate.summary.lcor.posterior <- function(x) {
    .validate.s3.contract(
        x, "summary.lcor.posterior",
        fields = c(
            "feature.summary", "prop.significant.overall",
            "n.features", "n.vertices", "n.samples", "credible.level"
        ),
        storage = "list"
    )
}

.validate.lcor.vector.matrix.result <- function(x) {
    .validate.s3.contract(
        x, "lcor_vector_matrix_result",
        attributes = c("n.vertices", "n.columns", "instrumented"),
        storage = "matrix_or_list"
    )
    if (isTRUE(attr(x, "instrumented")) &&
        (!is.list(x) || is.null(x$column.coefficients))) {
        stop("Instrumented lcor vector-matrix results require 'column.coefficients'.")
    }
    invisible(x)
}

.validate.lcor.matrix.matrix.result <- function(x) {
    .validate.s3.contract(
        x, "lcor_matrix_matrix_result",
        attributes = c("n.vertices", "symmetric", "instrumented"),
        storage = "array_or_list"
    )
    if (isTRUE(attr(x, "instrumented")) &&
        (!is.list(x) || is.null(x$pair.coefficients))) {
        stop("Instrumented lcor matrix-matrix results require 'pair.coefficients'.")
    }
    invisible(x)
}

.validate.lslope.gradient.result <- function(x) {
    .validate.s3.contract(
        x, "lslope_gradient_result",
        fields = c(
            "vertex.coefficients", "n.vertices", "n.local.maxima",
            "n.local.minima", "mean.coefficient", "median.coefficient"
        ),
        storage = "list"
    )
}

.validate.lslope.vector.matrix.result <- function(x) {
    .validate.s3.contract(
        x, "lslope_vector_matrix_result",
        attributes = c("n.vertices", "n.columns"),
        storage = "matrix"
    )
}

.validate.lslope.matrix.matrix.result <- function(x) {
    .validate.s3.contract(
        x, "lslope_matrix_matrix_result",
        attributes = c("n.vertices", "n.cols.y", "n.cols.z"),
        storage = "array"
    )
}

.validate.lslope.neighborhood.result <- function(x) {
    .validate.s3.contract(
        x, "lslope_neighborhood_result",
        fields = c("vertex.coefficients", "mean.coefficient", "median.coefficient", "lcor"),
        attributes = c("n.vertices", "weight.type"),
        storage = "list"
    )
}

.validate.madag <- function(x) {
    .validate.s3.contract(
        x, "madag",
        fields = c("source.vertex", "n.vertices", "cells"),
        storage = "list"
    )
}

.validate.weighted.p.summary <- function(x) {
    .validate.s3.contract(
        x, "weighted.p.summary",
        fields = c(
            "weighted.p.value", "classical.p.value", "summary.stats",
            "interpretation", "alternative", "null.parameters"
        ),
        storage = "list"
    )
}

.validate.se.tree <- function(x) {
    .validate.s3.contract(
        x, "se_tree",
        fields = c(
            "root.vertex", "root.is.maximum", "classification", "nodes",
            "ns.min.terminals", "ns.max.terminals"
        ),
        storage = "list"
    )
}
