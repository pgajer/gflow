#' Find Basin Index by Label
#'
#' @description
#' Finds the index of a basin in a basin list by its label.
#'
#' @param label Character string, basin label (e.g., "m1", "M3")
#' @param basins.obj Object of class "basins_of_attraction"
#' @param extrema.type Character string, either "min" or "max"
#'
#' @return Integer index of the basin, or NULL if not found
#'
find.basin.idx.by.label <- function(label, basins.obj, extrema.type) {
    if (!extrema.type %in% c("max", "min")) {
        stop("extrema.type must be either 'max' or 'min'")
    }
    
    if (extrema.type == "max") {
        basin.list <- basins.obj$lmax_basins
    } else {
        basin.list <- basins.obj$lmin_basins
    }
    
    for (i in seq_along(basin.list)) {
        ## Check if basin has label attribute
        if (!is.null(names(basin.list)[i]) && names(basin.list)[i] == label) {
            return(i)
        }
    }
    
    return(NULL)
}

#' Find Basin Index by Vertex
#'
#' @description
#' Finds the index of a basin in a basin list by its extremum vertex.
#'
#' @param vertex Integer vertex identifier (1-based)
#' @param basins.obj Object of class "basins_of_attraction"
#' @param extrema.type Character string, either "min" or "max"
#'
#' @return Integer index of the basin, or NULL if not found
#'
find.basin.idx.by.vertex <- function(vertex, basins.obj, extrema.type) {
    if (!extrema.type %in% c("max", "min")) {
        stop("extrema.type must be either 'max' or 'min'")
    }
    
    if (extrema.type == "max") {
        basin.list <- basins.obj$lmax_basins
    } else {
        basin.list <- basins.obj$lmin_basins
    }
    
    for (i in seq_along(basin.list)) {
        if (basin.list[[i]]$vertex == vertex) {
            return(i)
        }
    }
    
    return(NULL)
}

#' Find Basin by Label
#'
#' @description
#' Retrieves a basin structure by its label.
#'
#' @param label Character string, basin label (e.g., "m1", "M3")
#' @param basins.obj Object of class "basins_of_attraction"
#' @param extrema.type Character string, either "min" or "max"
#'
#' @return Basin structure (list), or NULL if not found
#'
find.basin.by.label <- function(label, basins.obj, extrema.type) {
    idx <- find.basin.idx.by.label(label, basins.obj, extrema.type)
    
    if (is.null(idx)) {
        return(NULL)
    }
    
    if (extrema.type == "max") {
        return(basins.obj$lmax_basins[[idx]])
    } else {
        return(basins.obj$lmin_basins[[idx]])
    }
}

#' Find Basin by Vertex
#'
#' @description
#' Retrieves a basin structure by its extremum vertex.
#'
#' @param vertex Integer vertex identifier (1-based)
#' @param basins.obj Object of class "basins_of_attraction"
#' @param extrema.type Character string, either "min" or "max"
#'
#' @return Basin structure (list), or NULL if not found
#'
find.basin.by.vertex.in.basins.of.attraction <- function(vertex, basins.obj, extrema.type) {
    idx <- find.basin.idx.by.vertex(vertex, basins.obj, extrema.type)
    
    if (is.null(idx)) {
        return(NULL)
    }
    
    if (extrema.type == "max") {
        return(basins.obj$lmax_basins[[idx]])
    } else {
        return(basins.obj$lmin_basins[[idx]])
    }
}

#' Get Basin Label from Vertex
#'
#' @description
#' Retrieves the label of a basin given its extremum vertex.
#'
#' @param vertex Integer vertex identifier (1-based)
#' @param basin.summary Data frame from summary.basins_of_attraction()
#'
#' @return Character string label, or NULL if not found
#'
get.basin.label.from.vertex <- function(vertex, basin.summary) {
    row <- basin.summary[basin.summary$vertex == vertex, ]
    
    if (nrow(row) == 0) {
        return(NULL)
    }
    
    return(row$label[1])
}

#' Get Basin Vertex from Label
#'
#' @description
#' Retrieves the vertex of a basin given its label.
#'
#' @param label Character string, basin label
#' @param basin.summary Data frame from summary.basins_of_attraction()
#'
#' @return Integer vertex identifier, or NULL if not found
#'
get.basin.vertex.from.label <- function(label, basin.summary) {
    row <- basin.summary[basin.summary$label == label, ]
    
    if (nrow(row) == 0) {
        return(NULL)
    }
    
    return(row$vertex[1])
}

#' Determine Extrema Type from Label
#'
#' @description
#' Determines whether a basin label corresponds to a minimum or maximum.
#'
#' @param label Character string, basin label (e.g., "m1" for min, "M3" for max)
#'
#' @return Character string "min" or "max", or NULL if cannot be determined
#'
get.extrema.type.from.label <- function(label) {
    if (grepl("^m[0-9]+$", label)) {
        return("min")
    } else if (grepl("^M[0-9]+$", label)) {
        return("max")
    } else {
        return(NULL)
    }
}
