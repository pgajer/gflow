#' Compute Basin and Cell Membership from Gradient Basins
#'
#' @description
#' Converts gradient basin structures into soft membership vectors suitable for
#' association analysis. This function handles basin multiplicity by assigning
#' normalized membership weights to each basin containing a vertex.
#'
#' @details
#' In discrete Morse theory on graphs, a vertex may belong to multiple basins
#' simultaneously. This occurs when the vertex has multiple neighbors with lower
#' (for descending basins) or higher (for ascending basins) function values,
#' leading to non-unique gradient trajectories.
#'
#' The membership vectors normalize this multiplicity. If vertex v belongs to
#' basins B_1, ..., B_k, the membership weight for each basin is 1/k (uniform
#' weighting). This soft assignment allows downstream computations to properly
#' account for ambiguity in basin boundaries.
#'
#' Cell membership is computed as the product of basin memberships. A Morse-Smale
#' cell \eqn{C_{ij}} is the intersection of maximum basin \eqn{B_i^+} and
#' minimum basin \eqn{B_j^-}. The cell membership \eqn{\gamma_{ij}(v)} is
#' normalized so that \eqn{\sum_{i,j}\gamma_{ij}(v)=1}.
#'
#' @param basins Object of class \code{"basins_of_attraction"} from
#'   \code{compute.basins.of.attraction}.
#'
#' @return An object of class \code{"gfassoc_membership"} containing:
#'   \item{max_basin_indices}{List of integer vectors. Element v contains 0-based
#'     indices of maximum basins containing vertex v.}
#'   \item{min_basin_indices}{List of integer vectors. Element v contains 0-based
#'     indices of minimum basins containing vertex v.}
#'   \item{max_membership}{List of numeric vectors. Element v contains normalized
#'     membership weights for maximum basins.}
#'   \item{min_membership}{List of numeric vectors. Element v contains normalized
#'     membership weights for minimum basins.}
#'   \item{max_vertices}{Integer vector of extremum vertex indices (1-based).}
#'   \item{min_vertices}{Integer vector of extremum vertex indices (1-based).}
#'   \item{max_values}{Numeric vector of function values at maxima.}
#'   \item{min_values}{Numeric vector of function values at minima.}
#'   \item{n_max_basins}{Number of maximum basins.}
#'   \item{n_min_basins}{Number of minimum basins.}
#'   \item{cell_indices}{List of 2-column matrices. Element v contains (max_idx, min_idx)
#'     pairs for cells containing vertex v (0-based indices).}
#'   \item{cell_membership}{List of numeric vectors. Element v contains normalized
#'     cell membership weights.}
#'
#' @examples
#' \dontrun{
#' ## Compute basins
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#'
#' ## Extract membership structure
#' membership <- gfassoc.membership(basins)
#'
#' ## Check multiplicity for a specific vertex
#' v <- 100
#' n.max.basins <- length(membership$max_basin_indices[[v]])
#' n.cells <- nrow(membership$cell_indices[[v]])
#' }
#'
#' @seealso
#' \code{\link{compute.basins.of.attraction}} for computing gradient basins,
#' \code{\link{gfassoc.polarity}} for computing polarity coordinates
#'
#' @export
gfassoc.membership <- function(basins) {

    if (!inherits(basins, "basins_of_attraction")) {
        stop("basins must be of class 'basins_of_attraction'")
    }

    result <- .Call(
        S_gfassoc_membership,
        basins$lmax_basins,
        basins$lmin_basins,
        as.integer(basins$n_vertices),
        PACKAGE = "gflow"
    )

    return(result)
}

#' Compute Polarity Coordinates for a Fitted Surface
#'
#' @description
#' Computes the polarity coordinate p(v) in \eqn{[-1, 1]} for each vertex, measuring
#' where the vertex sits within its accessible dynamic range defined by the
#' gradient flow structure.
#'
#' @details
#' The polarity coordinate captures the relative position of a vertex within
#' its gradient flow cell. For a vertex v in cell \eqn{C_{ij}} (between maximum M_i
#' and minimum m_j), the normalized height is:
#'
#'   \deqn{ theta(v) = (y(v) - y(m_j)) / (y(M_i) - y(m_j)) }
#'
#' This value lies in \eqn{[0, 1]} by the gradient flow property: any vertex reachable
#' from both M_i (via descending flow) and m_j (via ascending flow) must have
#' function value between y(m_j) and y(M_i).
#'
#' The polarity coordinate rescales to \eqn{[-1, 1]}:
#'
#'   \deqn{p(v) = 2 \cdot theta(v) - 1}
#'
#' So p(v) = +1 at maxima, p(v) = -1 at minima, and p(v) = 0 at the midpoint.
#'
#' When a vertex belongs to multiple cells (basin multiplicity), the polarity
#' is computed as the membership-weighted average of cell-specific polarities.
#'
#' @param y Numeric vector of function values at vertices.
#' @param membership Object of class \code{"gfassoc_membership"} from
#'   \code{gfassoc.membership}.
#' @param polarity.scale Character string: \code{"value"} for normalized height
#'   based on function values, \code{"rank"} for rank-based computation.
#' @param epsilon Numeric threshold for flat region detection.
#'
#' @return An object of class \code{"gfassoc_polarity"} containing:
#'   \item{theta}{Numeric vector of normalized heights in \eqn{[0, 1]}.}
#'   \item{polarity}{Numeric vector of polarity coordinates in \eqn{[-1, 1]}.}
#'   \item{range}{Numeric vector of dynamic ranges M(v) - m(v).}
#'   \item{is_valid}{Logical vector indicating valid polarity computation.}
#'   \item{epsilon}{The epsilon value used.}
#'
#' @examples
#' \dontrun{
#' ## Compute membership and polarity
#' basins <- compute.basins.of.attraction(adj.list, weight.list, y)
#' membership <- gfassoc.membership(basins)
#' polarity <- gfassoc.polarity(y, membership)
#'
#' ## Identify vertices near maxima (high polarity)
#' near.max <- which(polarity$polarity > 0.8 & polarity$is_valid)
#' }
#'
#' @seealso
#' \code{\link{gfassoc.membership}} for computing membership structure,
#' \code{\link{gfcor}} for full correlation analysis
#'
#' @export
gfassoc.polarity <- function(y,
                             membership,
                             polarity.scale = c("value", "rank"),
                             epsilon = 1e-10) {

    polarity.scale <- match.arg(polarity.scale)

    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }

    if (!inherits(membership, "gfassoc_membership")) {
        stop("membership must be of class 'gfassoc_membership'")
    }

    result <- .Call(
        S_gfassoc_polarity,
        as.numeric(y),
        membership,
        polarity.scale,
        as.numeric(epsilon),
        PACKAGE = "gflow"
    )

    return(result)
}

#' Compute Soft Overlap Matrices Between Basin Structures
#'
#' @description
#' Computes soft overlap matrices quantifying how the basin partitions of two
#' functions intersect. The overlap accounts for basin multiplicity through
#' soft membership weights.
#'
#' @details
#' For two functions y and z with basin structures, the overlap matrix entry
#' \eqn{O_{ij}^{++}} measures the effective intersection between y-maximum basin
#' i and z-maximum basin j:
#'
#'   \deqn{O_{ij}^{++} = \sum_v m_0(v)\,\mu_i^{y,+}(v)\,\mu_j^{z,+}(v)}
#'
#' where \eqn{m_0(v)} is vertex mass and \eqn{\mu} are normalized membership
#' weights.
#'
#' Four overlap matrices are computed:
#'   O_pp: y-max with z-max (positive-positive association regions)
#'   O_mm: y-min with z-min (negative-negative association regions)
#'   O_pm: y-max with z-min (opposite association regions)
#'   O_mp: y-min with z-max (opposite association regions)
#'
#' Large O_pp entries indicate regions where both functions tend toward high
#' values. Large O_pm entries indicate regions where y is high but z is low.
#'
#' @param y.membership Membership structure for function y.
#' @param z.membership Membership structure for function z.
#' @param vertex.mass Optional numeric vector of vertex weights.
#'
#' @return A list containing:
#'   \item{O_pp}{Matrix of y-max with z-max overlaps.}
#'   \item{O_mm}{Matrix of y-min with z-min overlaps.}
#'   \item{O_pm}{Matrix of y-max with z-min overlaps.}
#'   \item{O_mp}{Matrix of y-min with z-max overlaps.}
#'   \item{total_mass}{Sum of vertex masses.}
#'
#' @export
gfassoc.overlap <- function(y.membership,
                            z.membership,
                            vertex.mass = NULL) {

    if (!inherits(y.membership, "gfassoc_membership")) {
        stop("y.membership must be of class 'gfassoc_membership'")
    }

    if (!inherits(z.membership, "gfassoc_membership")) {
        stop("z.membership must be of class 'gfassoc_membership'")
    }

    result <- .Call(
        S_gfassoc_overlap,
        y.membership,
        z.membership,
        if (is.null(vertex.mass)) NULL else as.numeric(vertex.mass),
        PACKAGE = "gflow"
    )

    return(result)
}

#' Compute Deviation from Independence for Overlap Matrix
#'
#' @description
#' Analyzes whether basin pairs have more or less overlap than expected under
#' independence. Computes raw deviations and standardized Pearson residuals.
#'
#' @details
#' Under independence, the expected overlap is:
#'
#'   \deqn{E_{ij} = \frac{(\sum_k O_{ik})(\sum_l O_{lj})}{\sum_{k,l} O_{kl}}}
#'
#' The raw deviation \eqn{\delta_{ij} = O_{ij} - E_{ij}} measures the difference
#' from expected. The standardized deviation \eqn{\zeta_{ij}} follows the Pearson residual
#' form, providing a scale-free measure suitable for comparing across different
#' matrix sizes.
#'
#' Large positive zeta values indicate basin pairs that overlap more than
#' expected (attraction). Large negative values indicate pairs that overlap
#' less than expected (repulsion).
#'
#' @param overlap.matrix A numeric matrix (one of the overlap matrices from
#'   \code{gfassoc.overlap}).
#'
#' @return A list containing:
#'   \item{delta}{Matrix of raw deviations \eqn{O_{ij} - E_{ij}}.}
#'   \item{zeta}{Matrix of standardized Pearson residuals.}
#'   \item{expected}{Matrix of expected overlaps under independence.}
#'
#' @export
gfassoc.deviation <- function(overlap.matrix) {

    if (!is.matrix(overlap.matrix) || !is.numeric(overlap.matrix)) {
        stop("overlap.matrix must be a numeric matrix")
    }

    result <- .Call(
        S_gfassoc_deviation,
        overlap.matrix,
        PACKAGE = "gflow"
    )

    return(result)
}

#' Compute Basin Association Character
#'
#' @description
#' Computes the association character for each basin, measuring the average
#' polarity of the other function within that basin.
#'
#' @details
#' The association character of y-maximum basin i is:
#'
#'   \deqn{
#'   \chi_i^{y,+} =
#'   \frac{\sum_v m_0(v)\,\mu_i^{y,+}(v)\,p_z(v)}
#'        {\sum_v m_0(v)\,\mu_i^{y,+}(v)}
#'   }
#'
#' This measures the average z-polarity within the y-maximum basin. Values near
#' +1 indicate that the y-high region coincides with z-high (direct association).
#' Values near -1 indicate that y-high coincides with z-low (opposite association).
#'
#' Characters are computed for all four basin types: y-max, y-min, z-max, z-min.
#' The characters provide a basin-level summary of the association structure,
#' identifying which basins drive positive versus negative association.
#'
#' @param y.membership Membership structure for function y.
#' @param z.membership Membership structure for function z.
#' @param pol.y Polarity structure for function y.
#' @param pol.z Polarity structure for function z.
#' @param vertex.mass Optional numeric vector of vertex weights.
#'
#' @return A list containing:
#'   \item{chi_y_max}{Character of each y-maximum basin.}
#'   \item{chi_y_min}{Character of each y-minimum basin.}
#'   \item{chi_z_max}{Character of each z-maximum basin.}
#'   \item{chi_z_min}{Character of each z-minimum basin.}
#'   \item{mass_y_max}{Mass of each y-maximum basin.}
#'   \item{mass_y_min}{Mass of each y-minimum basin.}
#'   \item{mass_z_max}{Mass of each z-maximum basin.}
#'   \item{mass_z_min}{Mass of each z-minimum basin.}
#'
gfassoc.basin.character <- function(y.membership,
                                    z.membership,
                                    pol.y,
                                    pol.z,
                                    vertex.mass = NULL) {

    result <- .Call(
        S_gfassoc_basin_character,
        y.membership,
        z.membership,
        pol.y,
        pol.z,
        if (is.null(vertex.mass)) NULL else as.numeric(vertex.mass),
        PACKAGE = "gflow"
    )

    return(result)
}
