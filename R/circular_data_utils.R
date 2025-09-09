
#' Align Estimated Angles with True Angles
#'
#' @description
#' Aligns a set of estimated angles with reference angles by finding the optimal
#' rotation and handling circular boundaries.
#'
#' @details
#' This function seeks to align angles that represent positions on a circle by:
#' 1. Finding the optimal shift that maximizes correlation with true angles
#' 2. Handling boundary wrapping issues (when angles cross the \eqn{0/2pi} boundary)
#' 3. Optionally scaling the angles to match the range of true angles
#'
#' The function returns either a simple shifted version or a scaled version,
#' depending on which achieves higher correlation with the true angles.
#'
#' @param true.angles A numeric vector of reference angles (in radians)
#' @param est.angles A numeric vector of estimated angles (in radians) to be aligned
#'
#' @return A list containing:
#'   \item{angles}{The aligned angles after shifting and potential scaling}
#'   \item{correlation}{Pearson correlation between true angles and aligned angles}
#'   \item{type}{Character string indicating if the result is "shifted" or "scaled"}
#'
#' @examples
#' # Generate some true angles around a circle
#' true.angles <- seq(0, 2*pi - 0.1, length.out = 20)
#'
#' # Create estimated angles with an arbitrary shift and some noise
#' est.angles <- (true.angles + 1.5) %% (2*pi) + rnorm(20, 0, 0.1)
#'
#' # Align the estimated angles
#' aligned <- angle.alignment(true.angles, est.angles)
#'
#' # Plot results
#' plot(true.angles, ylim=c(-0.5, 2*pi+0.5), pch=19, col="blue",
#'      main="Angle Alignment")
#' points(est.angles, pch=19, col="red")
#' points(aligned$angles, pch=19, col="green")
#' legend("topright", legend=c("True", "Estimated", "Aligned"),
#'        col=c("blue", "red", "green"), pch=19)
#'
#' @export
angle.alignment <- function(true.angles, est.angles) {
    ## First find the approximate shift as before
    shifts <- seq(0, 2*pi, length.out = 500)
    cors <- numeric(length(shifts))

    for(i in seq_along(shifts)) {
        shift <- shifts[i]
        adjusted.angles <- (est.angles + shift) %% (2*pi)
        cors[i] <- cor(true.angles, adjusted.angles, method="pearson")
    }

    best_idx <- which.max(cors)
    best_shift <- shifts[best_idx]

    ## Apply the initial shift
    aligned.angles <- (est.angles + best_shift) %% (2*pi)

    ## Handle boundary wrapping - detect and fix any discontinuity
    ## Find potential wrapping points
    diff.angles <- diff(aligned.angles)
    wrap_points <- which(abs(diff.angles) > pi)

    if(length(wrap_points) > 0) {
        ## Adjust angles after the wrap point
        for(wp in wrap_points) {
            if(diff.angles[wp] > 0) {
                aligned.angles[(wp+1):length(aligned.angles)] <-
                    aligned.angles[(wp+1):length(aligned.angles)] - 2*pi
            } else {
                aligned.angles[(wp+1):length(aligned.angles)] <-
                    aligned.angles[(wp+1):length(aligned.angles)] + 2*pi
            }
        }
    }

    ## Scale to match the range of true angles if needed
    true_range <- max(true.angles) - min(true.angles)
    aligned_range <- max(aligned.angles) - min(aligned.angles)
    scale_factor <- true_range / aligned_range

    ## Center and scale
    aligned.angles_scaled <- (aligned.angles - min(aligned.angles)) * scale_factor + min(true.angles)

    ## Check both versions and return the better one
    cor1 <- cor(true.angles, aligned.angles, method="pearson")
    cor2 <- cor(true.angles, aligned.angles_scaled, method="pearson")

    if(cor2 > cor1) {
        return(list(angles = aligned.angles_scaled, correlation = cor2, type = "scaled"))
    } else {
        return(list(angles = aligned.angles, correlation = cor1, type = "shifted"))
    }
}
