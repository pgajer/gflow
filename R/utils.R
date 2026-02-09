#' Format and Print Elapsed Time
#'
#' This function calculates the elapsed time since a given start time,
#' formats it, and prints it with an optional message.
#'
#' @param start.time A proc.time() object representing the start time.
#' @param message Character string to be printed before the elapsed time.
#'   Default is "DONE".
#' @param with.brackets Logical, if TRUE (default), the time is enclosed
#'   in parentheses.
#'
#' @return None. This function is called for its side effect of printing.
#'
#' @examples
#' \dontrun{
#' start <- proc.time()
#' Sys.sleep(2)  # Do some work
#' elapsed.time(start, "Processing complete")
#' }
#' @export
elapsed.time <- function(start.time,
                         message = "DONE",
                         with.brackets = TRUE) {
    elapsed <- as.numeric(proc.time() - start.time)[3]
    minutes <- floor(elapsed / 60)
    seconds <- floor(elapsed %% 60)

    time.str <- sprintf("%d:%02d", minutes, seconds)

    if (with.brackets) {
        output <- sprintf("%s (%s)", message, time.str)
    } else {
        output <- sprintf("%s %s", message, time.str)
    }

    cat(output, "\n")
}

# ============================================================================
# UTILITY OPERATOR
# ============================================================================

# Default value operator (internal use)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Round Numeric Columns to a Fixed Number of Significant Digits
#'
#' @description
#' Returns a copy of a data.frame in which numeric columns are rounded to a specified
#' number of significant digits using \code{\link[base]{signif}}. Columns that are
#' integer-like (all finite values satisfy \code{x == floor(x)}) can optionally be
#' retained as integers.
#'
#' @param df A data.frame.
#' @param digits Integer. Number of significant digits to keep. Default is 2.
#' @param keep.integers Logical. If TRUE (default), numeric columns that are integer-like
#'   are coerced to integer and not rounded.
#' @param exclude Character vector of column names to leave unchanged (useful for IDs).
#'   Default is NULL.
#' @param include Character vector of column names to process. If provided, only these
#'   columns are processed. Default is NULL (process all numeric columns not excluded).
#'
#' @return A data.frame of the same dimensions as \code{df}, with selected numeric columns
#'   transformed.
#'
#' @examples
#' \dontrun{
#' df.2sig <- signif.df(df, digits = 2)
#' df.3sig <- signif.df(df, digits = 3, exclude = c("sector"))
#' }
#'
#' @export
signif.df <- function(df,
                      digits = 2L,
                      keep.integers = TRUE,
                      exclude = NULL,
                      include = NULL) {

    if (!is.data.frame(df)) {
        stop("df must be a data.frame.")
    }

    digits <- as.integer(digits)
    if (length(digits) != 1L || is.na(digits) || digits < 1L) {
        stop("digits must be a single integer >= 1.")
    }

    out <- df

    num.cols <- vapply(out, is.numeric, logical(1))

    if (!is.null(include)) {
        if (!is.character(include)) stop("include must be a character vector of column names.")
        num.cols <- num.cols & names(out) %in% include
    }

    if (!is.null(exclude)) {
        if (!is.character(exclude)) stop("exclude must be a character vector of column names.")
        num.cols <- num.cols & !(names(out) %in% exclude)
    }

    if (!any(num.cols)) {
        return(out)
    }

    out[num.cols] <- lapply(out[num.cols], function(x) {

        if (!isTRUE(keep.integers)) {
            return(signif(x, digits = digits))
        }

        ## Treat as integer-like if all non-NA values are whole numbers
        x.ok <- x[!is.na(x)]
        if (length(x.ok) > 0L && all(x.ok == floor(x.ok))) {
            return(as.integer(x))
        }

        signif(x, digits = digits)
    })

    out
}


#' Map numeric values to a fixed color palette (min -> first, max -> last)
#'
#' Linearly rescales `x` into \eqn{[0, 1]} using `limits`, then assigns each value
#' to a palette entry. By default, the smallest value maps to `palette[1]`
#' and the largest maps to `palette[length(palette)]`.
#'
#' @param x Numeric vector to map.
#' @param color.palette Character vector of colors (e.g., "#RRGGBB"), length >= 2.
#' @param limits Numeric length-2 vector giving c(min, max) for mapping. If NULL,
#'   uses range(x, na.rm=TRUE).
#' @param na.color Color to use for NA / non-finite values in `x`.
#' @param clip Logical; if TRUE, values outside `limits` are clamped to the
#'   nearest endpoint color.
#'
#' @return Character vector of colors, same length as `x`.
#' @examples
#' cols <- map.values.to.palette(rel.sptb.hat, blue.yellow.red.color.palette)
#' head(cols)
map.values.to.palette <- function(x,
                                  color.palette,
                                  limits = NULL,
                                  na.color = NA_character_,
                                  clip = TRUE) {
  ## Input validation
  if (!is.numeric(x)) stop("`x` must be numeric.")
  if (!is.character(color.palette) || length(color.palette) < 2L) {
    stop("`color.palette` must be a character vector of length >= 2.")
  }

  if (is.null(limits)) {
    limits <- range(x, na.rm = TRUE)
  }
  if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits)) {
    stop("`limits` must be a numeric vector of length 2 with no NA.")
  }

  lo <- limits[1]
  hi <- limits[2]
  if (!is.finite(lo) || !is.finite(hi)) stop("`limits` must be finite.")
  if (hi < lo) stop("`limits[2]` must be >= `limits[1]`.")

  n.colors <- length(color.palette)

  ## Degenerate case: all values equal
  if (hi == lo) {
    out <- rep(na.color, length(x))
    ok <- is.finite(x)
    out[ok] <- color.palette[1L]
    return(out)
  }

  ## Rescale to [0, 1]
  t <- (x - lo) / (hi - lo)

  ## Clamp if requested
  if (isTRUE(clip)) {
    t <- pmin(pmax(t, 0), 1)
  }

  ## Map to palette indices (min -> 1, max -> n.colors)
  idx <- floor(t * (n.colors - 1L)) + 1L
  idx <- pmin(pmax(idx, 1L), n.colors)

  ## Fill output
  out <- rep(na.color, length(x))
  ok <- is.finite(x)
  out[ok] <- color.palette[idx[ok]]
  out
}
