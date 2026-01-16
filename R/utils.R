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
