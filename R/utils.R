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
#' start <- proc.time()
#' Sys.sleep(2)  # Do some work
#' elapsed.time(start, "Processing complete")
#'
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
