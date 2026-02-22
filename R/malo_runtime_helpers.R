#' @keywords internal
.gflow.require.malo <- function(api = NULL) {
    if (!requireNamespace("malo", quietly = TRUE)) {
        if (is.null(api) || !nzchar(api)) {
            stop(
                "This functionality requires package 'malo'. Install it first (e.g. remotes::install_github('pgajer/malo')).",
                call. = FALSE
            )
        }
        stop(
            sprintf(
                "'%s' requires package 'malo'. Install it first (e.g. remotes::install_github('pgajer/malo')).",
                api
            ),
            call. = FALSE
        )
    }
}
