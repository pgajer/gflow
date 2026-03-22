#' @keywords internal
.gflow.require.namespace <- function(pkg, api = NULL, install_hint = NULL) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        return(invisible(TRUE))
    }

    if (is.null(install_hint) || !nzchar(install_hint)) {
        install_hint <- sprintf("Install package '%s' first.", pkg)
    }

    if (is.null(api) || !nzchar(api)) {
        stop(install_hint, call. = FALSE)
    }

    stop(sprintf("'%s' requires package '%s'. %s", api, pkg, install_hint), call. = FALSE)
}

#' @keywords internal
.gflow.get_namespace_export <- function(pkg, name, api = NULL, install_hint = NULL) {
    .gflow.require.namespace(pkg, api = api, install_hint = install_hint)
    getExportedValue(pkg, name)
}

#' @keywords internal
.gflow.require.malo <- function(api = NULL) {
    .gflow.require.namespace(
        pkg = "malo",
        api = api,
        install_hint = "Install it first (e.g. remotes::install_github('pgajer/malo'))."
    )
}
