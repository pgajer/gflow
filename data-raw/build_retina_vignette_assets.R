pkg.root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate package root containing DESCRIPTION.")
    }
    path <- parent
  }
}

root <- pkg.root()
source(file.path(root, "inst", "scripts", "build_retina_vignette_assets.R"))
