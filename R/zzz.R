# internal package state
if (!exists(".gflow_state", inherits = FALSE)) {
  .gflow_state <- new.env(parent = emptyenv())
}

# Minimal lazy S3 registration helper (no forced load of the other pkg)
.s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1L)
  parts <- strsplit(generic, "::", fixed = TRUE)[[1]]
  if (length(parts) == 2L) {
    package <- parts[[1]]
    gen     <- parts[[2]]
  } else {
    package <- NULL
    gen     <- parts[[1]]
  }

  if (is.null(method)) {
    # look up a method named e.g. plot3d.class in the caller's env
    method <- get(paste0(gen, ".", class), envir = parent.frame())
  }

  register <- function(pkg) {
    ns <- asNamespace(pkg)
    if (exists(gen, envir = ns, inherits = FALSE)) {
      registerS3method(gen, class, method, envir = ns)
    }
  }

  if (is.null(package)) {
    # base-only generics
    registerS3method(gen, class, method)
  } else {
    if (isNamespaceLoaded(package)) register(package)
    setHook(packageEvent(package, "onLoad"), function(...) register(package))
  }
}

.onLoad <- function(libname, pkgname) {
  # 1) Headless safety: only set rgl.useNULL if user hasn't set it
  .gflow_state$old_rgl_useNULL <- getOption("rgl.useNULL", NULL)
  headless <- (!interactive()) ||
              Sys.getenv("DISPLAY") == "" ||
              identical(Sys.getenv("RGL_USE_NULL"), "TRUE")
  if (headless && is.null(.gflow_state$old_rgl_useNULL)) {
    options(rgl.useNULL = TRUE)
  }
}

.onUnload <- function(libpath) {
  # restore user's prior rgl.useNULL option
  if (is.null(.gflow_state$old_rgl_useNULL)) {
    options(rgl.useNULL = NULL)   # remove the option if we introduced it
  } else {
    options(rgl.useNULL = .gflow_state$old_rgl_useNULL)
  }
}

is_valgrind <- function() nzchar(Sys.getenv("R_VALGRIND_OPTS"))
is_check    <- function() identical(Sys.getenv("_R_CHECK_PACKAGE_NAME_"), "gflow") ||
                           nzchar(Sys.getenv("_R_CHECK_LENGTH_1_LOGIC2_"))          # common CRAN flag
safe_cores  <- function(n) if (is_valgrind() || is_check()) 1L else as.integer(n)
