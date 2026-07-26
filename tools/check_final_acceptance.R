#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
script.path <- if (length(file.arg)) {
    sub("^--file=", "", file.arg[[1L]])
} else {
    normalizePath("tools/check_final_acceptance.R")
}
root <- normalizePath(file.path(dirname(script.path), ".."))

source(file.path(root, "tools/cleanup_ledger_lib.R"))
source(file.path(root, "tools/s3_namespace_lib.R"))

ledger <- read.csv(
    file.path(root, "split_audit/cleanup/api-ownership.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
dependencies <- read.csv(
    file.path(root, "split_audit/cleanup/dependency-ownership.csv"),
    stringsAsFactors = FALSE
)
native <- read.csv(
    file.path(root, "split_audit/cleanup/native-symbol-ownership.csv"),
    stringsAsFactors = FALSE
)
errors <- c(
    cleanup.validate.ledger(root, ledger),
    cleanup.check.protected.surface(root),
    s3.validate.namespace(root)
)

exports <- ledger[ledger$item_type == "export", , drop = FALSE]
allowed.public <- c("CORE-ANALYSIS", "CORE-ASSOCIATION", "PROTECTED")
bad.public <- exports$item_name[!exports$ownership %in% allowed.public]
if (length(bad.public)) {
    errors <- c(errors, paste(
        "Public functions have non-core ownership:",
        paste(bad.public, collapse = ", ")
    ))
}

in.scope <- ledger$ownership != "PROTECTED"
missing.docs <- ledger$item_name[
    in.scope & ledger$documentation == "none"
]
missing.tests <- ledger$item_name[
    in.scope & ledger$tests == "none"
]
if (length(missing.docs)) {
    errors <- c(errors, paste(
        "In-scope namespace items lack a documentation owner:",
        paste(missing.docs, collapse = ", ")
    ))
}
if (length(missing.tests)) {
    errors <- c(errors, paste(
        "In-scope namespace items lack a test owner:",
        paste(missing.tests, collapse = ", ")
    ))
}

description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, ]
declared <- paste(description[c(
    "Depends", "Imports", "Suggests", "LinkingTo", "VignetteBuilder"
)], collapse = ",")
banned.dependencies <- c(
    "gflowx", "transport", "smacof", "geometry", "plotly", "shiny",
    "clValid", "fpc"
)
present.dependencies <- banned.dependencies[vapply(
    banned.dependencies,
    function(package) grepl(
        paste0("(^|[,[:space:]])", package, "([,[:space:]()]|$)"),
        declared
    ),
    logical(1)
)]
if (length(present.dependencies)) {
    errors <- c(errors, paste(
        "Retired dependencies remain declared:",
        paste(present.dependencies, collapse = ", ")
    ))
}

retired.files <- c(
    "R/clustering.R",
    "R/diffusion_pseudotime_sparse.R",
    "R/phate_core.R",
    "R/quadform_geodesics.R",
    "R/graph_endpoint_geometry.R",
    "R/select_pts.R",
    "R/fit_rho_randomwalk.R",
    "R/distance_quantile_bin_analysis.R",
    "R/geodesics.R",
    "src/diffusion_pseudotime_sparse_r.cpp",
    "src/phate_core_r.cpp",
    "src/quadform_delaunay_edges.cpp",
    "src/quadform_grid_geodesics.cpp",
    "src/graph_endpoint_geometry_rcpp.cpp"
)
present.files <- retired.files[file.exists(file.path(root, retired.files))]
if (length(present.files)) {
    errors <- c(errors, paste(
        "Retired implementation files remain:",
        paste(present.files, collapse = ", ")
    ))
}

required.story <- c(
    "README.md",
    "REFERENCE.md",
    "NEWS.md",
    "man/gflow-package.Rd",
    "split_audit/cleanup/breaking-release-migration.md"
)
missing.story <- required.story[!file.exists(file.path(root, required.story))]
if (length(missing.story)) {
    errors <- c(errors, paste(
        "Public package story is incomplete:",
        paste(missing.story, collapse = ", ")
    ))
}

if (!identical(unname(description[["Version"]]), "0.2.0")) {
    errors <- c(errors, "Breaking release version must be 0.2.0.")
}
if (any(!nzchar(dependencies$ownership)) ||
    any(!nzchar(native$ownership)) ||
    any(native$r_callers == "none") ||
    any(native$native_sources == "none")) {
    errors <- c(
        errors,
        "Dependency/native ownership contains blank or uncalled entries."
    )
}

errors <- unique(errors)
if (length(errors)) {
    cat(paste0("ERROR: ", errors, "\n"), sep = "")
    quit(status = 1L)
}

cat(
    "Final acceptance passes for ",
    nrow(exports), " exports (",
    sum(exports$ownership == "PROTECTED"), " protected), ",
    sum(ledger$item_type == "s3"), " S3 methods, ",
    nrow(dependencies), " dependency declarations, and ",
    nrow(native), " native registrations.\n",
    sep = ""
)
