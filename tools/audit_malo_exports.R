#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
    stop("Usage: Rscript tools/audit_malo_exports.R [path_to_malo_repo]")
}

parse_exports <- function(namespace_path) {
    lines <- readLines(namespace_path, warn = FALSE)
    exports <- grep("^export\\(", lines, value = TRUE)
    sub("^export\\(([^)]+)\\)$", "\\1", exports)
}

parse_forwarder_targets <- function(forwarders_path) {
    lines <- readLines(forwarders_path, warn = FALSE)
    targets <- grep("malo::[A-Za-z0-9_.]+\\(", lines, value = TRUE)
    targets <- sub(".*malo::([A-Za-z0-9_.]+)\\(.*", "\\1", targets)
    sort(unique(targets))
}

gflow_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
default_malo_dir <- file.path(dirname(gflow_dir), "malo")
malo_dir_arg <- if (length(args) == 1L) args[[1L]] else Sys.getenv("MALO_DIR", default_malo_dir)
malo_dir <- normalizePath(malo_dir_arg, winslash = "/", mustWork = TRUE)

gflow_namespace_path <- file.path(gflow_dir, "NAMESPACE")
gflow_forwarders_path <- file.path(gflow_dir, "R", "malo_1d_forwarders.R")
malo_namespace_path <- file.path(malo_dir, "NAMESPACE")

if (!file.exists(gflow_namespace_path)) {
    stop("Missing gflow NAMESPACE at: ", gflow_namespace_path)
}
if (!file.exists(malo_namespace_path)) {
    stop("Missing malo NAMESPACE at: ", malo_namespace_path)
}

gflow_exports <- parse_exports(gflow_namespace_path)
malo_exports <- parse_exports(malo_namespace_path)
forwarder_targets <- if (file.exists(gflow_forwarders_path)) {
    parse_forwarder_targets(gflow_forwarders_path)
} else {
    character(0)
}

migrated_mab_family <- sort(grep("^get\\..*MAB", malo_exports, value = TRUE))

missing_forwarders_in_malo <- setdiff(forwarder_targets, malo_exports)
missing_forwarders_in_gflow <- setdiff(forwarder_targets, gflow_exports)
migrated_mab_still_in_gflow <- intersect(migrated_mab_family, gflow_exports)

cat("malo export audit\n")
cat(sprintf("- gflow dir: %s\n", gflow_dir))
cat(sprintf("- malo dir: %s\n", malo_dir))
cat(sprintf("- forwarder targets found: %d\n", length(forwarder_targets)))
cat(sprintf("- malo MAB-family exports found: %d\n", length(migrated_mab_family)))
if (!file.exists(gflow_forwarders_path)) {
    cat("- note: gflow forwarder file not present (hard-cut mode).\n")
}

ok <- TRUE

if (length(missing_forwarders_in_malo) > 0L) {
    ok <- FALSE
    cat("- ERROR: forwarders missing in malo exports:\n")
    for (nm in missing_forwarders_in_malo) {
        cat(sprintf("  * %s\n", nm))
    }
}

if (length(missing_forwarders_in_gflow) > 0L) {
    ok <- FALSE
    cat("- ERROR: forwarder targets not exported by gflow:\n")
    for (nm in missing_forwarders_in_gflow) {
        cat(sprintf("  * %s\n", nm))
    }
}

if (length(migrated_mab_still_in_gflow) > 0L) {
    ok <- FALSE
    cat("- ERROR: migrated MAB-family symbols still exported by gflow:\n")
    for (nm in migrated_mab_still_in_gflow) {
        cat(sprintf("  * %s\n", nm))
    }
}

if (!ok) {
    quit(save = "no", status = 1L)
}

if (length(forwarder_targets) > 0L) {
    cat("- OK: all gflow->malo forwarders are present in malo exports.\n")
} else {
    cat("- OK: no gflow 1D forwarders detected.\n")
}
cat("- OK: migrated MAB-family symbols are malo-only.\n")
quit(save = "no", status = 0L)
