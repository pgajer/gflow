args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("^--file=", args, value = TRUE)
script.path <- if (length(file.arg)) {
    sub("^--file=", "", file.arg[[1L]])
} else {
    normalizePath("tools/build_phase7_ownership.R")
}
root <- normalizePath(file.path(dirname(script.path), ".."))

split.dependencies <- function(value) {
    if (is.na(value) || !nzchar(value)) {
        return(character())
    }
    trimws(sub("\\s*\\(.*$", "", strsplit(value, ",", fixed = TRUE)[[1L]]))
}

description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, ]
dependency.fields <- c("Depends", "Imports", "Suggests", "LinkingTo")
declared <- do.call(
    rbind,
    lapply(dependency.fields, function(field) {
        data.frame(
            field = field,
            dependency = split.dependencies(description[[field]]),
            stringsAsFactors = FALSE
        )
    })
)
declared <- declared[nzchar(declared$dependency), , drop = FALSE]
declared <- rbind(
    declared,
    data.frame(
        field = "VignetteBuilder",
        dependency = split.dependencies(description[["VignetteBuilder"]]),
        stringsAsFactors = FALSE
    ),
    data.frame(
        field = "SystemRequirements",
        dependency = c("C++17", "GNU make", "OpenMP"),
        stringsAsFactors = FALSE
    )
)

dependency.metadata <- list(
    "Depends:R" = c("CORE-PRIVATE", "R runtime", "required"),
    "Imports:Rcpp" = c("CORE-PRIVATE", "registered native interface", "required"),
    "Imports:stats" = c("PROTECTED", "core statistics and protected construction", "required"),
    "Imports:graphics" = c("PROTECTED", "core and protected plotting", "required"),
    "Imports:grDevices" = c("PROTECTED", "core and protected graphics devices", "required"),
    "Imports:FNN" = c("PROTECTED", "nearest-neighbor support for protected construction", "required"),
    "Imports:doParallel" = c("DOMAIN", "legacy private clustering", "required-until-phase-8"),
    "Imports:foreach" = c("DOMAIN", "legacy private clustering", "required-until-phase-8"),
    "Imports:parallel" = c("CORE-ASSOCIATION", "parallel local association", "required"),
    "Imports:dgraphs" = c("PROTECTED", "graph substrate for construction and analysis", "required"),
    "Imports:igraph" = c("PROTECTED", "protected and core graph visualization", "required"),
    "Imports:utils" = c("PROTECTED", "package utilities used by protected and core code", "required"),
    "Suggests:rgl" = c("PROTECTED", "optional protected 3D visualization", "optional"),
    "Suggests:htmlwidgets" = c("CORE-ANALYSIS", "optional interactive core vignette output", "optional"),
    "Suggests:htmltools" = c("CORE-PRIVATE", "optional widget rendering support", "optional"),
    "Suggests:plotly" = c("CORE-PRIVATE", "legacy interactive selection", "remove-phase-8"),
    "Suggests:shiny" = c("CORE-PRIVATE", "legacy interactive selection", "remove-phase-8"),
    "Suggests:smacof" = c("DGRAPHS", "legacy PHATE embedding backend", "remove-phase-8"),
    "Suggests:grip" = c("EXAMPLE", "optional layout in core vignette", "optional"),
    "Suggests:dbscan" = c("CORE-ANALYSIS", "optional extrema clustering diagnostics", "optional"),
    "Suggests:clValid" = c("DOMAIN", "legacy private clustering validation", "remove-phase-8"),
    "Suggests:mclust" = c("CORE-PRIVATE", "optional private plotting palettes", "optional"),
    "Suggests:MASS" = c("DOMAIN", "legacy trajectory-report robust fit", "optional"),
    "Suggests:RColorBrewer" = c("DOMAIN", "legacy report palettes", "optional"),
    "Suggests:viridis" = c("DOMAIN", "legacy MGCP visualization", "optional"),
    "Suggests:Matrix" = c("DGRAPHS", "sparse diffusion and stability support", "optional"),
    "Suggests:arm" = c("DOMAIN", "legacy disk-association report", "optional"),
    "Suggests:brglm2" = c("DOMAIN", "legacy disk-association report", "optional"),
    "Suggests:circlize" = c("CORE-ASSOCIATION", "optional lcor heatmap colors", "optional"),
    "Suggests:cluster" = c("DOMAIN", "legacy trajectory and clustering analysis", "optional"),
    "Suggests:ComplexHeatmap" = c("CORE-ASSOCIATION", "optional lcor heatmap", "optional"),
    "Suggests:dtw" = c("DOMAIN", "legacy trajectory clustering", "optional"),
    "Suggests:fpc" = c("DOMAIN", "legacy private clustering validation", "remove-phase-8"),
    "Suggests:geometry" = c("DGRAPHS", "legacy quadform geodesic backend", "remove-phase-8"),
    "Suggests:ggplot2" = c("PROTECTED", "optional protected flow plotting", "optional"),
    "Suggests:hexbin" = c("DOMAIN", "legacy distance-diagnostic plot", "optional"),
    "Suggests:logistf" = c("DOMAIN", "legacy disk-association report", "optional"),
    "Suggests:proxy" = c("DOMAIN", "legacy trajectory distance", "optional"),
    "Suggests:scales" = c("DOMAIN", "legacy concordance visualization", "optional"),
    "Suggests:testthat" = c("CORE-PRIVATE", "package tests", "development"),
    "Suggests:knitr" = c("CORE-PRIVATE", "vignette build", "development"),
    "Suggests:rmarkdown" = c("CORE-PRIVATE", "vignette build", "development"),
    "Suggests:rstan" = c("CORE-PRIVATE", "legacy private statistical helper", "optional"),
    "Suggests:HDInterval" = c("CORE-PRIVATE", "legacy private interval helper", "optional"),
    "Suggests:rootSolve" = c("CORE-PRIVATE", "legacy private root helper", "optional"),
    "Suggests:lme4" = c("DOMAIN", "legacy two-factor analysis", "optional"),
    "Suggests:clue" = c("DGRAPHS", "legacy private graph matching", "remove-phase-8"),
    "Suggests:gflowx" = c("GFLOWX", "legacy endpoint-smoothing adapter", "remove-phase-8"),
    "LinkingTo:Rcpp" = c("CORE-PRIVATE", "compiled interface headers", "required"),
    "VignetteBuilder:knitr" = c("CORE-PRIVATE", "vignette build", "development"),
    "SystemRequirements:C++17" = c("CORE-PRIVATE", "compiled package baseline", "required"),
    "SystemRequirements:GNU make" = c("CORE-PRIVATE", "native build orchestration", "required"),
    "SystemRequirements:OpenMP" = c("CORE-PRIVATE", "optional native parallelism", "optional")
)

keys <- paste(declared$field, declared$dependency, sep = ":")
missing.metadata <- setdiff(keys, names(dependency.metadata))
if (length(missing.metadata)) {
    stop(
        "Missing dependency ownership metadata: ",
        paste(missing.metadata, collapse = ", ")
    )
}

scan.files <- c(
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$", full.names = TRUE),
    list.files(
        file.path(root, "src"),
        pattern = "\\.(c|cc|cpp|cxx|h|hpp)$",
        recursive = TRUE,
        full.names = TRUE
    ),
    list.files(
        file.path(root, "tests"),
        pattern = "\\.[Rr]$",
        recursive = TRUE,
        full.names = TRUE
    ),
    list.files(
        file.path(root, "vignettes"),
        pattern = "\\.(Rmd|qmd)$",
        recursive = TRUE,
        full.names = TRUE
    )
)
scan.text <- lapply(scan.files, readLines, warn = FALSE)
relative <- substring(scan.files, nchar(root) + 2L)

dependency.evidence <- function(package) {
    hits <- relative[
        vapply(
            scan.text,
            function(lines) any(grepl(package, lines, fixed = TRUE)),
            logical(1L)
        )
    ]
    hits <- sort(unique(hits))
    if (length(hits) > 12L) {
        hits <- c(hits[seq_len(12L)], paste0("+", length(hits) - 12L, " files"))
    }
    if (!length(hits)) "DESCRIPTION-only" else paste(hits, collapse = ";")
}

dependency.rows <- lapply(seq_len(nrow(declared)), function(i) {
    metadata <- dependency.metadata[[keys[[i]]]]
    data.frame(
        key = keys[[i]],
        field = declared$field[[i]],
        dependency = declared$dependency[[i]],
        ownership = metadata[[1L]],
        feature = metadata[[2L]],
        policy = metadata[[3L]],
        evidence = dependency.evidence(declared$dependency[[i]]),
        stringsAsFactors = FALSE
    )
})
dependency.ledger <- do.call(rbind, dependency.rows)
write.csv(
    dependency.ledger,
    file.path(root, "split_audit/cleanup/dependency-ownership.csv"),
    row.names = FALSE,
    quote = TRUE
)

init.path <- file.path(root, "src/init.c")
init.lines <- readLines(init.path)
registration.lines <- grep(
    "^[[:space:]]*\\{\"[^\"]+\"",
    init.lines
)
registrations <- sub(
    "^[[:space:]]*\\{\"([^\"]+)\".*",
    "\\1",
    init.lines[registration.lines]
)
if (anyDuplicated(registrations)) {
    stop(
        "Duplicate native registrations: ",
        paste(unique(registrations[duplicated(registrations)]), collapse = ", ")
    )
}

r.files <- list.files(
    file.path(root, "R"),
    pattern = "\\.[Rr]$",
    full.names = TRUE
)
r.text <- lapply(r.files, readLines, warn = FALSE)
r.relative <- substring(r.files, nchar(root) + 2L)
api <- read.csv(
    file.path(root, "split_audit/cleanup/api-ownership.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
protected.files <- sub(
    "^FILE\\t([^\\t]+)\\t.*$",
    "\\1",
    grep(
        "^FILE\\t",
        readLines(file.path(root, "split_audit/cleanup/protected-basin-surface.txt")),
        value = TRUE
    )
)

owner.priority <- c(
    PROTECTED = 1L,
    `CORE-ASSOCIATION` = 2L,
    `CORE-ANALYSIS` = 3L,
    DGRAPHS = 4L,
    GFLOWX = 5L,
    DOMAIN = 6L,
    EXAMPLE = 7L,
    `CORE-PRIVATE` = 8L
)

native.owner <- function(files, symbol) {
    if (identical(symbol, "S_gflow_openmp_diag")) {
        return("CORE-PRIVATE")
    }
    if (any(files %in% protected.files)) {
        return("PROTECTED")
    }
    owners <- api$ownership[api$source_file %in% files]
    owners <- owners[owners %in% names(owner.priority)]
    if (!length(owners)) {
        return("CORE-PRIVATE")
    }
    owners[[which.min(owner.priority[owners])]]
}

src.files <- list.files(
    file.path(root, "src"),
    pattern = "\\.(c|cc|cpp|cxx)$",
    recursive = TRUE,
    full.names = TRUE
)
src.files <- setdiff(src.files, init.path)
src.text <- lapply(src.files, readLines, warn = FALSE)
src.relative <- substring(src.files, nchar(root) + 2L)

native.rows <- lapply(seq_along(registrations), function(i) {
    symbol <- registrations[[i]]
    callers <- r.relative[
        vapply(r.text, function(lines) any(grepl(symbol, lines, fixed = TRUE)), logical(1L))
    ]
    if (identical(symbol, "S_gflow_openmp_diag")) {
        callers <- "tests/testthat/test-phase7-native-dependencies.R"
    }
    sources <- src.relative[
        vapply(src.text, function(lines) any(grepl(symbol, lines, fixed = TRUE)), logical(1L))
    ]
    method <- if (
        registration.lines[[i]] <
            grep("static const R_CallMethodDef CallMethods", init.lines)[[1L]]
    ) "C" else "Call"
    data.frame(
        symbol = symbol,
        interface = method,
        ownership = native.owner(callers, symbol),
        r_callers = if (length(callers)) paste(callers, collapse = ";") else "none",
        native_sources = if (length(sources)) paste(sources, collapse = ";") else "none",
        status = "retained",
        stringsAsFactors = FALSE
    )
})
native.ledger <- do.call(rbind, native.rows)
write.csv(
    native.ledger,
    file.path(root, "split_audit/cleanup/native-symbol-ownership.csv"),
    row.names = FALSE,
    quote = TRUE
)

cat(
    "Wrote", nrow(dependency.ledger), "dependency declarations and",
    nrow(native.ledger), "native registrations.\n"
)
