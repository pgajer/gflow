#!/usr/bin/env Rscript
# Build against an installed artifact with the current generated code manifest.
if (!requireNamespace("pkgdown", quietly = TRUE)) stop("Install pkgdown to build the website.")
installed <- system.file("extdata", "gflow-code-manifest.tsv", package = "gflow")
if (!nzchar(installed) || !identical(unname(tools::md5sum(installed)),
                                    unname(tools::md5sum("inst/extdata/gflow-code-manifest.tsv")))) {
    stop("Install the current Makefile-built archive in the selected R library before building the website.")
}
# Reuse the guide's task categories, and resolve multiple aliases to one help page.
db <- tools::Rd_db("gflow")
aliases <- lapply(db, function(rd) vapply(Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), rd),
                                        function(x) paste(unlist(x), collapse = ""), character(1)))
alias.topic <- setNames(rep(sub("[.]Rd$", "", names(db)), lengths(aliases)), unlist(aliases))
lines <- readLines("vignettes/function-guide.Rmd")
lines <- lines[(match("<!-- export-catalog-start -->", lines) + 1L):(match("<!-- export-catalog-end -->", lines) - 1L)]
groups <- list(); assigned <- character(); heading <- "Getting started"
for (line in lines) {
    if (grepl("^### ", line)) heading <- sub("^### ", "", line)
    if (!grepl("^\\| `[^`]+\\(\\)` ", line)) next
    name <- sub("^\\| `([^`]+)\\(\\)`.*", "\\1", line)
    topic <- unname(alias.topic[name])
    if (is.na(topic) || topic %in% assigned) next
    groups[[heading]] <- c(groups[[heading]], topic)
    assigned <- c(assigned, topic)
}
for (topic in c("gflow-package", "gflow-migration", "summary.basin_complex", "plot.basin_complex")) {
    if (topic %in% sub("[.]Rd$", "", names(db)) && !topic %in% assigned) {
        groups[["Package help and basin displays"]] <- c(groups[["Package help and basin displays"]], topic)
        assigned <- c(assigned, topic)
    }
}
reference <- lapply(names(groups), function(name) list(title = name, contents = groups[[name]]))
topics <- pkgdown::as_pkgdown(".")$topics
additional <- setdiff(topics$name[!topics$internal], assigned)
if (length(additional)) {
    reference <- c(reference, list(list(
        title = "Additional object methods and historical helpers",
        desc = paste("Existing method and helper documentation is retained here.",
                     "For supported public entry points, use the task sections above."),
        contents = additional
    )))
}
pkgdown::build_site(install = FALSE, new_process = FALSE, examples = FALSE,
                    preview = FALSE, override = list(reference = reference))
# Evaluate the small supported entry-point examples; retain source-only examples
# on specialized/archive reference pages that require saved or optional inputs.
pkgdown::build_reference(topics = c("gflow-package", "lcor", "create.basin.complex"),
                         examples = TRUE, lazy = FALSE, devel = FALSE,
                         override = list(reference = reference))
# Rd example plots have no caption channel; describe this known entry-point plot.
page <- "build/site/reference/create.basin.complex.html"
html <- xml2::read_html(page)
plots <- xml2::xml_find_all(html, ".//img[not(@alt) or @alt='']")
stopifnot(length(plots) == 1L)
xml2::xml_set_attr(plots, "alt", paste(
    "Merge tree for the seven-vertex path: the plateau peak at vertices 2–3",
    "joins the higher peak at vertex 6 at field height 1."
))
xml2::write_html(html, page)
dir.create("build/site/vignettes/figures", recursive = TRUE, showWarnings = FALSE)
file.copy("vignettes/figures/basin-showcase.png", "build/site/vignettes/figures", overwrite = TRUE)
cat("Website built at build/site/index.html from the installed package.\n")
