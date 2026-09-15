#!/usr/bin/env Rscript
# Resolve catalog links from generated aliases; safe to rerun after roxygen.
files <- list.files("man", "[.]Rd$", full.names = TRUE)
alias.topic <- unlist(lapply(files, function(file) {
    rd <- tools::parse_Rd(file)
    aliases <- vapply(Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), rd),
        function(x) paste(unlist(x), collapse = ""), character(1))
    setNames(rep(sub("[.]Rd$", "", basename(file)), length(aliases)), aliases)
}))
path <- "vignettes/function-guide.Rmd"
lines <- readLines(path)
for (i in grep("^[|]", lines)) {
    line <- gsub("\\[(`[^`]+`)\\]\\([^)]*\\)", "\\1", lines[i])
    tokens <- unique(regmatches(line, gregexpr("`[[:alnum:]_.]+[(][)]`", line))[[1L]])
    for (token in tokens) {
        name <- sub("[(][)]`$", "", substring(token, 2L))
        topic <- unname(alias.topic[name])
        if (length(topic) && !is.na(topic)) line <- gsub(token, paste0("[",token,"](",
            "https://pgajer.github.io/gflow/reference/",topic,".html)"),line,fixed=TRUE)
    }
    lines[i] <- line
}
writeLines(lines, path)
