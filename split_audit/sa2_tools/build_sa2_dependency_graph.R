#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run from the gflow repository root.")
}

audit_dir <- file.path(repo_root, "split_audit")
out_dir <- file.path(audit_dir, "sa2_dependency_graph")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_lines <- function(path) readLines(file.path(repo_root, path), warn = FALSE)

list_r_files <- function() {
  file.path("R", list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = FALSE))
}

list_src_files <- function() {
  file.path("src", list.files(file.path(repo_root, "src"), pattern = "\\.(c|cpp|h|hpp)$", full.names = FALSE))
}

extract_defs <- function(path) {
  txt <- read_lines(path)
  rx <- "^([A-Za-z.][A-Za-z0-9._]*)\\s*(<-|=)\\s*function\\s*\\(.*$"
  hit <- grep(rx, txt, value = TRUE)
  if (!length(hit)) {
    return(data.frame(function_name = character(), file = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    function_name = sub(rx, "\\1", hit),
    file = path,
    stringsAsFactors = FALSE
  )
}

call_name <- function(expr) {
  if (!is.call(expr)) return(NA_character_)
  f <- expr[[1L]]
  if (is.symbol(f)) return(as.character(f))
  if (is.call(f) && identical(as.character(f[[1L]]), "::")) {
    return(paste(as.character(f[[2L]]), as.character(f[[3L]]), sep = "::"))
  }
  if (is.call(f) && identical(as.character(f[[1L]]), ":::")) {
    return(paste(as.character(f[[2L]]), as.character(f[[3L]]), sep = ":::"))
  }
  NA_character_
}

collect_calls <- function(expr) {
  out <- character()
  walk <- function(x) {
    if (is.call(x)) {
      nm <- call_name(x)
      if (!is.na(nm)) out <<- c(out, nm)
      for (ii in seq_along(x)[-1L]) walk(x[[ii]])
    } else if (is.expression(x) || is.pairlist(x)) {
      for (ii in seq_along(x)) walk(x[[ii]])
    }
  }
  walk(expr)
  out
}

extract_calls <- function(path) {
  parsed <- tryCatch(parse(file.path(repo_root, path), keep.source = FALSE),
                     error = function(e) expression())
  calls <- collect_calls(parsed)
  calls <- calls[!is.na(calls)]
  if (!length(calls)) {
    return(data.frame(file = character(), call = character(), count = integer(), stringsAsFactors = FALSE))
  }
  tab <- sort(table(calls), decreasing = TRUE)
  data.frame(
    file = path,
    call = names(tab),
    count = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

extract_dot_call_symbols <- function(path) {
  txt <- read_lines(path)
  idx <- grep("\\.Call\\s*\\(", txt)
  idx <- idx[!grepl("^\\s*#", txt[idx])]
  if (!length(idx)) {
    return(data.frame(file = character(), symbol = character(), raw = character(), stringsAsFactors = FALSE))
  }
  rows <- lapply(idx, function(ii) {
    window <- txt[seq.int(ii, min(length(txt), ii + 4L))]
    joined <- paste(window, collapse = " ")
    sym <- NA_character_
    quoted <- regexec("\\.Call\\s*\\(\\s*[`'\"]?([A-Za-z0-9_\\.]+)[`'\"]?", joined)
    m <- regmatches(joined, quoted)[[1L]]
    if (length(m) >= 2L && !identical(m[[2L]], "PACKAGE")) {
      sym <- m[[2L]]
    } else {
      next_line <- trimws(window[min(2L, length(window))])
      bare <- regexec("^([A-Za-z0-9_\\.]+)\\s*,", next_line)
      bm <- regmatches(next_line, bare)[[1L]]
      if (length(bm) >= 2L) sym <- bm[[2L]]
    }
    data.frame(
      file = path,
      symbol = ifelse(is.na(sym), "UNRESOLVED", sym),
      raw = paste(window, collapse = " "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

classify_call <- function(call, def_map) {
  syntax <- c("{", "(", "if", "for", "while", "function", "<-", "=", "::",
              ":::", "$", "[", "[[", "@", "return", "next", "break",
              "+", "-", "*", "/", "^", "%%", "%/%", "%in%", ":", "~",
              "!", "&&", "||", "&", "|", "==", "!=", "<", "<=", ">", ">=")
  if (call %in% syntax) return(c(package = "syntax", target_file = NA_character_))
  if (call %in% def_map$function_name) {
    files <- paste(unique(def_map$file[def_map$function_name == call]), collapse = ";")
    return(c(package = "gflow_internal", target_file = files))
  }
  if (grepl("::", call, fixed = TRUE)) {
    return(c(package = sub("::.*", "", call), target_file = NA_character_))
  }
  c(package = "external_or_base", target_file = NA_character_)
}

inventory_path <- file.path(audit_dir, "gflow_split_inventory.csv")
if (!file.exists(inventory_path)) {
  stop("Missing SA1 inventory. Run split_audit/scripts/build_sa1_inventory.R first.")
}
inventory <- read.csv(inventory_path, stringsAsFactors = FALSE)
geosmooth_r_files <- inventory$asset[inventory$asset_type == "R" &
                                       inventory$proposed_package == "geosmooth"]
pilot_files <- c("R/kernel_local_polynomial_cv.R")

defs <- do.call(rbind, lapply(list_r_files(), extract_defs))
defs <- defs[order(defs$file, defs$function_name), ]
write.csv(defs, file.path(out_dir, "sa2_all_function_defs.csv"), row.names = FALSE)

target_calls <- do.call(rbind, lapply(geosmooth_r_files, extract_calls))
if (nrow(target_calls)) {
  classes <- t(vapply(target_calls$call, classify_call, character(2), def_map = defs))
  target_calls$call_package <- classes[, "package"]
  target_calls$target_file <- classes[, "target_file"]
  target_calls$target_proposed_package <- NA_character_
  for (ii in seq_len(nrow(target_calls))) {
    tf <- target_calls$target_file[[ii]]
    if (!is.na(tf) && nzchar(tf)) {
      parts <- strsplit(tf, ";", fixed = TRUE)[[1L]]
      pkgs <- inventory$proposed_package[match(parts, inventory$asset)]
      pkgs <- unique(pkgs[!is.na(pkgs)])
      target_calls$target_proposed_package[[ii]] <- paste(pkgs, collapse = ";")
    }
  }
}
write.csv(target_calls, file.path(out_dir, "sa2_geosmooth_r_dependency_edges.csv"), row.names = FALSE)

pilot_calls <- target_calls[target_calls$file %in% pilot_files, , drop = FALSE]
write.csv(pilot_calls, file.path(out_dir, "sa2_pilot_kernel_local_polynomial_edges.csv"), row.names = FALSE)

native_calls <- do.call(rbind, lapply(list_r_files(), extract_dot_call_symbols))
write.csv(native_calls, file.path(out_dir, "sa2_all_native_symbol_calls.csv"), row.names = FALSE)

geosmooth_native <- native_calls[native_calls$file %in% geosmooth_r_files, , drop = FALSE]
write.csv(geosmooth_native, file.path(out_dir, "sa2_geosmooth_native_symbol_calls.csv"), row.names = FALSE)

rcpp_wrappers <- defs[defs$file == "R/RcppExports.R", , drop = FALSE]
geosmooth_rcpp_wrapper_calls <- target_calls[
  target_calls$call %in% rcpp_wrappers$function_name,
  c("file", "call", "count", "target_file"),
  drop = FALSE
]
write.csv(
  geosmooth_rcpp_wrapper_calls,
  file.path(out_dir, "sa2_geosmooth_rcpp_wrapper_calls.csv"),
  row.names = FALSE
)

src_text <- unlist(lapply(list_src_files(), function(path) {
  paste(path, read_lines(path), sep = "\t")
}), use.names = FALSE)
symbols <- unique(c(native_calls$symbol, paste0("_gflow_", rcpp_wrappers$function_name)))
symbols <- symbols[nzchar(symbols) & !is.na(symbols)]
src_hits <- do.call(rbind, lapply(symbols, function(sym) {
  hits <- grep(sym, src_text, fixed = TRUE, value = TRUE)
  if (!length(hits)) return(NULL)
  data.frame(
    symbol = sym,
    source_hit = hits,
    stringsAsFactors = FALSE
  )
}))
if (is.null(src_hits)) {
  src_hits <- data.frame(symbol = character(), source_hit = character())
}
write.csv(src_hits, file.path(out_dir, "sa2_native_symbol_source_hits.csv"), row.names = FALSE)

dep_summary <- aggregate(
  count ~ file + target_proposed_package,
  target_calls[target_calls$call_package == "gflow_internal" &
                 nzchar(target_calls$target_proposed_package), ],
  sum
)
dep_summary <- dep_summary[order(dep_summary$file, dep_summary$target_proposed_package), ]
write.csv(dep_summary, file.path(out_dir, "sa2_geosmooth_internal_dependency_summary.csv"), row.names = FALSE)

top_external <- target_calls[target_calls$call_package == "external_or_base", ]
top_external <- aggregate(count ~ file + call, top_external, sum)
top_external <- top_external[order(top_external$file, -top_external$count, top_external$call), ]
write.csv(top_external, file.path(out_dir, "sa2_geosmooth_external_or_base_calls.csv"), row.names = FALSE)

cat("SA2 dependency graph written to ", out_dir, "\n", sep = "")
cat("Geosmooth R files: ", length(geosmooth_r_files), "\n", sep = "")
cat("R-level geosmooth call rows: ", nrow(target_calls), "\n", sep = "")
cat("Pilot call rows: ", nrow(pilot_calls), "\n", sep = "")
cat("Geosmooth native .Call rows: ", nrow(geosmooth_native), "\n", sep = "")
cat("Geosmooth Rcpp wrapper call rows: ", nrow(geosmooth_rcpp_wrapper_calls), "\n", sep = "")
