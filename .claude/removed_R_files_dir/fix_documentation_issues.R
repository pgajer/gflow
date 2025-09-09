#!/usr/bin/env Rscript

# Script to fix R documentation issues for CRAN compliance
# This script fixes:
# 1. Bracket notation issues in roxygen comments
# 2. Non-ASCII characters
# 3. Missing quotes in .C() and .Call() functions

fix_bracket_notation <- function(file_path) {
    lines <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
    modified <- FALSE
    
    for (i in seq_along(lines)) {
        if (grepl("^#'", lines[i])) {
            original_line <- lines[i]
            
            # Fix various bracket notation patterns in roxygen comments
            # [i,j] -> \code{[i,j]} (but avoid already escaped ones)
            new_line <- gsub("([^\\\\]|^)(\\[[^]]*,[^]]*\\])", "\\1\\\\code{\\2}", lines[i])
            
            # Fix double backslashes that might have been created
            new_line <- gsub("\\\\\\\\code\\{", "\\\\code{", new_line)
            
            # Fix specific problematic patterns
            new_line <- gsub("\\\\\\([^)]*\\[i,j\\][^)]*\\\\\\)", "\\\\code{[i,j]}", new_line)
            new_line <- gsub("element \\[i,j\\]", "element \\\\code{[i,j]}", new_line)
            new_line <- gsub("Matrix where \\[i,j\\]", "Matrix where \\\\code{[i,j]}", new_line)
            new_line <- gsub("Matrix where \\\\\\([^)]*\\[i,j\\][^)]*\\\\\\)", "Matrix where \\\\code{[i,j]}", new_line)
            
            # Handle escaped brackets that should be code
            new_line <- gsub("\\\\\\[([^\\]]*,[^\\]]*)\\\\\\]", "\\\\code{[\\1]}", new_line)
            
            if (new_line != original_line) {
                lines[i] <- new_line
                modified <- TRUE
                cat("Fixed bracket notation in", basename(file_path), "line", i, "\n")
                cat("  Before:", original_line, "\n")
                cat("  After: ", new_line, "\n\n")
            }
        }
    }
    
    if (modified) {
        writeLines(lines, file_path, useBytes = TRUE)
        return(TRUE)
    }
    return(FALSE)
}

fix_call_quotes <- function(file_path) {
    content <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
    modified <- FALSE
    
    # Pattern to match .Call or .C with unquoted first argument
    # This is more conservative - only fix obvious cases
    for (i in seq_along(content)) {
        line <- content[i]
        original_line <- line
        
        # Fix .Call("S_function_name" -> .Call("S_function_name"
        # Only fix if clearly unquoted
        if (grepl("\\.Call\\s*\\(\\s*[^\"']\\w", line)) {
            line <- gsub("(\\.Call\\s*\\(\\s*)([A-Za-z_][A-Za-z0-9_]*)", "\\1\"\\2\"", line)
        }
        
        if (grepl("\\.C\\s*\\(\\s*[^\"']\\w", line)) {
            line <- gsub("(\\.C\\s*\\(\\s*)([A-Za-z_][A-Za-z0-9_]*)", "\\1\"\\2\"", line)
        }
        
        if (line != original_line) {
            content[i] <- line
            modified <- TRUE
            cat("Fixed .Call/.C quotes in", basename(file_path), "line", i, "\n")
            cat("  Before:", original_line, "\n")
            cat("  After: ", line, "\n\n")
        }
    }
    
    if (modified) {
        writeLines(content, file_path, useBytes = TRUE)
        return(TRUE)
    }
    return(FALSE)
}

fix_non_ascii <- function(file_path) {
    content <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
    modified <- FALSE
    
    # Common non-ASCII replacements
    replacements <- list(
        "\u2018" = "'", # left single quote
        "\u2019" = "'", # right single quote
        "\u201C" = "\"", # left double quote
        "\u201D" = "\"", # right double quote
        "\u00D7" = "x", # multiplication
        "\u2264" = "<=", # less than or equal
        "\u2265" = ">=", # greater than or equal
        "\u2260" = "!=", # not equal
        "\u03B1" = "\\\\alpha",
        "\u03B2" = "\\\\beta", 
        "\u03B3" = "\\\\gamma",
        "\u03B4" = "\\\\delta",
        "\u03B5" = "\\\\epsilon",
        "\u03B8" = "\\\\theta",
        "\u03BB" = "\\\\lambda",
        "\u03BC" = "\\\\mu",
        "\u03BD" = "\\\\nu",
        "\u03C0" = "\\\\pi",
        "\u03C1" = "\\\\rho",
        "\u03C3" = "\\\\sigma",
        "\u03C4" = "\\\\tau",
        "\u03C6" = "\\\\phi",
        "\u03C7" = "\\\\chi",
        "\u03C8" = "\\\\psi",
        "\u03C9" = "\\\\omega"
    )
    
    for (i in seq_along(content)) {
        line <- content[i]
        original_line <- line
        
        for (char in names(replacements)) {
            if (grepl(char, line, fixed = TRUE)) {
                line <- gsub(char, replacements[[char]], line, fixed = TRUE)
            }
        }
        
        if (line != original_line) {
            content[i] <- line
            modified <- TRUE
            cat("Fixed non-ASCII in", basename(file_path), "line", i, "\n")
            cat("  Before:", original_line, "\n")
            cat("  After: ", line, "\n\n")
        }
    }
    
    if (modified) {
        writeLines(content, file_path, useBytes = TRUE)
        return(TRUE)
    }
    return(FALSE)
}

# Main execution
cat("=== Fixing R Documentation Issues ===\n\n")

# Get all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)

# Track results
bracket_fixes <- 0
call_fixes <- 0 
ascii_fixes <- 0

for (file in r_files) {
    cat("Processing:", basename(file), "\n")
    
    # Fix bracket notation
    if (fix_bracket_notation(file)) {
        bracket_fixes <- bracket_fixes + 1
    }
    
    # Fix .Call/.C quotes
    if (fix_call_quotes(file)) {
        call_fixes <- call_fixes + 1
    }
    
    # Fix non-ASCII
    if (fix_non_ascii(file)) {
        ascii_fixes <- ascii_fixes + 1
    }
}

cat("\n=== Summary ===\n")
cat("Files processed:", length(r_files), "\n")
cat("Bracket notation fixes:", bracket_fixes, "files\n")
cat("Call quote fixes:", call_fixes, "files\n")
cat("Non-ASCII fixes:", ascii_fixes, "files\n")
cat("\nDone!\n")