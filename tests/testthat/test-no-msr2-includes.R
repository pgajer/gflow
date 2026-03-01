test_that("src does not include legacy msr2.h", {
  src_dir <- testthat::test_path("..", "..", "src")
  src_files <- list.files(src_dir, pattern = "\\.(c|cc|cpp|cxx|h|hpp)$", full.names = TRUE)

  offenders <- character()
  for (f in src_files) {
    lines <- readLines(f, warn = FALSE)
    if (any(grepl('^\\s*#\\s*include\\s*"msr2\\.h"', lines))) {
      offenders <- c(offenders, f)
    }
  }

  if (length(offenders) > 0) {
    fail(paste("Offending files:", paste(offenders, collapse = ", ")))
  }
})
