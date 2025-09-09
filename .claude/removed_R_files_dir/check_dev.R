# Quick development check script
library(devtools)

# Document and update NAMESPACE
devtools::document()

# Run tests
devtools::test()

# Quick check (faster than R CMD check)
devtools::check(args = "--no-examples --no-tests")

# Full check when ready
# devtools::check(args = "--as-cran")
