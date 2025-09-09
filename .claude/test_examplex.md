# # Method 1: Test specific example manually
test_get_path_data_example <- function() {
  # Create a simple graph with 5 vertices
  adj.list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2), c(3))
  weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1), c(1))
  y <- c(1.5, 2.0, 0.5, 1.0, 1.5)
  
  # Find paths centered around vertex 2
  tryCatch({
    paths <- get.path.data(adj.list, weight.list, y,
                          ref.vertex = 2, bandwidth = 2)
    print("Example ran successfully!")
    return(paths)
  }, error = function(e) {
    print(paste("Example failed with error:", e$message))
    return(NULL)
  })
}

# Method 2: Use example() function after loading the package
test_with_example_function <- function(pkg_name) {
  library(pkg_name, character.only = TRUE)
  
  # Test specific function example
  tryCatch({
    example(get.path.data, package = pkg_name, give.lines = TRUE)
    print("Example ran successfully!")
  }, error = function(e) {
    print(paste("Example failed:", e$message))
  })
}

# Method 3: Run all examples in a package using tools::Rd2ex
test_all_examples <- function(pkg_dir) {
  library(tools)
  
  # Get all Rd files
  rd_files <- list.files(file.path(pkg_dir, "man"), 
                        pattern = "\\.Rd$", 
                        full.names = TRUE)
  
  results <- list()
  
  for (rd_file in rd_files) {
    func_name <- gsub("\\.Rd$", "", basename(rd_file))
    
    # Extract examples from Rd file
    tmp_file <- tempfile(fileext = ".R")
    tryCatch({
      tools::Rd2ex(rd_file, tmp_file)
      
      # Run the extracted examples
      result <- tryCatch({
        source(tmp_file, echo = FALSE)
        list(status = "success", func = func_name)
      }, error = function(e) {
        list(status = "error", func = func_name, error = e$message)
      })
      
      results[[func_name]] <- result
    }, error = function(e) {
      results[[func_name]] <- list(status = "no_examples", func = func_name)
    }, finally = {
      unlink(tmp_file)
    })
  }
  
  # Print summary
  success_count <- sum(sapply(results, function(x) x$status == "success"))
  error_count <- sum(sapply(results, function(x) x$status == "error"))
  no_examples <- sum(sapply(results, function(x) x$status == "no_examples"))
  
  cat("Example Test Summary:\n")
  cat(sprintf("  Success: %d\n", success_count))
  cat(sprintf("  Errors: %d\n", error_count))
  cat(sprintf("  No examples: %d\n", no_examples))
  
  # Show errors
  if (error_count > 0) {
    cat("\nErrors found in:\n")
    for (res in results) {
      if (res$status == "error") {
        cat(sprintf("  - %s: %s\n", res$func, res$error))
      }
    }
  }
  
  return(results)
}

# Method 4: Use devtools for comprehensive testing
test_with_devtools <- function(pkg_dir) {
  library(devtools)
  
  # Run all examples
  tryCatch({
    devtools::run_examples(pkg_dir)
    print("All examples passed!")
  }, error = function(e) {
    print(paste("Example check failed:", e$message))
  })
  
  # You can also use check() for full package check
  # devtools::check(pkg_dir, args = "--examples")
}

# Method 5: Quick test during development
quick_test_example <- function() {
  # Source the fixed function (if not in package yet)
  # source("path/to/get.path.data.R")
  
  # Run the problematic example
  adj.list <- list(c(2,3), c(1,3,4), c(1,2,5), c(2), c(3))
  weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1), c(1))
  y <- c(1.5, 2.0, 0.5, 1.0, 1.5)
  
  # Test edge case that caused the error (5 vertices, so 5/2 = 2.5)
  tryCatch({
    paths <- get.path.data(adj.list, weight.list, y,
                          ref.vertex = 2, bandwidth = 2,
                          diff.threshold = 2.4)  # Should work now
    print("Fixed! Function handles non-integer division correctly.")
  }, error = function(e) {
    print(paste("Still has error:", e$message))
  })
}

# Usage examples:
# test_get_path_data_example()
# test_with_example_function("gflow")
# test_all_examples("/path/to/gflow/package")
# test_with_devtools("/path/to/gflow/package")
# quick_test_example()

