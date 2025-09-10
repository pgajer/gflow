# Global variables used in the package
# This file declares global variables to avoid R CMD check NOTEs

utils::globalVariables(c(
  # Variables used in various functions
  "x",
  "y", 
  "n",
  "L",
  "R",
  
  # Function names that are dynamically called or missing
  "IW.kNN.graph",
  "rm.duplicate.grid.elements"
  ))
