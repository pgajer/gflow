# pkg building/installation instruction

## After adding a C++ file run:

Rcpp::compileAttributes()  # generates R/RcppExports.R + src/RcppExports.cpp
devtools::document()       # regenerates NAMESPACE + Rd from roxygen


## devtools tools 

devtools::load_all()
devtools::clean_dll(); devtools::load_all()
devtools::install()
devtools::document()



