.PHONY: clean build build-verbose check check-fast check-correctness-report install document attrs audit-malo-exports
VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflow
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
HOMEBREW_BIN := /opt/homebrew/bin
GCC_BIN := /opt/homebrew/opt/gcc/bin
TIDY_BIN := $(shell if [ -x "$(HOMEBREW_BIN)/tidy" ]; then echo "$(HOMEBREW_BIN)/tidy"; elif command -v tidy >/dev/null 2>&1; then command -v tidy; else echo "$(HOMEBREW_BIN)/tidy"; fi)

clean:
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	rm -f src/RcppExports.cpp
	rm -rf $(PKGNAME).Rcheck
	rm -f $(TARBALL)
	rm -rf .claude

# 1) Always (re)generate Rcpp glue first
attrs:
	R -q -e "Rcpp::compileAttributes()"

# 2) Then regenerate NAMESPACE + Rd via roxygen (through devtools::document)
document: attrs
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R -q -e "roxygen2::roxygenise(load = 'source')"

build: clean document
	R CMD build .

build-verbose: clean document
	R CMD build .

build-log: clean document
	R CMD build .

check: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" R CMD check $(TARBALL) --as-cran

check-fast: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" R CMD check $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check-examples: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" R CMD check $(TARBALL) --as-cran --examples

check-correctness-report:
	Rscript tests/correctness/rdgraph_regression_report.R

install: build
	R CMD INSTALL $(TARBALL)

rchk:
	@tools/check_rchk.sh

audit-malo-exports:
	@Rscript tools/audit_malo_exports.R
