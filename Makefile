.PHONY: clean build build-verbose check check-fast install document attrs audit-malo-exports
VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflow
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
LOGDIR := .claude
HOMEBREW_BIN := /opt/homebrew/bin
GCC_BIN := /opt/homebrew/opt/gcc/bin

clean:
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	rm -f src/RcppExports.cpp
	rm -rf $(PKGNAME).Rcheck
	rm -f $(TARBALL)
	rm -f $(LOGDIR)/*.log

# 1) Always (re)generate Rcpp glue first
attrs:
	@mkdir -p $(LOGDIR)
	@echo "Running Rcpp::compileAttributes()..."
	@R -q -e "Rcpp::compileAttributes()" > $(LOGDIR)/$(PKGNAME)_rcppattrs.log 2>&1
	@echo "RcppExports regenerated (log: $(LOGDIR)/$(PKGNAME)_rcppattrs.log)"

# 2) Then regenerate NAMESPACE + Rd via roxygen (through devtools::document)
document: attrs
	@mkdir -p $(LOGDIR)
	@echo "Running devtools::document()..."
	@PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R -q -e "roxygen2::roxygenise(load = 'source')" > $(LOGDIR)/$(PKGNAME)_document.log 2>&1
	@echo "Documentation generated (log: $(LOGDIR)/$(PKGNAME)_document.log)"

build: clean document
	@mkdir -p $(LOGDIR)
	@echo "Building package..."
	@R CMD build . --output=. > $(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Package built successfully (log: $(LOGDIR)/$(PKGNAME)_build.log)"

build-verbose: clean document
	R CMD build . --output=.

build-log: clean document
	@mkdir -p $(LOGDIR)
	R CMD build . --output=. > $(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Build output saved to $(LOGDIR)/$(PKGNAME)_build.log"

check: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran

check-fast: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check-examples: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(HOMEBREW_BIN)/tidy" R CMD check $(TARBALL) --as-cran --examples

install: build
	R CMD INSTALL $(TARBALL)

rchk:
	@tools/check_rchk.sh

audit-malo-exports:
	@Rscript tools/audit_malo_exports.R
