.PHONY: clean build build-verbose check check-fast install document attrs check-r-toolchain audit-malo-exports audit-s3-namespace audit-phase7-ownership audit-cleanup-boundary
VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflow
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
HOMEBREW_BIN := /opt/homebrew/bin
GCC_BIN := /opt/homebrew/opt/gcc/bin
TIDY_BIN := $(shell if [ -x "$(HOMEBREW_BIN)/tidy" ]; then echo "$(HOMEBREW_BIN)/tidy"; elif command -v tidy >/dev/null 2>&1; then command -v tidy; else echo "$(HOMEBREW_BIN)/tidy"; fi)
R_ENV := env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE
R_RUN := $(R_ENV) R
RSCRIPT_RUN := $(R_ENV) Rscript

clean:
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	rm -f src/RcppExports.cpp
	rm -rf $(PKGNAME).Rcheck
	rm -f $(TARBALL)
	rm -rf .claude

check-r-toolchain:
	$(R_RUN) --vanilla --slave -f tools/check_r_toolchain.R

# 1) Always (re)generate Rcpp glue first
attrs: check-r-toolchain
	$(R_RUN) -q -e "Rcpp::compileAttributes()"

# 2) Then regenerate NAMESPACE + Rd via roxygen (through devtools::document)
document: attrs
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" $(R_RUN) -q -e "roxygen2::roxygenise(load = 'source')"

build: clean document
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	$(R_RUN) CMD build .

build-verbose: clean document
	$(R_RUN) CMD build .

build-log: clean document
	$(R_RUN) CMD build .

check: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" $(R_RUN) CMD check $(TARBALL) --as-cran

check-fast: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" $(R_RUN) CMD check $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check-examples: build
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" $(R_RUN) CMD check $(TARBALL) --as-cran --examples

install: build
	$(R_RUN) CMD INSTALL $(TARBALL)

rchk:
	@tools/check_rchk.sh

audit-malo-exports:
	@$(RSCRIPT_RUN) tools/audit_malo_exports.R

audit-s3-namespace:
	@$(RSCRIPT_RUN) tools/check_s3_namespace.R

audit-phase7-ownership:
	@$(RSCRIPT_RUN) tools/check_phase7_ownership.R

audit-cleanup-boundary: audit-s3-namespace audit-phase7-ownership
	@$(RSCRIPT_RUN) tools/check_cleanup_guardrails.R
