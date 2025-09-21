.PHONY: clean build build-verbose check check-fast install document attrs
VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflow
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
LOGDIR := .claude
export PATH := /opt/homebrew/bin:$(PATH)

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
	@R -q -e "devtools::document()" > $(LOGDIR)/$(PKGNAME)_document.log 2>&1
	@echo "Documentation generated (log: $(LOGDIR)/$(PKGNAME)_document.log)"

build: clean document
	@mkdir -p $(LOGDIR)
	@echo "Building package..."
	@cd .. && R CMD build $(PKGNAME) > $(PKGNAME)/$(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Package built successfully (log: $(LOGDIR)/$(PKGNAME)_build.log)"

build-verbose: clean document
	cd .. && R CMD build $(PKGNAME)

build-log: clean document
	@mkdir -p $(LOGDIR)
	cd .. && R CMD build $(PKGNAME) > $(PKGNAME)/$(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Build output saved to $(LOGDIR)/$(PKGNAME)_build.log"

check: build
	cd .. && R CMD check $(TARBALL) --as-cran

check-fast: build
	cd .. && R CMD check $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check-examples: build
	cd .. && R CMD check $(TARBALL) --as-cran --examples

install: build
	cd .. && R CMD INSTALL $(TARBALL)

rchk:
	@tools/check_rchk.sh
