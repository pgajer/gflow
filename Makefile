.PHONY: clean build build-verbose check check-fast install document

VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflow
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
LOGDIR := .claude

clean:
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	rm -rf $(PKGNAME).Rcheck
	rm -f $(TARBALL)
	rm -f $(LOGDIR)/*.log

document:
	@mkdir -p $(LOGDIR)
	@echo "Running roxygen2..."
	@R -e "roxygen2::roxygenise()" > $(LOGDIR)/$(PKGNAME)_roxygen2.log 2>&1
	@echo "Documentation generated (log: $(LOGDIR)/$(PKGNAME)_roxygen2.log)"

build: clean document
	@mkdir -p $(LOGDIR)
	@echo "Building package..."
	@cd .. && R CMD build $(PKGNAME) > $(PKGNAME)/$(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Package built successfully (log: $(LOGDIR)/$(PKGNAME)_build.log)"

build-verbose: clean document
	cd .. && R CMD build $(PKGNAME)

# Build with output logged to file
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
