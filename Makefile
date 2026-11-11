.PHONY: clean build build-verbose check check-fast install document manifest attrs check-r-toolchain audit-malo-exports audit-s3-namespace audit-phase7-ownership audit-cleanup-boundary audit-final-acceptance audit-api-guide preview-user-guides preview-knn-basins install-user showcase website
VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflow
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
CHECK_DIR ?= $(shell python3 -c "from pathlib import Path; print(Path.home() / '.codex' / 'private' / 'gflow' / 'checks')")
HOMEBREW_BIN := /opt/homebrew/bin
GCC_BIN := /opt/homebrew/opt/gcc/bin
TIDY_BIN := $(shell if [ -x "$(HOMEBREW_BIN)/tidy" ]; then echo "$(HOMEBREW_BIN)/tidy"; elif command -v tidy >/dev/null 2>&1; then command -v tidy; else echo "$(HOMEBREW_BIN)/tidy"; fi)
R_ENV := env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE
R_RUN := $(R_ENV) R
RSCRIPT_RUN := $(R_ENV) Rscript

clean:
	rm -f src/Makevars.local src/Makevars.win.local
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
	$(RSCRIPT_RUN) tools/update_guide_links.R

manifest:
	$(RSCRIPT_RUN) tools/update_gflow_build_manifest.R

.PHONY: check-public-content test-public-content
check-public-content:
	@python3 tools/check-public-content.py

test-public-content:
	@python3 tools/tests/test-public-content.py

build: check-public-content clean document manifest
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	$(R_RUN) CMD build .
	python3 tools/check-public-content.py --archive "$(TARBALL)"
	python3 tools/check-package-size.py "$(TARBALL)"

build-verbose: check-public-content clean document
	$(R_RUN) CMD build .
	python3 tools/check-public-content.py --archive "$(TARBALL)"
	python3 tools/check-package-size.py "$(TARBALL)"

build-log: check-public-content clean document
	$(R_RUN) CMD build .
	python3 tools/check-public-content.py --archive "$(TARBALL)"
	python3 tools/check-package-size.py "$(TARBALL)"

check: build
	mkdir -p "$(CHECK_DIR)"
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" $(R_RUN) CMD check --output="$(CHECK_DIR)" $(TARBALL) --as-cran

check-fast: build
	mkdir -p "$(CHECK_DIR)"
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" $(R_RUN) CMD check --output="$(CHECK_DIR)" $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check-examples: build
	mkdir -p "$(CHECK_DIR)"
	PATH="$(GCC_BIN):$(HOMEBREW_BIN):$$PATH" R_TIDYCMD="$(TIDY_BIN)" $(R_RUN) CMD check --output="$(CHECK_DIR)" $(TARBALL) --as-cran --examples

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

audit-final-acceptance: audit-cleanup-boundary
	@$(RSCRIPT_RUN) tools/check_final_acceptance.R

audit-api-guide:
	@$(RSCRIPT_RUN) tools/audit_api_guide.R

preview-user-guides:
	@$(RSCRIPT_RUN) tools/render_user_guides.R

preview-knn-basins:
	@$(RSCRIPT_RUN) tools/render_knn_basins_vignette.R

# User installation keeps the same libraries as the dependency bootstrap.
install-user:
	Rscript --vanilla tools/install_build_dependencies.R
	$(MAKE) install R_ENV="env -u R_HOME"
	Rscript --vanilla tools/check_installed_guides.R

showcase:
	$(RSCRIPT_RUN) tools/render_basin_showcase.R

website: check-public-content
	$(RSCRIPT_RUN) tools/build_website.R
	python3 tools/check-public-content.py --tree build/site
