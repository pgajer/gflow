#!/bin/bash
# CRAN check script for gflow package

# Exit on error
set -e

echo "=== Starting CRAN check process for gflow ==="

# Clean previous builds
echo "1. Cleaning previous builds..."
cd ~/current_projects/gflow
find src -name "*.o" -delete
find src -name "*.so" -delete
rm -f src/*.dll  # For Windows compatibility

# Build the package
echo "2. Building package..."
cd ~/current_projects
R CMD build gflow --no-build-vignettes

# Get the version dynamically
VERSION=$(grep "^Version:" gflow/DESCRIPTION | sed 's/Version: //')
TARBALL="gflow_${VERSION}.tar.gz"

# Full CRAN check
echo "3. Running CRAN checks..."
R CMD check $TARBALL --as-cran --no-manual

# Optional: Run additional checks
echo "4. Running additional platform checks..."
# Uncomment these as needed:
# R -e "rhub::check_for_cran('$TARBALL')"
# R -e "devtools::check_win_devel('gflow')"

echo "=== CRAN check complete ==="
echo "Check results in: gflow.Rcheck/"
