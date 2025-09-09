#!/bin/bash

# Script: check_migration_status.sh
# Purpose: Analyze migration status of R files from msr2 to gflow package
# Created: $(date '+%Y-%m-%d')
# Description: This script extracts filenames from msr2_dependencies.csv and compares
#              them with existing files in gflow/R to determine which files still need migration

# Set working directory
CLAUDE_DIR="$HOME/current_projects/gflow/.claude"
GFLOW_R_DIR="$HOME/current_projects/gflow/R"

echo "=== R File Migration Status Checker ==="
echo "Analyzing migration from msr2 to gflow..."
echo ""

# Step 1: Extract filename column from CSV (excluding header)
# tail -n +2: Skip the header row
# cut -d',' -f1: Extract first column (filename) using comma as delimiter
echo "Step 1: Extracting filenames from msr2_dependencies.csv..."
tail -n +2 "$CLAUDE_DIR/msr2_dependencies.csv" | cut -d',' -f1 > "$CLAUDE_DIR/msr2_filenames.txt"
MSR2_COUNT=$(wc -l < "$CLAUDE_DIR/msr2_filenames.txt")
echo "  Found $MSR2_COUNT files in msr2_dependencies.csv"

# Step 2: List all R files currently in gflow/R directory
# ls *.R: List all R files
# xargs -n1 basename: Extract just the filename (remove path)
# sort: Sort for later comparison
echo ""
echo "Step 2: Listing R files in gflow/R directory..."
ls "$GFLOW_R_DIR"/*.R 2>/dev/null | xargs -n1 basename | sort > "$CLAUDE_DIR/gflow_R_files.txt"
GFLOW_COUNT=$(wc -l < "$CLAUDE_DIR/gflow_R_files.txt")
echo "  Found $GFLOW_COUNT files in gflow/R"

# Step 3: Find files that are in msr2 but NOT in gflow (files to migrate)
# comm -23: Compare sorted files, show lines only in first file
# <(sort ...): Process substitution to sort files on the fly
echo ""
echo "Step 3: Identifying files that need migration..."
comm -23 <(sort "$CLAUDE_DIR/msr2_filenames.txt") <(sort "$CLAUDE_DIR/gflow_R_files.txt") > "$CLAUDE_DIR/msr2_files_to_migrate.txt"
TO_MIGRATE_COUNT=$(wc -l < "$CLAUDE_DIR/msr2_files_to_migrate.txt")
echo "  $TO_MIGRATE_COUNT files need to be migrated"

# Step 4: Find files that have been successfully migrated (in both lists)
# comm -12: Show lines that appear in both files
echo ""
echo "Step 4: Identifying already migrated files..."
comm -12 <(sort "$CLAUDE_DIR/msr2_filenames.txt") <(sort "$CLAUDE_DIR/gflow_R_files.txt") > "$CLAUDE_DIR/already_migrated_files.txt"
MIGRATED_COUNT=$(wc -l < "$CLAUDE_DIR/already_migrated_files.txt")
echo "  $MIGRATED_COUNT files from msr2 have been migrated"

# Step 5: Find files in gflow that are NOT from msr2 (native gflow files)
# comm -13: Show lines only in second file
echo ""
echo "Step 5: Identifying native gflow files (not from msr2)..."
comm -13 <(sort "$CLAUDE_DIR/msr2_filenames.txt") <(sort "$CLAUDE_DIR/gflow_R_files.txt") > "$CLAUDE_DIR/gflow_native_files.txt"
NATIVE_COUNT=$(wc -l < "$CLAUDE_DIR/gflow_native_files.txt")
echo "  $NATIVE_COUNT files are native to gflow"

# Calculate migration percentage
if [ $MSR2_COUNT -gt 0 ]; then
    MIGRATION_PERCENT=$((MIGRATED_COUNT * 100 / MSR2_COUNT))
else
    MIGRATION_PERCENT=0
fi

# Print summary report
echo ""
echo "============================================"
echo "           MIGRATION STATUS SUMMARY         "
echo "============================================"
echo "Total files in msr2_dependencies.csv: $MSR2_COUNT"
echo "Files successfully migrated:          $MIGRATED_COUNT"
echo "Files still to migrate:               $TO_MIGRATE_COUNT"
echo "Native gflow files:                   $NATIVE_COUNT"
echo "Total files in gflow/R:               $GFLOW_COUNT"
echo "Migration progress:                   ${MIGRATION_PERCENT}% complete"
echo ""
echo "Generated files:"
echo "  - $CLAUDE_DIR/msr2_filenames.txt (all msr2 files)"
echo "  - $CLAUDE_DIR/gflow_R_files.txt (all gflow/R files)"
echo "  - $CLAUDE_DIR/msr2_files_to_migrate.txt (files needing migration)"
echo "  - $CLAUDE_DIR/already_migrated_files.txt (successfully migrated)"
echo "  - $CLAUDE_DIR/gflow_native_files.txt (native gflow files)"

# Show sample of files to migrate
if [ $TO_MIGRATE_COUNT -gt 0 ]; then
    echo ""
    echo "Next files to migrate (showing first 10):"
    echo "-------------------------------------------"
    head -10 "$CLAUDE_DIR/msr2_files_to_migrate.txt" | sed 's/^/  - /'
fi

# Show recently migrated files
if [ $MIGRATED_COUNT -gt 0 ]; then
    echo ""
    echo "Recently migrated files (showing first 10):"
    echo "-------------------------------------------"
    head -10 "$CLAUDE_DIR/already_migrated_files.txt" | sed 's/^/  - /'
fi

echo ""
echo "Script completed successfully!"