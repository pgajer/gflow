# R Package Migration Instructions: msr2 to gflow

## Overview
You are tasked with systematically migrating R files from the `msr2` package to the `gflow` package. The `gflow` package currently contains 16 R files from `msr2`, which has approximately 100 R files total. Your goal is to intelligently migrate the remaining ~84 files in stages.

## Initial Analysis Tasks

### 1. Dependency Mapping
First, analyze the dependency structure of all R files in msr2:

```bash
# Navigate to msr2 package directory
cd /path/to/msr2

# Create a dependency analysis report
```

Generate a CSV file (`msr2_dependencies.csv`) with the following columns:
- `filename`: Name of the R file
- `imports`: Functions/files this file imports or sources
- `exports`: Functions this file exports
- `dependencies`: Other R files in msr2 this file depends on
- `dependents`: Other R files that depend on this file
- `external_deps`: External package dependencies
- `complexity`: Low/Medium/High based on LOC and dependencies

### 2. Current State Assessment
Analyze the 16 files already in gflow:
- List all files currently in `gflow/R/`
- Identify which functions from msr2 they might be calling
- Note any missing dependencies

### 3. File Categorization
Create a migration plan (`migration_plan.csv`) categorizing the remaining files:

| Category | Priority | Description | Example Files |
|----------|----------|-------------|---------------|
| core_utils | 1 | Utility functions with no dependencies | utils.R, constants.R |
| data_structures | 1 | Core classes and data types | classes.R, validators.R |
| standalone_features | 2 | Self-contained feature modules | plot_functions.R, report_gen.R |
| integration_layer | 3 | Files connecting multiple components | workflow.R, pipeline.R |
| high_dependency | 4 | Files with complex interdependencies | main_analysis.R |
| optional | 5 | Legacy or rarely used functions | deprecated.R |

## Migration Stages

### Stage 1: Foundation (Priority 1 files)
1. Identify all files with zero dependencies on other msr2 R files
2. Among these, prioritize files that:
   - Define S3/S4 classes or methods
   - Contain utility functions
   - Define package-wide constants
3. For each file:
   ```bash
   # Copy file to gflow
   cp msr2/R/[filename] gflow/R/
   
   # Update NAMESPACE if needed
   # Update documentation
   # Run tests
   cd gflow && Rscript -e "devtools::test()"
   ```

### Stage 2: Feature Modules (Priority 2 files)
1. Identify standalone feature sets (groups of 2-5 related files)
2. Verify they only depend on Stage 1 files or external packages
3. Migrate each feature set as a unit

### Stage 3: Integration (Priority 3-4 files)
1. Sort remaining files by number of dependencies (ascending)
2. Migrate in dependency order
3. Handle circular dependencies by:
   - Identifying the cycle
   - Determining if refactoring is needed
   - Moving files together if they're tightly coupled

## Specific Instructions

### For Each File Migration:

1. **Pre-migration checks:**
   ```r
   # Check if file depends on unmigrated files
   # List all functions called from other msr2 files
   # Verify external package dependencies are in gflow DESCRIPTION
   ```

2. **Migration steps:**
   ```bash
   # Copy file
   cp msr2/R/[filename] gflow/R/
   
   # If file has corresponding tests
   cp msr2/tests/testthat/test-[filename] gflow/tests/testthat/
   
   # If file has man pages
   cp msr2/man/[related_functions].Rd gflow/man/
   ```

3. **Post-migration validation:**
   ```r
   # In gflow directory
   devtools::load_all()
   devtools::check()
   devtools::test()
   ```

4. **Update tracking:**
   - Add entry to `migration_log.csv` with: date, filename, status, issues
   - Update `migration_plan.csv` status column

### Handling Common Issues:

1. **Namespace conflicts:**
   - Check for function name collisions
   - Use `gflow::` prefix if needed
   - Update NAMESPACE file

2. **Missing dependencies:**
   - Add to DESCRIPTION file
   - Document minimum version requirements

3. **Circular dependencies:**
   - Document the cycle in `circular_deps.md`
   - Propose refactoring approach
   - Get approval before breaking cycles

## Output Deliverables

Create these files in a `migration_analysis` directory:

1. **msr2_dependencies.csv** - Complete dependency analysis
2. **migration_plan.csv** - Staged migration plan with priorities
3. **migration_log.csv** - Track progress and issues
4. **circular_deps.md** - Document any circular dependencies
5. **missing_deps.md** - List any missing external dependencies
6. **refactoring_notes.md** - Suggested improvements during migration

## Automation Script

Create an R script `migrate_files.R` that:
```r
# Function to migrate a single file with all checks
migrate_file <- function(filename, source_pkg = "msr2", target_pkg = "gflow") {
  # Pre-checks
  # Copy file
  # Update NAMESPACE
  # Run tests
  # Log results
}

# Function to migrate a stage
migrate_stage <- function(stage_number) {
  # Read migration_plan.csv
  # Filter for stage
  # Call migrate_file for each
}
```

## Final Notes

- Commit after each successful file migration
- Create a branch for each stage: `migration-stage-1`, etc.
- If a file breaks tests, revert and document the issue
- Prioritize maintaining backward compatibility
- Consider creating adapter functions for breaking changes

Begin by running the dependency analysis and creating the initial CSV files. Then proceed with Stage 1 migration.
