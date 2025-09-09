# Single File Migration from msr2 to gflow

## Task Overview
Migrate ONE file from msr2/R to gflow/R 

## Migration Steps

### 1. Select and Copy specified File
```bash
# Copy the file from msr2/R to gflow/R
cp ~/current_projects/msr2/R/[filename] ~/current_projects/gflow/R/
```

### 2. Fix Documentation Issues

#### 2.1 Non-ASCII Characters
Search for non-ASCII characters (especially Greek symbols) in documentation and replace them with ASCII codes:
- α → \alpha
- β → \beta  
- γ → \gamma
- λ → \lambda
- μ → \mu
- σ → \sigma
- Σ → \Sigma
- ∞ → \infty
- ≤ → \le
- ≥ → \ge
- ≠ → \neq

Check all roxygen2 documentation blocks (@param, @return, @description, @details, @examples).

#### 2.2 Square Brackets Interpreted as Links
Find patterns like `[i]`, `[i,j]`, `[,j]` in documentation that should represent indices, not links.
Replace with:
- `[i]` → `\[i\]`
- `[i,j]` → `\[i,j\]`
- `[,j]` → `\[,j\]`
- `[[i]]` → `\[\[i\]\]`

Common locations: @param descriptions for matrix/vector arguments, @return descriptions, @examples.

#### 2.3 Itemize Environments
Convert any `\itemize{}` blocks in documentation to proper roxygen2 format:

**From:**
```r
#' @param x \itemize{
#'   \item First item
#'   \item Second item
#' }
```

**To:**
```r
#' @param x A parameter with options:
#' \itemize{
#'   \item First item
#'   \item Second item
#' }
```
Ensure itemize blocks are properly formatted with correct indentation.

#### R Package C Function Call Syntax
When reviewing R package code, check that all calls to C functions via `.C()`
use string literals for the C function name. A common error is writing
`.C(C_function_name, ...)` instead of the correct `.C("C_function_name", ...)`.
The C function name must be passed as a string to the `.C()` function, not as a
bare symbol. This error will cause runtime failures with messages like "object
'C_function_name' not found" during package checks or when running examples.
When you encounter this pattern, immediately fix it by adding quotes around the
C function name. This same principle applies to other R-to-C interface functions
like `.Call()` and `.Fortran()`.


### 3. Verify Documentation Consistency

#### 3.1 Function-Documentation Alignment
For each function in the file:
1. Check that all function parameters have corresponding @param entries
2. Verify parameter names match exactly
3. Ensure @return documentation exists if function returns a value
4. Check that @export is present if function should be exported

#### 3.2 Common Documentation Issues
- Missing @param entries for function arguments
- @param entries for removed arguments
- Mismatched argument names between function and documentation
- Missing @examples or broken example code
- Incorrect @importFrom statements

### 4. Update Import Directives
Check if the file uses functions from other packages and ensure proper roxygen2 directives:

#### 4.1 Identify External Function Calls
Look for:
- `package::function()` calls in the code
- Direct function calls that might be from external packages

#### 4.2 Add Missing @importFrom Directives
For each external function used, ensure the roxygen2 block has:
```r
#' @importFrom package function1 function2
```

#### 4.3 Check for Missing Package Dependencies
If using functions from packages not in gflow/DESCRIPTION, add them to Imports:
- Edit gflow/DESCRIPTION
- Add missing packages to the Imports: field

#### 4.4 Regenerate NAMESPACE
After adding all necessary @importFrom directives:
```bash
cd ~/current_projects/gflow
Rscript -e "devtools::document()"
```

**IMPORTANT**: Never manually edit the NAMESPACE file. All imports/exports must be managed through roxygen2 directives in the R files.

### 5. Run Initial Check
```bash
cd ~/current_projects/gflow
R CMD build . && R CMD check gflow_*.tar.gz --as-cran
```

### 6. Fix Issues Iteratively
If check reveals issues:

#### 6.1 Documentation Issues
- Fix any "undocumented arguments" warnings
- Resolve "missing documentation" errors
- Correct malformed documentation syntax

#### 6.2 Code Issues  
- Fix any undefined global variables
- Resolve namespace conflicts
- Correct S3 method consistency

#### 6.3 Example Issues
- Ensure examples run without errors
- Add \dontrun{} for long-running examples
- Fix any missing data or function calls in examples

### 7. Iterate Until Clean
Repeat the cycle:
1. Fix identified issues
2. Run `R CMD check gflow_*.tar.gz --as-cran`
3. Continue until no warnings or errors remain

### 8. Final Verification
Once check passes:
1. Document the migrated file in a log:
   ```
   echo "[DATE] [FILENAME] migrated successfully" >> ~/current_projects/gflow/.claude/migration_log.txt
   ```
2. Update the CSV to mark this file as processed
3. Report completion with summary of changes made

## Output Format
Please provide:
1. The name of the file being migrated
2. List of documentation fixes applied (greek symbols, brackets, itemize, etc.)
3. Any function/documentation mismatches found and fixed
4. Initial R CMD check issues encountered
5. How those issues were resolved
6. Final R CMD check status (should be clean)

## Important Notes
- Only migrate ONE file per run
- Ensure R CMD check passes completely before considering the migration complete
- If circular dependencies are discovered, stop and report the issue
- Preserve all existing functionality - do not modify function behavior
