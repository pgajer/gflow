# PROTECT/UNPROTECT Fixes for SEXP Conversion Functions

## Summary
This document describes the refactoring of conversion functions in `SEXP_cpp_conversion_utils.cpp` to fix PROTECT/UNPROTECT issues identified by `rhub::rhub_check(platforms = c("rchk"))`.

## Issues Identified by rchk

### 1. Stack-depth accounting errors (CRITICAL - CRAN auto-fail)
- **Issue**: "negative depth", "attempt to unprotect more items than protected"
- **Cause**: Mismatched PROTECT/UNPROTECT counts, conditional returns, early errors
- **Affected functions**: Multiple conversion functions with inconsistent protection patterns

### 2. Unsupported UNPROTECT forms (CRITICAL - CRAN auto-fail) 
- **Issue**: "unsupported form of unprotect … results will be incomplete"
- **Cause**: Using variables or expressions instead of literal constants in UNPROTECT
- **Example**: `UNPROTECT(n+m)` or `UNPROTECT(some_var)` where some_var changes

### 3. Protection stack imbalance (HIGH PRIORITY)
- **Issue**: Can't prove correctness of protection patterns
- **Affected**: Vector conversion helpers, nested structure builders

### 4. Unprotected live objects across allocating calls (HIGH PRIORITY)
- **Issue**: "unprotected variable X while calling allocating function"
- **Cause**: Creating SEXP then calling allocating function before protecting
- **Example**: 
  ```cpp
  SEXP names = Rf_allocVector(STRSXP, m);  // unprotected
  SET_STRING_ELT(names, 0, Rf_mkChar("a")); // mkChar allocates -> GC risk
  ```

## Refactoring Principles Applied

### A. One Disciplined Counter Per Function
```cpp
int nprot = 0;
SEXP ans = PROTECT(Rf_allocVector(VECSXP, k)); ++nprot;
// ... whenever you PROTECT(x), ++nprot
UNPROTECT(nprot);  // single exit point
return ans;
```

### B. Protect Immediately
**Bad:**
```cpp
SEXP names = Rf_allocVector(STRSXP, m);  // unprotected
SET_STRING_ELT(names, 0, Rf_mkChar("a")); // GC risk
```

**Good:**
```cpp
SEXP names = PROTECT(Rf_allocVector(STRSXP, m)); ++nprot;
SET_STRING_ELT(names, 0, Rf_mkChar("a")); // safe
```

### C. Container-First Pattern for Lists
1. Allocate and protect container first
2. Allocate and protect element
3. Fill element
4. Set element in container (container now holds reference)
5. Unprotect element (safe because container holds it)

### D. Loop Protection Pattern
```cpp
SEXP list = PROTECT(Rf_allocVector(VECSXP, n)); ++nprot;
for (int i = 0; i < n; ++i) {
    SEXP elem = PROTECT(Rf_allocVector(REALSXP, m)); ++nprot;
    // fill elem
    SET_VECTOR_ELT(list, i, elem);
    UNPROTECT(1); --nprot;  // elem now held by list
}
```

## Functions Refactored

### 1. `convert_vector_double_to_R`
- **Issue**: Returned protected SEXP without clear contract
- **Fix**: Single nprot counter, unprotect before return

### 2. `convert_vector_int_to_R`
- **Issue**: Same as above
- **Fix**: Same pattern applied

### 3. `convert_vector_bool_to_R`
- **Issue**: Same as above, plus std::vector<bool> special handling
- **Fix**: Explicit TRUE/FALSE assignment

### 4. `convert_vector_vector_double_to_R`
- **Issue**: PROTECT/UNPROTECT imbalance in loop
- **Fix**: Proper loop pattern with per-iteration unprotect

### 5. `convert_vector_vector_int_to_R`
- **Issue**: Same as above
- **Fix**: Same pattern

### 6. `convert_vector_vector_bool_to_R`
- **Issue**: Same as above
- **Fix**: Same pattern with bool handling

### 7. `convert_vector_vector_double_to_matrix`
- **Issue**: Single PROTECT without matching UNPROTECT
- **Fix**: Proper nprot counter and UNPROTECT

### 8. `convert_map_int_vector_int_to_R`
- **Issue**: Names created unprotected, then mkChar called
- **Fix**: Protect names before mkChar, proper ordering

### 9. `convert_wgraph_to_R`
- **Issue**: Complex protection with variable nprot increments
- **Fix**: Disciplined counter with proper increment/decrement

## Integration Instructions

### Step 1: Backup Original
```bash
cp src/SEXP_cpp_conversion_utils.cpp src/SEXP_cpp_conversion_utils.cpp.backup
```

### Step 2: Replace Functions
Replace each function in `src/SEXP_cpp_conversion_utils.cpp` with its corresponding version from `cran_docs/maintenance/fixed_converts.cpp`.

### Step 3: Test Compilation
```bash
R CMD build .
R CMD INSTALL gflow_*.tar.gz
```

### Step 4: Verify with rchk
```r
rhub::rhub_check(platforms = c("rchk"))
```

### Step 5: Run Package Tests
```r
devtools::test()
devtools::check()
```

## Key Changes Summary

| Function | Before | After |
|----------|--------|-------|
| All functions | Mixed PROTECT patterns | Single nprot counter |
| All functions | UNPROTECT at various points | UNPROTECT(nprot) at end |
| Loop functions | PROTECT without UNPROTECT in loop | PROTECT/UNPROTECT per iteration |
| Functions with names | Names unprotected during mkChar | Names protected before mkChar |
| Matrix function | Missing UNPROTECT | Proper UNPROTECT |
| Map function | Complex protection pattern | Ordered: container, names, elements |

## Testing Checklist

- [ ] All functions compile without warnings
- [ ] Package builds successfully
- [ ] Package installs successfully
- [ ] All existing tests pass
- [ ] rchk shows no PROTECT errors for these functions
- [ ] No memory leaks in valgrind
- [ ] CRAN check passes

## Notes

1. **Return Value Protection**: The refactored functions return UNPROTECTED SEXPs. Callers must PROTECT if needed.

2. **Empty Input Handling**: All functions properly handle empty inputs with balanced PROTECT/UNPROTECT.

3. **Error Paths**: All error paths (via Rf_error) occur before any PROTECTs or after proper cleanup.

4. **Performance**: The refactoring should have negligible performance impact as it only changes protection patterns, not algorithms.

5. **Compatibility**: The function signatures and behavior remain unchanged, ensuring backward compatibility.

## Remaining Work

After fixing these conversion functions, similar patterns should be applied to:
- Entry point functions mentioned in rchk report (S_agemalo, S_ray_agemalo, etc.)
- Other helper functions with protection issues
- Any new functions added to the codebase

## References

- [Writing R Extensions - Section 5.9.1](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Garbage-Collection)
- [rchk documentation](https://github.com/kalibera/rchk)
- [Tomas Kalibera's blog on PROTECT errors](https://blog.r-project.org/2019/04/18/common-protect-errors/)