# CRAN Submission Checklist

## Local Checks
- [ ] `R CMD check --as-cran` passes with no ERRORs, WARNINGs, or NOTEs
- [ ] Examples run in < 5 seconds each
- [ ] All functions have examples
- [ ] Package size < 5MB

## Remote Checks
- [ ] win-builder: `devtools::check_win_devel()`
- [ ] R-hub: `rhub::check_for_cran()`
- [ ] macOS: `rhub::check_on_macos()`

## Documentation
- [ ] NEWS.md updated
- [ ] DESCRIPTION version bumped
- [ ] All exported functions documented
- [ ] Vignettes build successfully

## Final Steps
- [ ] `devtools::spell_check()`
- [ ] `goodpractice::gp()`
- [ ] `urlchecker::url_check()`
