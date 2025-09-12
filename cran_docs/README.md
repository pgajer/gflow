# `cran_docs/` – CRAN Submission, Maintenance & Dev Playbooks

This folder centralizes everything *not shipped* with the package build but crucial for
releasing to CRAN, keeping the package healthy, and preserving institutional memory.

> Add `^cran_docs/` to `.Rbuildignore` to keep these files out of the build tarball.

## Structure

```
cran_docs/
├─ README.md                  # this file
├─ submission/                # material you prepare to talk to CRAN (or attach in emails)
│  ├─ cran_comments.md        # editable template for the "cran-comments.md" CRAN expects
│  ├─ rhub_results.md         # summarized outputs from rhub checks (incl. rchk)
│  ├─ winbuilder_results.md   # win-builder runs and notes
│  └─ reverse_depends.md      # notes from revdep checks
├─ maintenance/               # recurring hygiene & incident response logs
│  ├─ rchk_issues.md          # protection/GC analysis and remediation
│  ├─ sanitizer_reports.md    # ASAN/UBSAN/valgrind summaries
│  ├─ check_results.md        # `R CMD check` across platforms
│  └─ portability.md          # OS/toolchain quirks and workarounds
└─ development/               # internal engineering guides
   ├─ protect_patterns.md     # PROTECT/UNPROTECT coding standards
   ├─ rcpp_migration.md       # plan & conventions for migrating wrappers to Rcpp
   ├─ testing.md              # how to run tests locally/CI; seeds; reproducibility
   └─ release_checklist.md    # end-to-end steps for a clean release
```

## Workflow Tips

- Keep each CRAN interaction reproducible: paste exact messages and your response in `submission/cran_comments.md`.
- After every major CI run (rhub/win-builder), paste *summaries* here, link to raw logs.
- When you fix an issue class (e.g., a PROTECT pattern), document it once in `development/`, link it from incident notes in `maintenance/`.
- Treat this folder as your long-lived memory: it saves hours every release.
