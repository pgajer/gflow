# Public content checks

Run `make check-public-content` before publishing source changes. Package build
targets also check the built source archive, including generated vignette HTML.
Run `make test-public-content` to exercise the detector's regression fixtures.
Python 3.9 or later and Git are required; no Python packages are needed.

For an archive built outside Make, run:

```sh
python3 tools/check-public-content.py --archive path/to/package_version.tar.gz
```

For generated website or documentation output, use `--tree path/to/output`.
Every invocation also checks the working checkout, including tracked, untracked
and ignored files. R build exclusions and Git ignore rules do not remove local
documents, old source snapshots or generated pages from this review. Git's own
metadata is not scanned. Use `--tree` for output outside the checkout or to
record a separate rendered-output result.

The checks identify personal filesystem references, task-coordination prose,
prompt fragments and misplaced working records. They inspect text, HTML comments
and entity/URL-escaped links. Ordinary scientific discussion of benchmarks,
audits, literature reviews and limitations is allowed. Reports distinguish
reference errors from wording that needs contextual review; both produce a
nonzero exit status until addressed.

Detailed JSON reports are written outside the checkout in the user's private
report directory, grouped by package. Set `R_PACKAGE_CONTENT_REPORT_DIR` or pass
`--report-dir` to choose another external directory. Reports inside the checkout
are rejected. CI does not upload these reports as public artifacts.

The `check`, `check-fast` and `check-examples` Make targets also write their R
check directories outside the checkout. Set `CHECK_DIR` to another external
directory when needed.

Fix the maintained source and regenerate its derived output. When a match is a
legitimate test fixture or documented tooling behavior, add a narrowly scoped
entry to `tools/public-content-allowlist.json`: the report's exact `scope`,
`path`, `rule` and `fingerprint`, plus a meaningful `reason`. Do not copy the
matched text into the exception file. Wildcard paths are rejected; changed or
removed exception content fails the next check and requires a fresh review.
Exceptions do not transfer between the repository, archive and rendered tree.

This is a targeted text check, not a replacement for scientific or editorial
review. Binary formats such as PDF and serialized R objects are counted but
not decoded; retain the package's format-specific checks and inspect published
figures and documents. A missing input, unreadable archive, external symbolic
link or incomplete scan fails the check rather than being treated as a pass.
