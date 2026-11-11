# Checking size-type includes

The helper `tools/ensure_cstddef_size_t.sh` checks placement of `<cstddef>` and
`std::size_t` declarations. It requires Bash 4 or later (`mapfile` is used);
on macOS, use an installed modern Bash rather than `/bin/bash`. Run from the
package root, first in dry-run mode:

```sh
DRY_RUN=1 bash tools/ensure_cstddef_size_t.sh
```

Inspect the proposed changes before applying them:

```sh
bash tools/ensure_cstddef_size_t.sh
```
