#!/usr/bin/env /opt/homebrew/bin/bash

set -euo pipefail

ROOT="${1:-inst/include}"

# Process all headers that contain 'using std::size_t;'
grep -RIl --include='*.hpp' -e 'using[[:space:]]\+std::size_t[[:space:]]*;' "$ROOT" \
| while IFS= read -r f; do
  tmp="$(mktemp)"

  awk '
    BEGIN { seen_inc=0; inserted=0 }
    # If we hit include<cstddef>, note it. If we already inserted one earlier, skip duplicates.
    /^[[:space:]]*#include[[:space:]]*<cstddef>[[:space:]]*$/ {
      if (inserted) { next } else { seen_inc=1; print; next }
    }
    # When we hit the using line but haven’t seen an include yet, insert it right above.
    /^[[:space:]]*using[[:space:]]+std::size_t[[:space:]]*;[[:space:]]*$/ {
      if (!seen_inc) { print "#include <cstddef>"; inserted=1; seen_inc=1 }
      print; next
    }
    { print }
  ' "$f" > "$tmp"

  # Only overwrite if changed
  if ! cmp -s "$f" "$tmp"; then
    mv "$tmp" "$f"
    echo "FIXED order: $f"
  else
    rm -f "$tmp"
  fi
done
