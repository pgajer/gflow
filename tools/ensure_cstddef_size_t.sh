#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-inst/include/gflow}"

# List all headers that use size_t
mapfile -t FILES < <(grep -RIl --include='*.hpp' -e '\<size_t\>' "$ROOT")

for f in "${FILES[@]}"; do
  has_inc=0
  has_using=0

  grep -Eq '^[[:space:]]*#include[[:space:]]*<cstddef>'   "$f" && has_inc=1
  grep -Eq '^[[:space:]]*#include[[:space:]]*<stddef\.h>' "$f" && has_inc=1
  grep -Eq '^[[:space:]]*using[[:space:]]+std::size_t[[:space:]]*;' "$f" && has_using=1

  (( has_inc==1 && has_using==1 )) && { echo "OK   $f"; continue; }

  # Decide insertion point:
  #   1) after #pragma once, or
  #   2) after #ifndef ... / #define ... guard, or
  #   3) after top include block, else at top of file.
  insert_after=0

  guard_line=$(awk '
    /^#pragma[[:space:]]+once/            { print NR; exit }
    /^#ifndef[[:space:]]/                 { g=NR; next }
    g && /^#define[[:space:]]/            { print NR; exit }
  ' "$f")
  if [[ -n "${guard_line:-}" ]]; then
    insert_after="$guard_line"
  else
    top_includes=$(awk '
      BEGIN{last=0; started=0}
      NR<=200 {
        if ($0 ~ /^[[:space:]]*#include[[:space:]]*</) { last=NR; started=1 }
        else if (started && $0 !~ /^[[:space:]]*#include[[:space:]]*</) { print last; exit }
      }
      END{ if (last) print last }
    ' "$f")
    insert_after="${top_includes:-0}"
  fi

  lines=""
  (( has_inc==0   )) && lines+="#include <cstddef>\n"
  (( has_using==0 )) && lines+="using std::size_t;\n"

  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    echo "WOULD edit $f after line $insert_after:"
    printf "%b" "$lines"
  else
    tmp="$(mktemp)"
    if [[ "$insert_after" -eq 0 ]]; then
      printf "%b" "$lines" | cat - "$f" > "$tmp"
    else
      awk -v n="$insert_after" -v ins="$lines" 'NR==n{print; printf "%s", ins; next}1' "$f" > "$tmp"
    fi
    mv "$tmp" "$f"
    echo "EDIT $f"
  fi
done
