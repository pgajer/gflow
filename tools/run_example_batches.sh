#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: tools/run_example_batches.sh <manifest>

Manifest format:
  One batch per line, using tab-separated fields:

    <commit message><TAB><path1 path2 ...>

  Blank lines and lines starting with '#' are ignored.

Behavior:
  - For each batch, the script checks whether any listed paths have staged or
    unstaged changes.
  - If no listed paths changed, the batch is skipped.
  - If changes are present, the script runs:
      1. make document
      2. Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
      3. make check
  - On success, it stages only the listed paths, commits with the provided
    message, and pushes to origin/main.

Notes:
  - This script automates validation/checkpointing. It does not generate the
    source edits for a batch.
  - Stop conditions are strict: the script aborts on the first failing batch.
EOF
}

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

manifest=$1

if [[ ! -f "$manifest" ]]; then
  echo "Manifest not found: $manifest" >&2
  exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "This script must be run inside a git repository." >&2
  exit 1
fi

run_qa() {
  make document
  Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
  make check
}

batch_num=0

while IFS=$'\t' read -r commit_message path_spec rest; do
  if [[ -n "${rest:-}" ]]; then
    echo "Invalid manifest line: expected exactly two tab-separated fields." >&2
    exit 1
  fi

  if [[ -z "${commit_message//[[:space:]]/}" ]]; then
    continue
  fi

  if [[ ${commit_message:0:1} == "#" ]]; then
    continue
  fi

  if [[ -z "${path_spec//[[:space:]]/}" ]]; then
    echo "Missing path list for batch: $commit_message" >&2
    exit 1
  fi

  # shellcheck disable=SC2206
  paths=( $path_spec )
  batch_num=$((batch_num + 1))

  echo
  echo "==> Batch $batch_num: $commit_message"
  echo "    Paths: ${paths[*]}"

  if git diff --quiet -- "${paths[@]}" && git diff --cached --quiet -- "${paths[@]}"; then
    echo "    No changes detected for this batch. Skipping."
    continue
  fi

  run_qa

  git add -- "${paths[@]}"

  if git diff --cached --quiet -- "${paths[@]}"; then
    echo "    Nothing staged after git add. Skipping commit."
    continue
  fi

  git commit -m "$commit_message"
  git push origin HEAD:main
done < "$manifest"

echo
echo "Finished processing manifest: $manifest"
