#!/usr/bin/env bash
# ============================================================================
# Generic Shell Templates — copy/paste friendly utilities and examples
# Location: cran_docs/cheatsheets/shell_templates.sh
#
# Usage patterns:
#   1) Run commands directly by copy/paste (comments are always on their own lines).
#   2) Or source this file and call helper functions:
#        source shell_templates.sh
#        help
#        git_tag 0.3.1
#
# Conventions:
#   - Strict mode for safer scripts.
#   - No inline comments on command lines. Explanations are above each command.
#   - Defaults are conservative. Adjust paths and names as needed.
# ============================================================================

set -euo pipefail

# -----------------------------
# Top-level helper (function)
# -----------------------------
help() {
  cat <<'EOF'
Available sections (search in this file):
  [PKG]   Package dev (R build/check, Rcpp glue)
  [GIT]   Git workflows (status, branches, tags, bisect)
  [FS]    Filesystem ops (find, rsync, tar, zip, permissions)
  [PROC]  Processes (ps, kill, nohup, nice)
  [NET]   Networking (curl, wget, ports, ssh, scp)
  [TEXT]  Text processing (grep, sed, awk, sort, uniq, tr, split, paste)
  [JSON]  JSON with jq
  [SYS]   System info (disk, memory, cpu)
  [DOCKER]Docker basics (images, containers)
  [PY]    Python venv and pip
  [MISC]  Date/time, checksum, random data, temporary dirs

Helper functions (call after 'source shell_templates.sh'):
  help
  git_tag VERSION
  git_undo_last_commit_soft
  r_build
  r_check_as_cran
  r_valgrind_script path/to/script.R
  link_many target_dir file1 file2 ...
EOF
}

# ---------------------------------------------------------------------------
# [PKG] Package development templates (R-focused)
# ---------------------------------------------------------------------------

r_load_all() {
  Rscript -e 'devtools::load_all()'
}

r_compile_attributes() {
  Rscript -e 'Rcpp::compileAttributes()'
}

r_build() {
  R CMD build .
}

r_check_latest() {
  local tgz
  tgz=$(ls -t ./*.tar.gz | head -n1)
  R CMD check "${tgz}"
}

r_check_as_cran() {
  local tgz
  tgz=$(ls -t ./*.tar.gz | head -n1)
  R CMD check --as-cran "${tgz}"
}

r_valgrind_script() {
  local script="${1:-}"
  if [[ -z "${script}" ]]; then
    echo "Usage: r_valgrind_script path/to/script.R" >&2
    return 2
  fi
  R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes" -f "${script}"
}

clean_check_dir() {
  rm -rf ./*Rcheck/
}

clean_build_artifacts() {
  find src -name '*.o' -delete
  find src -name '*.so' -delete
  find . -name '*.dll' -delete
}

# ---------------------------------------------------------------------------
# [GIT] Git workflows
# ---------------------------------------------------------------------------

git_status_short() {
  git status -sb
  git log --oneline -n 10
}

git_tag() {
  local version="${1:-}"
  if [[ -z "${version}" ]]; then
    echo "Usage: git_tag VERSION" >&2
    return 2
  fi
  git tag -a "v${version}" -m "Release v${version}"
  git push --tags
}

git_undo_last_commit_soft() {
  git reset --soft HEAD~1
}

git_fixup_last() {
  git commit --amend --no-edit
  git push --force-with-lease
}

git_new_branch() {
  local name="${1:-feature/branch}"
  git switch -c "${name}"
}

git_find_large_files() {
  git rev-list --objects --all | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' |   awk '$1=="blob"{print $3,$4}' | sort -nr | head -n 20
}

# ---------------------------------------------------------------------------
# [FS] Filesystem operations
# ---------------------------------------------------------------------------

link_many() {
  local target_dir="${1:-}"
  shift || true
  if [[ -z "${target_dir}" || "$#" -eq 0 ]]; then
    echo "Usage: link_many TARGET_DIR FILE1 [FILE2 ...]" >&2
    return 2
  fi
  mkdir -p "${target_dir}"
  for f in "$@"; do
    ln -s "$(realpath -m "${f}")" "${target_dir}/$(basename "${f}")"
  done
}

# Find files by name under current directory
# (example: find all .cpp files)
# find . -type f -name '*.cpp'

# Find files modified in last 2 days
# find . -type f -mtime -2

# Copy preserving attrs
# cp -a src/ dest/

# Move with rsync merge semantics (directory into directory)
rsync_merge_dir() {
  local src="${1:-}"
  local dst="${2:-}"
  if [[ -z "${src}" || -z "${dst}" ]]; then
    echo "Usage: rsync_merge_dir SRC DST" >&2
    return 2
  fi
  mkdir -p "${dst}"
  rsync -a "${src}/" "${dst}/"
}

# Create tar.gz archive of current directory (excluding VCS)
tar_cwd() {
  local name="${1:-archive}"
  tar --exclude-vcs -czf "${name}.tar.gz" .
}

# Extract tar.gz
# tar -xzf file.tar.gz

# Zip directory recursively
zip_dir() {
  local src="${1:-}"
  local out="${2:-archive.zip}"
  if [[ -z "${src}" ]]; then
    echo "Usage: zip_dir SRC_DIR [OUT.zip]" >&2
    return 2
  fi
  (cd "${src}" && zip -r "../${out}" .)
}

# Unzip
# unzip archive.zip -d outdir

# Permissions
# chmod -R u+rwX,go-rwx path
# chown -R user:group path

# ---------------------------------------------------------------------------
# [PROC] Processes
# ---------------------------------------------------------------------------

# List processes with tree
# pstree -p

# Find process on a port (Linux/macOS)
# lsof -i :8787

# Kill by name
# pkill -f 'process-name'

# Background process with nohup
# nohup long_task.sh > long_task.log 2>&1 &

# Lower priority run
# nice -n 10 command args

# ---------------------------------------------------------------------------
# [NET] Networking / SSH
# ---------------------------------------------------------------------------

# HTTP GET with headers
# curl -i https://example.com

# Download to file
# wget -O out.bin https://example.com/file.bin

# Test TCP port
# nc -zv localhost 5432

# SSH
# ssh -i ~/.ssh/id_ed25519 user@host

# Copy file via scp
# scp -i ~/.ssh/id_ed25519 local.file user@host:/remote/path/

# Port forward local:8080 -> remote:80
# ssh -L 8080:localhost:80 user@host

# ---------------------------------------------------------------------------
# [TEXT] Text processing
# ---------------------------------------------------------------------------

# Grep recursive for pattern
# rg -n 'pattern' .
# or: grep -RIn 'pattern' .

# Sed in-place replace (backup .bak)
# sed -i.bak 's/old/new/g' file.txt

# Awk: print 1st and 3rd columns
# awk '{print $1, $3}' file.txt

# Sort unique count
# sort file.txt | uniq -c | sort -nr

# Translate tabs to commas
# tr '\t' ',' < in.tsv > out.csv

# Split file into chunks of 1000 lines
# split -l 1000 big.txt chunk_

# Paste columns from files side-by-side
# paste file1 file2 > merged.txt

# ---------------------------------------------------------------------------
# [JSON] jq basics
# ---------------------------------------------------------------------------

# Pretty print JSON
# jq . data.json

# Extract field
# jq -r '.items[].name' data.json

# Filter by condition
# jq '.[] | select(.size > 100)' data.json

# ---------------------------------------------------------------------------
# [SYS] System info
# ---------------------------------------------------------------------------

# Disk usage of current dir
# du -sh .

# Top largest items
# du -sh * | sort -hr | head -n 20

# Memory and CPU (Linux)
# free -h
# nproc

# Uptime
# uptime

# ---------------------------------------------------------------------------
# [DOCKER] Docker essentials
# ---------------------------------------------------------------------------

# List images
# docker images

# List containers
# docker ps -a

# Build image
# docker build -t myimage:latest .

# Run with volume
# docker run --rm -it -v "$(pwd)":/work -w /work myimage:latest bash

# ---------------------------------------------------------------------------
# [PY] Python virtualenv + pip
# ---------------------------------------------------------------------------

# Create venv
# python3 -m venv .venv

# Activate venv
# source .venv/bin/activate

# Upgrade pip and install
# python -m pip install --upgrade pip
# pip install -r requirements.txt

# Freeze deps
# pip freeze > requirements.txt

# ---------------------------------------------------------------------------
# [MISC] Date/time, checksums, random, temp dirs
# ---------------------------------------------------------------------------

# ISO timestamp
# date -u +%Y-%m-%dT%H:%M:%SZ

# MD5 / SHA256
# md5sum file.bin
# shasum -a 256 file.bin

# Random 32 bytes hex
# openssl rand -hex 32

# Temp dir pattern
# tmpdir=$(mktemp -d)
# echo "${tmpdir}"
# rm -rf "${tmpdir}"
