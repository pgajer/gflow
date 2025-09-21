#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Total timer
# -----------------------------
SCRIPT_START="$(date +%s)"

# -----------------------------
# Helpers
# -----------------------------
format_time() {
    local seconds="$1"
    if (( seconds < 60 )); then
        echo "${seconds}s"
    else
        echo "$((seconds / 60))m $((seconds % 60))s"
    fi
}

show_elapsed() {
    local start="$1"
    local label="$2"
    local end
    end="$(date +%s)"
    local elapsed=$((end - start))
    echo "⏱️  ${label}: $(format_time "$elapsed")" | tee -a "$LOG"
}

# -----------------------------
# Configuration
# -----------------------------
PACKAGE_NAME="gflow"
REMOTE_HOST="pawel@192.168.4.45"

# Path to R on the remote host (adjust if needed)
REMOTE_R_BIN="${REMOTE_R_BIN:-~/R-devel/bin/R}"

# Per-run output under build/valgrind/<UTC timestamp>, plus a 'latest' symlink
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOCAL_OUT_DIR="build/valgrind/${RUN_STAMP}"
LOG="${LOCAL_OUT_DIR}/valgrind_stdout.txt"  # local log file

# Rsync excludes
RSYNC_EXCLUDES=(
  "--exclude=.git"
  "--exclude=*.Rcheck"
  "--exclude=.Rproj.user"
  "--exclude=*.o"
  "--exclude=*.so"
  "--exclude=.DS_Store"
  "--exclude=__pycache__"
  "--exclude=build/"
)

echo "🔍 Valgrind check on remote host: ${REMOTE_HOST}  —  package: ${PACKAGE_NAME}"

# Sanity: run from package root
if [[ ! -f DESCRIPTION ]]; then
  echo "❌ Please run from the package root (no DESCRIPTION found)." >&2
  exit 2
fi

# Warn if working tree is dirty (we rsync what we see)
if [[ -n "$(git status --porcelain 2>/dev/null || true)" ]]; then
  echo "ℹ️  Working tree has uncommitted changes; we will rsync the current files as-is."
fi

# 1) Remote temp dir
REMOTE_TMP="$(ssh "${REMOTE_HOST}" 'mktemp -d')"
echo "📁 Remote temp: ${REMOTE_TMP}"

# Ensure the remote tmp gets cleaned even if the script fails later
cleanup() {
  local ec=$?
  if [[ -n "${REMOTE_TMP:-}" ]]; then
    ssh -o BatchMode=yes -o ConnectTimeout=5 "${REMOTE_HOST}" "rm -rf '${REMOTE_TMP}'" >/dev/null 2>&1 || true
  fi
  exit "$ec"
}
trap cleanup EXIT

# 2) Send snapshot
echo "📤 Sending snapshot to remote…"
RSYNC_START="$(date +%s)"

# Make sure local output dir exists before we first tee -a into $LOG
mkdir -p "${LOCAL_OUT_DIR}"

# Ensure remote target dir exists before rsync; then push
rsync -az --delete "${RSYNC_EXCLUDES[@]}" \
  --rsync-path="mkdir -p '${REMOTE_TMP}/${PACKAGE_NAME}' && rsync" \
  ./ "${REMOTE_HOST}:${REMOTE_TMP}/${PACKAGE_NAME}/"

show_elapsed "$RSYNC_START" "File transfer"

# 3) Build + Valgrind check on remote
echo "🏃 Running build + Valgrind check on remote…"
REMOTE_START="$(date +%s)"

ssh "${REMOTE_HOST}" bash -s -- "${REMOTE_TMP}" "${PACKAGE_NAME}" "${REMOTE_R_BIN}" << 'EOF'
set -euo pipefail
shopt -s nullglob

REMOTE_TMP="$1"
PACKAGE_NAME="$2"
R_BIN="$3"

# Where to store full stdout on the remote
LOG="$REMOTE_TMP/valgrind_stdout.txt"

# --- Build the package tarball in the snapshot dir ---
cd "$REMOTE_TMP/$PACKAGE_NAME"

echo "🛠️  Building package…" | tee "$LOG"
BUILD_START="$(date +%s)"
"$R_BIN" CMD build . | tee -a "$LOG"
echo "⏱️  Package build: $(($(date +%s) - BUILD_START))s" | tee -a "$LOG"

TARBALL="$(ls -t *.tar.gz | head -n1)"
echo "📦 Tarball: $TARBALL" | tee -a "$LOG"

# --- Prepare Valgrind options (close to what rhub valgrind environments use) ---
# You can tweak these via REMOTE_VALGRIND_OPTS in your environment if needed.
VALGRIND_OPTS_DEFAULT='--tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes --num-callers=30 --run-libc-freeres=no'
VALGRIND_OPTS="${REMOTE_VALGRIND_OPTS:-$VALGRIND_OPTS_DEFAULT}"

# R's default suppression file (if present)
RHOME="$("$R_BIN" RHOME)"
SUPP_FILE="$RHOME/etc/valgrind.supp"
if [[ -f "$SUPP_FILE" ]]; then
  VALGRIND_OPTS="$VALGRIND_OPTS --suppressions=$SUPP_FILE"
fi

# Export so R CMD check sees it
export R_VALGRIND_OPTS="$VALGRIND_OPTS"

# Optional: be kind to CRAN-like checks (you can drop/modify these)
export _R_CHECK_FORCE_SUGGESTS_=false
export _R_CHECK_CRAN_INCOMING_=false

# --- Run R CMD check under Valgrind ---
echo "🧪 Running: R CMD check --use-valgrind ${TARBALL}" | tee -a "$LOG"
CHECK_START="$(date +%s)"

# Notes:
#  - --use-valgrind: runs examples, tests, vignettes under Valgrind
#  - --no-manual: avoid LaTeX dependency on the remote, tweak as you like
#  - Add --as-cran if you want CRAN-flavor checks (slower)
"$R_BIN" CMD check \
  --use-valgrind \
  --no-manual \
  "$TARBALL" | tee -a "$LOG" || true

echo "⏱️  R CMD check (valgrind): $(($(date +%s) - CHECK_START))s" | tee -a "$LOG"

# Figure out the *.Rcheck dir name (R names it <pkg>.<version>.Rcheck)
RCHECK_DIR="$(ls -d ${PACKAGE_NAME}*.Rcheck | head -n1 || true)"
if [[ -z "${RCHECK_DIR}" || ! -d "${RCHECK_DIR}" ]]; then
  echo "❗ Could not find .Rcheck directory" | tee -a "$LOG"
else
  echo "📁 Rcheck dir: ${RCHECK_DIR}" | tee -a "$LOG"
fi

# Quick inline summary: pull out Valgrind "ERROR SUMMARY" lines if any
if [[ -n "${RCHECK_DIR:-}" && -d "${RCHECK_DIR}" ]]; then
  echo "" | tee -a "$LOG"
  echo "================= Valgrind ERROR SUMMARY (quick scan) =================" | tee -a "$LOG"
  # Look in tests/, examples/, and vignettes/ subtrees for valgrind logs
  # Grep won't fail the script if nothing is found
  grep -R --include='*.log' -n "ERROR SUMMARY" "$RCHECK_DIR" | tee -a "$LOG" || true
fi

# Collect artifacts for retrieval
mkdir -p "$REMOTE_TMP/artifacts"
cp -r "$RCHECK_DIR" "$REMOTE_TMP/artifacts/" 2>/dev/null || true
cp "$TARBALL" "$REMOTE_TMP/artifacts/" 2>/dev/null || true
EOF

show_elapsed "$REMOTE_START" "Total remote execution"

# 4) Pull logs/artifacts
echo "📥 Pulling logs/artifacts to ${LOCAL_OUT_DIR} …"
PULL_START="$(date +%s)"
rsync -az "${REMOTE_HOST}:${REMOTE_TMP}/" "${LOCAL_OUT_DIR}/"
show_elapsed "$PULL_START" "Results retrieval"
echo "ℹ️  Logs: ${LOCAL_OUT_DIR}/valgrind_stdout.txt"

# 5) Local quick summary of ERROR SUMMARY lines (if any)
echo "📊 Scanning Valgrind logs locally (if present)…" | tee -a "$LOG"
if compgen -G "${LOCAL_OUT_DIR}/artifacts/"'*'.Rcheck  > /dev/null; then
  for d in "${LOCAL_OUT_DIR}/artifacts/"*.Rcheck; do
    echo "— $(basename "$d") —" | tee -a "$LOG"
    grep -R --include='*.log' -n "ERROR SUMMARY" "$d" | tee -a "$LOG" || echo "  (no ERROR SUMMARY lines found)" | tee -a "$LOG"
  done
else
  echo "ℹ️  No .Rcheck directory found under artifacts/." | tee -a "$LOG"
fi

# 6) Create latest symlink
mkdir -p build/valgrind
ln -sfn "${RUN_STAMP}" build/valgrind/latest

# (Remote temp cleanup handled by trap)

# Show total time
echo ""
echo "✨ Done! Total time: $(format_time "$(( $(date +%s) - SCRIPT_START ))")"
echo "📋 Results in: build/valgrind/latest/"
