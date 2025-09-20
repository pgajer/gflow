#!/usr/bin/env bash
set -euo pipefail

# Start total timer
SCRIPT_START=$(date +%s)

# Function to format elapsed time
format_time() {
    local seconds=$1
    if (( seconds < 60 )); then
        echo "${seconds}s"
    else
        echo "$((seconds / 60))m $((seconds % 60))s"
    fi
}

# Function to show elapsed time for a step
show_elapsed() {
    local start=$1
    local label=$2
    local end=$(date +%s)
    local elapsed=$((end - start))
    echo "⏱️  $label: $(format_time $elapsed)" | tee -a "$LOG"
}

# ---------- Configuration ----------
PACKAGE_NAME="gflow"
REMOTE_HOST="pawel@192.168.4.45"

# Per-run output under build/rchk/<UTC timestamp>, plus a 'latest' symlink
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOCAL_OUT_DIR="build/rchk/${RUN_STAMP}"
LOG="${LOCAL_OUT_DIR}/rchk_stdout.txt"  # Define LOG for local use

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

echo "🔍 rchk on remote host (no GitHub): ${REMOTE_HOST}  —  package: ${PACKAGE_NAME}"

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

# 2) Send snapshot
echo "📤 Sending snapshot to remote…"
RSYNC_START=$(date +%s)
rsync -az --delete "${RSYNC_EXCLUDES[@]}" ./ "${REMOTE_HOST}:${REMOTE_TMP}/${PACKAGE_NAME}/"
mkdir -p "${LOCAL_OUT_DIR}"  # Create local dir for logging
show_elapsed $RSYNC_START "File transfer"

# 3) Build + rchk on remote
echo "🏃 Running build + rchk on remote…"
REMOTE_START=$(date +%s)
ssh "${REMOTE_HOST}" bash -s -- "${REMOTE_TMP}" "${PACKAGE_NAME}" << 'EOF'
set -euo pipefail

REMOTE_TMP="$1"
PACKAGE_NAME="$2"

# Source rchk configuration
source ~/rchk/scripts/config.inc
source ~/rchk/scripts/cmpconfig.inc

# Where to store full stdout on the remote
LOG="$REMOTE_TMP/rchk_stdout.txt"

# --- Build the package tarball in the snapshot dir ---
cd "$REMOTE_TMP/$PACKAGE_NAME"

echo "🛠️  Building package…" | tee "$LOG"
BUILD_START=$(date +%s)
~/R-devel/bin/R CMD build . | tee -a "$LOG"
echo "⏱️  Package build: $(($(date +%s) - BUILD_START))s" | tee -a "$LOG"

TARBALL="$(ls -t *.tar.gz | head -n1)"
echo "📦 Tarball: $TARBALL" | tee -a "$LOG"

# --- Install the package to R-devel's library ---
echo "📦 Installing package..." | tee -a "$LOG"
mkdir -p ~/R-devel/packages/lib
INSTALL_START=$(date +%s)
~/R-devel/bin/R CMD INSTALL --library=~/R-devel/packages/lib "$TARBALL" | tee -a "$LOG"
echo "⏱️  Package install: $(($(date +%s) - INSTALL_START))s" | tee -a "$LOG"

# --- Run rchk checks ---
echo "🔬 Running rchk…" | tee -a "$LOG"
cd ~/R-devel
RCHK_START=$(date +%s)
~/rchk/scripts/check_package.sh "$PACKAGE_NAME" | tee -a "$LOG" || true
echo "⏱️  rchk analysis: $(($(date +%s) - RCHK_START))s" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "================= RCHK RESULTS =================" | tee -a "$LOG"
FOUND=0

# Look for check files in the correct location
for check in ~/R-devel/packages/lib/${PACKAGE_NAME}/libs/*.so.*check; do
  if [ -f "$check" ]; then
    FOUND=1
    echo "--- $(basename "$check") ---" | tee -a "$LOG"
    grep -v "ERROR: too many states" "$check" | grep -v "not found in module" | grep -v "Cannot find function R_gc_internal" | tee -a "$LOG"
    echo "" | tee -a "$LOG"
  fi
done

if [[ $FOUND -eq 0 ]]; then
  echo "ℹ️  No *.so.*check artifacts found for $PACKAGE_NAME." | tee -a "$LOG"
fi

# Copy check results to temp dir for retrieval
mkdir -p "$REMOTE_TMP/check_results"
cp ~/R-devel/packages/lib/${PACKAGE_NAME}/libs/*.so.*check "$REMOTE_TMP/check_results/" 2>/dev/null || true
EOF
show_elapsed $REMOTE_START "Total remote execution"

# 4) Pull logs/artifacts
echo "📥 Pulling logs/artifacts to ${LOCAL_OUT_DIR} …"
PULL_START=$(date +%s)
rsync -az "${REMOTE_HOST}:${REMOTE_TMP}/" "${LOCAL_OUT_DIR}/"
show_elapsed $PULL_START "Results retrieval"
echo "ℹ️  Logs: ${LOCAL_OUT_DIR}/rchk_stdout.txt"

# 5) Process rchk errors into organized reports
if [ -f "${LOCAL_OUT_DIR}/gflow/libsonly/gflow/libs/gflow.so.bcheck" ]; then
  echo "📊 Processing rchk errors..."
  python3 tools/process_rchk_errors.py "${LOCAL_OUT_DIR}/gflow/libsonly/gflow/libs/gflow.so.bcheck"

  # Move reports into the timestamped directory
  if [ -d "rchk_reports" ]; then
    mv rchk_reports "${LOCAL_OUT_DIR}/"
    echo "📁 Error reports saved in: ${LOCAL_OUT_DIR}/rchk_reports/"
  fi
fi

# 6) Create latest symlink
ln -sfn "${RUN_STAMP}" build/rchk/latest

# 7) Cleanup remote temp
echo "🧹 Cleaning remote temp…"
ssh "${REMOTE_HOST}" "rm -rf '${REMOTE_TMP}'"

# Show total time
echo ""
echo "✨ Done! Total time: $(format_time $(($(date +%s) - SCRIPT_START)))"
echo "📋 Results in: build/rchk/latest/"
