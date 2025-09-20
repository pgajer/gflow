#!/usr/bin/env bash
set -euo pipefail

# ---------- Configuration ----------
PACKAGE_NAME="gflow"
REMOTE_HOST="pawel@192.168.4.45"

# Per-run output under build/rchk/<UTC timestamp>, plus a 'latest' symlink
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOCAL_OUT_DIR="build/rchk/${RUN_STAMP}"

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
rsync -az --delete "${RSYNC_EXCLUDES[@]}" ./ "${REMOTE_HOST}:${REMOTE_TMP}/${PACKAGE_NAME}/"

# 3) Build + rchk on remote
echo "🏃 Running build + rchk on remote…"
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
~/R-devel/bin/R CMD build . | tee -a "$LOG"

TARBALL="$(ls -t *.tar.gz | head -n1)"
echo "📦 Tarball: $TARBALL" | tee -a "$LOG"

# --- Install the package to R-devel's library ---
echo "📦 Installing package..." | tee -a "$LOG"
# Install to the R-devel library path that rchk expects
mkdir -p ~/R-devel/packages/lib
~/R-devel/bin/R CMD INSTALL --library=~/R-devel/packages/lib "$TARBALL" | tee -a "$LOG"

# --- Run rchk checks ---
echo "🔬 Running rchk…" | tee -a "$LOG"
cd ~/R-devel
~/rchk/scripts/check_package.sh "$PACKAGE_NAME" | tee -a "$LOG" || true

echo "" | tee -a "$LOG"
echo "================= RCHK RESULTS =================" | tee -a "$LOG"
FOUND=0

# Look for check files in the correct location (where check_package.sh puts them)
for check in ~/R-devel/packages/lib/${PACKAGE_NAME}/libs/*.so.*check; do
  if [ -f "$check" ]; then
    FOUND=1
    echo "--- $(basename "$check") ---" | tee -a "$LOG"
    # Show actual content, filtering out the "too many states" messages
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

# 4) Pull logs/artifacts + latest symlink
if [[ -n "${LOCAL_OUT_DIR}" ]]; then
  mkdir -p "${LOCAL_OUT_DIR}"
  echo "📥 Pulling logs/artifacts to ${LOCAL_OUT_DIR} …"
  rsync -az "${REMOTE_HOST}:${REMOTE_TMP}/" "${LOCAL_OUT_DIR}/"
  echo "ℹ️  Logs: ${LOCAL_OUT_DIR}/rchk_stdout.txt"
  ln -sfn "${RUN_STAMP}" build/rchk/latest
fi

# 5) Cleanup remote temp
echo "🧹 Cleaning remote temp…"
ssh "${REMOTE_HOST}" "rm -rf '${REMOTE_TMP}'"

echo "✨ Done. Inspect build/rchk/latest/rchk_stdout.txt"
