#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TEX_FILE="transported_graph_hessian_trend_filtering.tex"
BUILD_INFO_TEX="$SCRIPT_DIR/transported_graph_hessian_trend_filtering_build_info.tex"

escape_for_tex() {
  local s="$1"
  s=${s//_/\\_}
  printf '%s\n' "$s"
}

BUILD_DATETIME="$(date '+%Y-%m-%d %H:%M:%S %Z')"

cat > "$BUILD_INFO_TEX" <<EOF
\renewcommand{\reportbuilddatetime}{$(escape_for_tex "$BUILD_DATETIME")}
EOF

cd "$SCRIPT_DIR"
pdflatex -interaction=nonstopmode "$TEX_FILE"
pdflatex -interaction=nonstopmode "$TEX_FILE"

printf '%s\n' "$SCRIPT_DIR/${TEX_FILE:r}.pdf"
