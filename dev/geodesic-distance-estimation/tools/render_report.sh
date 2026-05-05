#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 input.md [output-stem]" >&2
    exit 2
fi

input_md="$1"
if [ ! -f "$input_md" ]; then
    echo "Input Markdown file not found: $input_md" >&2
    exit 2
fi

if [ "$#" -eq 2 ]; then
    output_stem="$2"
else
    output_stem="${input_md%.md}"
fi

output_html="${output_stem}.html"
output_pdf="${output_stem}.pdf"
chrome_pdf="${output_stem}.chrome.pdf"
gs_tmp="${output_stem}.print-safe.tmp.pdf"
title="$(awk '/^# / {sub(/^# /, ""); print; exit}' "$input_md")"
if [ -z "$title" ]; then
    title="$(basename "$output_stem" | sed 's/_/ /g')"
fi

pandoc "$input_md" \
    --standalone \
    --metadata "title=${title}" \
    --from=markdown+tex_math_dollars+tex_math_single_backslash \
    --to=html5 \
    --mathjax \
    --toc \
    --toc-depth=3 \
    -o "$output_html"

chrome_bin="${CHROME_BIN:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"
if [ ! -x "$chrome_bin" ]; then
    echo "Chrome executable not found. Set CHROME_BIN or install Google Chrome." >&2
    exit 1
fi

abs_html="$(cd "$(dirname "$output_html")" && pwd)/$(basename "$output_html")"
"$chrome_bin" \
    --headless \
    --disable-gpu \
    --no-sandbox \
    --run-all-compositor-stages-before-draw \
    --virtual-time-budget=10000 \
    --print-to-pdf="$chrome_pdf" \
    "file://${abs_html}"

if ! command -v gs >/dev/null 2>&1; then
    echo "Ghostscript is required for print-safe PDF output." >&2
    echo "Chrome PDF was left at: $chrome_pdf" >&2
    exit 1
fi

gs \
    -dSAFER \
    -dBATCH \
    -dNOPAUSE \
    -sDEVICE=pdfwrite \
    -dCompatibilityLevel=1.4 \
    -dWriteObjStms=false \
    -dWriteXRefStm=false \
    -dPDFSETTINGS=/prepress \
    -dEmbedAllFonts=true \
    -dSubsetFonts=true \
    -dCompressFonts=true \
    -dDetectDuplicateImages=true \
    -dAutoRotatePages=/None \
    -sOutputFile="$gs_tmp" \
    "$chrome_pdf" >/dev/null

mv "$gs_tmp" "$output_pdf"
rm -f "$chrome_pdf"

echo "Wrote HTML: $output_html"
echo "Wrote print-safe PDF: $output_pdf"
