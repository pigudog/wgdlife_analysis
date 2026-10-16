#!/usr/bin/env bash
# Render the figure panels into results/. Run from anywhere: bash run.sh [figure2 figure3 figure4 figure5]
set -euo pipefail
cd "$(dirname "$0")"
figs=("$@"); [ ${#figs[@]} -eq 0 ] && figs=(figure2 figure3 figure4 figure5)
for f in "${figs[@]}"; do [ -d "$f" ] || { echo "no such figure directory: $f" >&2; exit 1; }; done
for f in "${figs[@]}"; do
  for s in "$f"/fig*.R; do
    echo "== $s"; Rscript "$s"
  done
done
echo "done -> results/"
