#!/usr/bin/env bash
# Copy the JSON-only bundles from a sibling darboux-igph checkout into
# this repo's artifacts_out/ so GitHub Pages can serve them.
#
# Usage:  ./scripts/sync_artifacts.sh [path/to/darboux-igph]
# Default source: ../darboux-igph

set -euo pipefail

SRC_ROOT="${1:-../darboux-igph}"
SRC="$SRC_ROOT/artifacts_out"
DST="$(cd "$(dirname "$0")/.." && pwd)/artifacts_out"

if [ ! -d "$SRC" ]; then
  echo "error: $SRC does not exist." >&2
  echo "Run \`darboux-run\` in $SRC_ROOT first, then \`darboux-export-web artifacts_out\`." >&2
  exit 1
fi

if [ ! -f "$SRC/manifest.json" ]; then
  echo "error: $SRC/manifest.json missing. Run \`darboux-export-web artifacts_out\` in $SRC_ROOT." >&2
  exit 1
fi

mkdir -p "$DST"
cp "$SRC/manifest.json" "$DST/manifest.json"

count=0
for d in "$SRC"/*/; do
  [ -d "$d" ] || continue
  run="$(basename "$d")"
  bundle="$d/web_bundle.json"
  if [ -f "$bundle" ]; then
    mkdir -p "$DST/$run"
    cp "$bundle" "$DST/$run/web_bundle.json"
    count=$((count + 1))
  fi
done

echo "Synced $count run bundles + manifest into $DST"
