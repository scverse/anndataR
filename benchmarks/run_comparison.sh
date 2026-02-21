#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Running performance benchmark with installed package..."
Rscript "$SCRIPT_DIR/performance_benchmark.R" --source installed \
  > "$SCRIPT_DIR/results_before.txt"

echo "Running performance benchmark with local source..."
Rscript "$SCRIPT_DIR/performance_benchmark.R" --source local \
  > "$SCRIPT_DIR/results_after.txt"

echo ""
echo "Results written to:"
echo "  $SCRIPT_DIR/results_before.txt"
echo "  $SCRIPT_DIR/results_after.txt"
