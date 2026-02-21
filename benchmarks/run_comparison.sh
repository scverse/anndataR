#!/bin/bash
set -euo pipefail

echo "Running performance benchmark with installed package..."
Rscript "benchmarks/run_benchmark.R" \
  --source installed \
  --times 5 \
  > "benchmarks/results_before.txt"

echo "Running performance benchmark with local source..."
Rscript "benchmarks/run_benchmark.R" \
  --source local \
  --times 5 \
  > "benchmarks/results_after.txt"

echo ""
echo "Results written to:"
echo "  benchmarks/results_before.txt"
echo "  benchmarks/results_after.txt"
