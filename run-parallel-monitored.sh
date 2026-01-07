#!/bin/bash

# Usage: ./run-parallel-monitored.sh <labels_file> [num_jobs] [magma_script]
# Example: ./run-parallel-monitored.sh 12-unstable-labels-5-1000.txt 6 label-search-12.m

if [ $# -lt 1 ]; then
    echo "Usage: $0 <labels_file> [num_jobs] [magma_script]"
    echo "  labels_file  : Text file with one chi_label per line"
    echo "  num_jobs     : Number of parallel jobs (default: 6)"
    echo "  magma_script : Magma script to run (default: label-search-12.m)"
    exit 1
fi

LABELS_FILE="$1"
JOBS="${2:-6}"
MAGMA_SCRIPT="${3:-label-search-12.m}"
LOGDIR="parallel-logs"

if [ ! -f "$LABELS_FILE" ]; then
    echo "Error: Labels file '$LABELS_FILE' not found"
    exit 1
fi

if [ ! -f "$MAGMA_SCRIPT" ]; then
    echo "Error: Magma script '$MAGMA_SCRIPT' not found"
    exit 1
fi

mkdir -p "$LOGDIR"

# Record the start time of this parallel run
date +%s > "$LOGDIR/run_start_time.txt"

# Run parallel with joblog, building commands from the labels file
cat "$LABELS_FILE" | parallel --group -j "$JOBS" \
    --joblog "$LOGDIR/joblog.txt" \
    --tagstring '{= s/.*chi_label:=(\S+)\s.*/\1/ =}' \
    "magma -b chi_label:={} $MAGMA_SCRIPT"
