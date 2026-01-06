#!/bin/bash

LOGDIR="parallel-logs"

if [ ! -f "$LOGDIR/joblog.txt" ]; then
    echo "No joblog found. Run the parallel job first."
    exit 1
fi

# Show header
echo "=== Currently Running Jobs ==="
echo ""

# Get current time
current_time=$(date +%s)

# Parse joblog and show running jobs
tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' -v now="$current_time" '
    $7 == "" && $3 != "" {
        chi_label = $0
        gsub(/.*chi_label:=/, "", chi_label)
        gsub(/[[:space:]].*/, "", chi_label)
        
        elapsed = now - int($3)
        
        printf "%-50s %4ds\n", chi_label, elapsed
    }
' | tail -n 6

echo ""
# Show completion stats
completed=$(tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' '$7 != "" {count++} END {print count+0}')
total=$(wc -l < label-commands.sh)
echo "Progress: $completed / $total completed"
