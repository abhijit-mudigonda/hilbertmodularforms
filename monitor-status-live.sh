#!/bin/bash

# ANSI color codes
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Monitor currently running magma processes
echo "=== Currently Running Jobs ==="
echo ""

# Get list of running magma processes with their chi_labels and start times
ps -eo pid,etime,cmd | grep "magma.exe -b chi_label" | grep -v grep | while read pid elapsed rest; do
    # Extract chi_label from the rest of the command line
    chi_label=$(echo "$rest" | sed -n 's/.*chi_label:=\([^ ]*\).*/\1/p')
    
    # Convert elapsed time to seconds for consistent display
    # etime format can be: SS, MM:SS, HH:MM:SS, or DD-HH:MM:SS
    # Strip leading zeros to avoid octal interpretation
    if [[ $elapsed =~ ^([0-9]+)-([0-9]+):([0-9]+):([0-9]+)$ ]]; then
        # DD-HH:MM:SS
        days=$((10#${BASH_REMATCH[1]}))
        hours=$((10#${BASH_REMATCH[2]}))
        mins=$((10#${BASH_REMATCH[3]}))
        secs=$((10#${BASH_REMATCH[4]}))
        total_secs=$((days * 86400 + hours * 3600 + mins * 60 + secs))
    elif [[ $elapsed =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
        # HH:MM:SS
        hours=$((10#${BASH_REMATCH[1]}))
        mins=$((10#${BASH_REMATCH[2]}))
        secs=$((10#${BASH_REMATCH[3]}))
        total_secs=$((hours * 3600 + mins * 60 + secs))
    elif [[ $elapsed =~ ^([0-9]+):([0-9]+)$ ]]; then
        # MM:SS
        mins=$((10#${BASH_REMATCH[1]}))
        secs=$((10#${BASH_REMATCH[2]}))
        total_secs=$((mins * 60 + secs))
    else
        # Just seconds
        total_secs=$((10#$elapsed))
    fi
    
    printf "[PID %6s] %-40s %5ds\n" "$pid" "$chi_label" "$total_secs"
done | sort -k3

echo ""

# Show completion stats if joblog exists
LOGDIR="parallel-logs"
if [ -f "$LOGDIR/joblog.txt" ]; then
    completed=$(tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' '$7 != "" {count++} END {print count+0}')
    running=$(ps aux | grep "magma.exe -b chi_label" | grep -v grep | wc -l)
    total=$(wc -l < label-commands.sh)
    echo "Progress: $completed completed, $running running, $((total - completed - running)) pending (total: $total)"
else
    running=$(ps aux | grep "magma.exe -b chi_label" | grep -v grep | wc -l)
    total=$(wc -l < label-commands.sh)
    echo "Progress: $running running (total: $total)"
fi

echo ""
echo "=== Completed Jobs ==="
echo ""

# Build sets of labels from each log file
declare -A zero_labels
declare -A nonzero_labels
declare -A error_labels
declare -A all_completed

# Read zero_sk.log
if [ -f "logging/zero_sk.log" ] && [ -s "logging/zero_sk.log" ]; then
    while read label; do
        zero_labels["$label"]=1
        all_completed["$label"]=1
    done < logging/zero_sk.log
fi

# Read nonzero_sk.log
if [ -f "logging/nonzero_sk.log" ] && [ -s "logging/nonzero_sk.log" ]; then
    while read label; do
        nonzero_labels["$label"]=1
        all_completed["$label"]=1
    done < logging/nonzero_sk.log
fi

# Read errors.log
if [ -f "logging/errors.log" ] && [ -s "logging/errors.log" ]; then
    while read label; do
        error_labels["$label"]=1
        all_completed["$label"]=1
    done < logging/errors.log
fi

# Get completed jobs from joblog
if [ -f "$LOGDIR/joblog.txt" ]; then
    tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' '$7 != "" {print $0}' | while IFS=$'\t' read -r seq host starttime runtime send receive exitval signal command rest; do
        label=$(echo "$command" | sed -n 's/.*chi_label:=\([^ ]*\).*/\1/p')
        if [ -n "$label" ]; then
            all_completed["$label"]=1
        fi
    done
fi

# Display completed jobs with appropriate colors (green, then red, then yellow)
if [ ${#all_completed[@]} -eq 0 ]; then
    echo "No completed jobs yet"
else
    # First show green (nonzero)
    for label in "${!all_completed[@]}"; do
        if [ -n "${nonzero_labels[$label]}" ]; then
            echo -e "1${GREEN}${label}${NC} (nonzero)"
        fi
    done
    
    # Then show red (errors and missing)
    for label in "${!all_completed[@]}"; do
        if [ -n "${error_labels[$label]}" ]; then
            echo -e "2${RED}${label}${NC} (error)"
        elif [ -z "${zero_labels[$label]}" ] && [ -z "${nonzero_labels[$label]}" ]; then
            echo -e "2${RED}${label}${NC} (missing)"
        fi
    done
    
    # Finally show yellow (zero)
    for label in "${!all_completed[@]}"; do
        if [ -n "${zero_labels[$label]}" ]; then
            echo -e "3${YELLOW}${label}${NC} (zero)"
        fi
    done
fi | sort | sed 's/^[0-9]//'
