#!/bin/bash

# ANSI color codes
# Bright colors for new completions
BRIGHT_GREEN='\033[1;32m'
BRIGHT_YELLOW='\033[1;33m'
BRIGHT_RED='\033[1;31m'
BRIGHT_BLUE='\033[1;34m'
# Dim colors for old completions
DIM_GREEN='\033[2;32m'
DIM_YELLOW='\033[2;33m'
DIM_RED='\033[2;31m'
DIM_BLUE='\033[2;34m'
DIM_GRAY='\033[2;37m'
NC='\033[0m' # No Color

# Track when the parallel run started
RUN_START_FILE="parallel-logs/run_start_time.txt"

# Monitor currently running magma processes
echo "=== Currently Running Jobs ==="
echo ""

# Build quick lookup for labels in log files
declare -A log_status
if [ -f "logging/zero_sk.log" ]; then
    while IFS= read -r label; do
        [ -n "$label" ] && log_status["$label"]="zero"
    done < logging/zero_sk.log
fi
if [ -f "logging/nonzero_sk.log" ]; then
    while IFS= read -r label; do
        [ -n "$label" ] && log_status["$label"]="nonzero"
    done < logging/nonzero_sk.log
fi
if [ -f "logging/squarechecked_nonzero_sk.log" ]; then
    while IFS= read -r label; do
        [ -n "$label" ] && log_status["$label"]="squarechecked"
    done < logging/squarechecked_nonzero_sk.log
fi
if [ -f "logging/errors.log" ]; then
    while IFS= read -r label; do
        if [[ "$label" =~ ^-[0-9]+\.[0-9]+\.[0-9]+_ ]]; then
            log_status["$label"]="error"
        fi
    done < logging/errors.log
fi

# Get list of running magma processes with their chi_labels, start times, and memory usage
ps -eo pid,etime,rss,cmd | grep "magma.exe -b chi_label" | grep -v grep | while read pid elapsed rss rest; do
    # Extract chi_label from the rest of the command line
    chi_label=$(echo "$rest" | sed -n 's/.*chi_label:=\([^ ]*\).*/\1/p')
    
    # Check if this label already exists in a log file
    status_note=""
    if [ -n "${log_status[$chi_label]}" ]; then
        status_note=" (${log_status[$chi_label]})"
    fi
    
    # Convert RSS (in KB) to MB or GB for display
    mem_mb=$((rss / 1024))
    if [ $mem_mb -ge 1024 ]; then
        mem_gb=$(echo "scale=1; $mem_mb / 1024" | bc)
        mem_display="${mem_gb}G"
    else
        mem_display="${mem_mb}M"
    fi
    
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
    
    printf "[PID %6s] %-35s %5ds %6s%s\n" "$pid" "$chi_label" "$total_secs" "$mem_display" "$status_note"
done | sort -k3

echo ""

# Show completion stats if joblog exists
LOGDIR="parallel-logs"
if [ -f "$LOGDIR/joblog.txt" ]; then
    completed=$(tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' '$7 != "" {count++} END {print count+0}')
    running=$(ps aux | grep "magma.exe -b chi_label" | grep -v grep | wc -l)
    # Get total from joblog - only count lines that have been started (have a starttime)
    total=$(tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' '$3 != "" {count++} END {print count+0}')
    if [ $total -gt 0 ]; then
        pending=$((total - completed - running))
        # Ensure pending is not negative
        if [ $pending -lt 0 ]; then
            pending=0
        fi
        echo "Progress: $completed completed, $running running, $pending pending (total: $total)"
    else
        echo "Progress: $running running"
    fi
else
    running=$(ps aux | grep "magma.exe -b chi_label" | grep -v grep | wc -l)
    echo "Progress: $running running"
fi

echo ""
echo "=== Completed Jobs ==="
echo ""

# Build sets of labels from each log file
declare -A zero_labels
declare -A nonzero_labels
declare -A squarechecked_labels
declare -A error_labels
declare -A all_completed

# Get the run start time
RUN_START_TIME=0
if [ -f "$RUN_START_FILE" ]; then
    RUN_START_TIME=$(cat "$RUN_START_FILE" 2>/dev/null || echo "0")
else
    # If no run start time file exists, create one with current time
    # This way, any jobs that complete from now on will be bolded
    date +%s > "$RUN_START_FILE" 2>/dev/null || true
    RUN_START_TIME=$(cat "$RUN_START_FILE" 2>/dev/null || echo "0")
fi

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

# Read squarechecked_nonzero_sk.log
if [ -f "logging/squarechecked_nonzero_sk.log" ] && [ -s "logging/squarechecked_nonzero_sk.log" ]; then
    while read label; do
        squarechecked_labels["$label"]=1
        all_completed["$label"]=1
    done < logging/squarechecked_nonzero_sk.log
fi

# Read errors.log (only lines that look like chi_labels)
if [ -f "logging/errors.log" ] && [ -s "logging/errors.log" ]; then
    while read label; do
        # Only process lines that look like chi_labels (start with - and contain underscores and dots)
        if [[ "$label" =~ ^-[0-9]+\.[0-9]+\.[0-9]+_ ]]; then
            error_labels["$label"]=1
            all_completed["$label"]=1
        fi
    done < logging/errors.log
fi

# Get completed jobs from joblog and track which are new (completed after run start)
declare -A is_new_completion
if [ -f "$LOGDIR/joblog.txt" ] && [ "$RUN_START_TIME" -gt 0 ]; then
    while IFS=$'\t' read -r seq host starttime runtime send receive exitval signal command rest; do
        if [ -n "$exitval" ] && [ -n "$starttime" ]; then
            label=$(echo "$command" | sed -n 's/.*chi_label:=\([^ ]*\).*/\1/p')
            if [ -n "$label" ]; then
                all_completed["$label"]=1
                # Compare full starttime (with decimals) to run start time
                # Use bc for floating point comparison
                is_new=$(echo "$starttime >= $RUN_START_TIME" | bc)
                if [ "$is_new" = "1" ]; then
                    is_new_completion["$label"]=1
                fi
            fi
        fi
    done < <(tail -n +2 "$LOGDIR/joblog.txt")
fi

# Display completed jobs with appropriate colors (blue, green, red, yellow)
# Sorted by the number after the last underscore
# Bold new jobs that completed after run started
if [ ${#all_completed[@]} -eq 0 ]; then
    echo "No completed jobs yet"
else
    # First show blue (squarechecked)
    for label in "${!all_completed[@]}"; do
        if [ -n "${squarechecked_labels[$label]}" ]; then
            # Extract the number after the last underscore before the first dot for sorting
            # Handles both -5.0.1_XXX and -2.0.1_XXX formats
            sort_key=$(echo "$label" | sed -n 's/.*_\([0-9]*\)\..*/\1/p')
            # Default to 0 if extraction fails
            sort_key=${sort_key:-0}
            # Check if this is a new completion (completed after run started)
            if [ -n "${is_new_completion[$label]}" ]; then
                printf "0 %010d ${BRIGHT_BLUE}%s${NC} (squarechecked)\n" "$sort_key" "$label"
            else
                printf "0 %010d ${DIM_BLUE}%s ${DIM_GRAY}(squarechecked)${NC}\n" "$sort_key" "$label"
            fi
        fi
    done
    
    # Then show green (nonzero)
    for label in "${!all_completed[@]}"; do
        if [ -n "${nonzero_labels[$label]}" ] && [ -z "${squarechecked_labels[$label]}" ]; then
            sort_key=$(echo "$label" | sed -n 's/.*_\([0-9]*\)\..*/\1/p')
            sort_key=${sort_key:-0}
            if [ -n "${is_new_completion[$label]}" ]; then
                printf "1 %010d ${BRIGHT_GREEN}%s${NC} (nonzero)\n" "$sort_key" "$label"
            else
                printf "1 %010d ${DIM_GREEN}%s ${DIM_GRAY}(nonzero)${NC}\n" "$sort_key" "$label"
            fi
        fi
    done
    
    # Then show red (errors and missing)
    for label in "${!all_completed[@]}"; do
        if [ -n "${error_labels[$label]}" ]; then
            sort_key=$(echo "$label" | sed -n 's/.*_\([0-9]*\)\..*/\1/p')
            sort_key=${sort_key:-0}
            if [ -n "${is_new_completion[$label]}" ]; then
                printf "2 %010d ${BRIGHT_RED}%s${NC} (error)\n" "$sort_key" "$label"
            else
                printf "2 %010d ${DIM_RED}%s ${DIM_GRAY}(error)${NC}\n" "$sort_key" "$label"
            fi
        elif [ -z "${zero_labels[$label]}" ] && [ -z "${nonzero_labels[$label]}" ] && [ -z "${squarechecked_labels[$label]}" ]; then
            sort_key=$(echo "$label" | sed -n 's/.*_\([0-9]*\)\..*/\1/p')
            sort_key=${sort_key:-0}
            if [ -n "${is_new_completion[$label]}" ]; then
                printf "2 %010d ${BRIGHT_RED}%s${NC} (missing)\n" "$sort_key" "$label"
            else
                printf "2 %010d ${DIM_RED}%s ${DIM_GRAY}(missing)${NC}\n" "$sort_key" "$label"
            fi
        fi
    done
    
    # Finally show yellow (zero)
    for label in "${!all_completed[@]}"; do
        if [ -n "${zero_labels[$label]}" ]; then
            sort_key=$(echo "$label" | sed -n 's/.*_\([0-9]*\)\..*/\1/p')
            sort_key=${sort_key:-0}
            if [ -n "${is_new_completion[$label]}" ]; then
                printf "3 %010d ${BRIGHT_YELLOW}%s${NC} (zero)\n" "$sort_key" "$label"
            else
                printf "3 %010d ${DIM_YELLOW}%s ${DIM_GRAY}(zero)${NC}\n" "$sort_key" "$label"
            fi
        fi
    done
fi | sort -k1,1n -k2,2nr | sed 's/^[0-9]* [0-9]* //' | {
    # Read into array and display in two columns
    mapfile -t lines
    total=${#lines[@]}
    half=$(( (total + 1) / 2 ))
    
    for ((i=0; i<half; i++)); do
        left="${lines[i]}"
        right_idx=$((i + half))
        if [ $right_idx -lt $total ]; then
            right="${lines[right_idx]}"
            # Calculate visible length by removing ANSI codes
            left_visible=$(echo -e "$left" | sed 's/\x1b\[[0-9;]*m//g')
            left_len=${#left_visible}
            # Pad to 60 characters of visible text
            padding=$((60 - left_len))
            if [ $padding -gt 0 ]; then
                printf "%s%*s%s\n" "$left" "$padding" "" "$right"
            else
                printf "%s %s\n" "$left" "$right"
            fi
        else
            echo -e "$left"
        fi
    done
}
