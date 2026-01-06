#!/bin/bash

JOBS=6
LOGDIR="parallel-logs"
mkdir -p "$LOGDIR"

cleanup() {
    tput cnorm
    stty echo
    exit
}
trap cleanup EXIT INT TERM

tput civis
stty -echo

# Start parallel in background with joblog
parallel --line-buffer -j "$JOBS" \
    --joblog "$LOGDIR/joblog.txt" \
    --results "$LOGDIR" \
    --tagstring '{= s/.*chi_label:=(\S+)\s.*/\1/ =}' \
    < label-commands.sh > "$LOGDIR/combined.log" 2>&1 &

PARALLEL_PID=$!

# Monitor loop
while kill -0 $PARALLEL_PID 2>/dev/null; do
    tput cup 0 0
    tput ed
    
    echo "=== Parallel Job Monitor - $(date +%H:%M:%S) ==="
    echo ""
    
    # Initialize slot array
    declare -A slots
    for i in $(seq 1 $JOBS); do
        slots[$i]="[Slot $i] Idle"
    done
    
    if [ -f "$LOGDIR/joblog.txt" ]; then
        # Parse joblog to find running jobs
        tail -n +2 "$LOGDIR/joblog.txt" | while IFS=$'\t' read -r seq host starttime runtime send receive exitval signal command rest; do
            # Running job: exitval is empty and starttime exists
            if [ -z "$exitval" ] && [ -n "$starttime" ]; then
                chi_label=$(echo "$command" | sed -n 's/.*chi_label:=\(\S\+\).*/\1/p')
                current_time=$(date +%s.%N)
                elapsed=$(echo "$current_time - $starttime" | bc)
                elapsed_int=${elapsed%.*}
                
                # Determine which slot (simple round-robin based on seq)
                slot=$(( (seq - 1) % JOBS + 1 ))
                
                printf "[Slot %d] %s (running %ds)\n" "$slot" "$chi_label" "$elapsed_int"
            fi
        done | tail -n $JOBS
    fi
    
    # Show completed count
    if [ -f "$LOGDIR/joblog.txt" ]; then
        completed=$(tail -n +2 "$LOGDIR/joblog.txt" | awk -F'\t' '$7 != "" {count++} END {print count+0}')
        total=$(wc -l < label-commands.sh)
        echo ""
        echo "Progress: $completed / $total completed"
    fi
    
    sleep 0.5
done

wait $PARALLEL_PID
tput cnorm
stty echo
tput ed
echo ""
echo "All jobs completed!"
