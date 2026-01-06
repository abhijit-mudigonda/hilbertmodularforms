#!/bin/bash

JOBS=6
LOGDIR="parallel-logs"
mkdir -p "$LOGDIR"

cleanup() {
    tput cnorm
    exit
}
trap cleanup EXIT INT TERM

tput civis

parallel --line-buffer -j "$JOBS" \
    --joblog "$LOGDIR/joblog.txt" \
    --results "$LOGDIR" \
    --tagstring '{= s/.*chi_label:=(\S+)\s.*/\1/ =}' \
    < label-commands.sh &

PARALLEL_PID=$!

while kill -0 $PARALLEL_PID 2>/dev/null; do
    clear
    echo "=== Parallel Job Monitor ($(date +%H:%M:%S)) ==="
    echo ""
    
    if [ -f "$LOGDIR/joblog.txt" ]; then
        tail -n +2 "$LOGDIR/joblog.txt" | while IFS=$'\t' read -r seq host starttime runtime send receive exitval signal command rest; do
            if [ -z "$exitval" ] || [ "$exitval" = "0" ]; then
                chi_label=$(echo "$command" | sed -n 's/.*chi_label:=\(\S\+\).*/\1/p')
                
                if [ -n "$runtime" ] && [ "$runtime" != "0" ]; then
                    status="RUNNING"
                    elapsed="${runtime}s"
                elif [ -n "$exitval" ]; then
                    status="DONE"
                    elapsed="${runtime}s"
                else
                    status="RUNNING"
                    if [ -n "$starttime" ]; then
                        current_time=$(date +%s)
                        elapsed=$((current_time - ${starttime%.*}))
                        elapsed="${elapsed}s"
                    else
                        elapsed="0s"
                    fi
                fi
                
                printf "Job %2d: %-40s [%s] %s\n" "$seq" "$chi_label" "$status" "$elapsed"
            fi
        done | tail -n "$JOBS"
    fi
    
    sleep 1
done

wait $PARALLEL_PID
tput cnorm
echo ""
echo "All jobs completed!"
