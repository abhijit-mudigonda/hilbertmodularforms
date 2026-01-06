#!/bin/bash

JOBS=6
LOGDIR="parallel-logs"
mkdir -p "$LOGDIR"

# Run parallel with joblog (using --group instead of --line-buffer for better joblog updates)
parallel --group -j "$JOBS" \
    --joblog "$LOGDIR/joblog.txt" \
    --tagstring '{= s/.*chi_label:=(\S+)\s.*/\1/ =}' \
    < label-commands.sh
