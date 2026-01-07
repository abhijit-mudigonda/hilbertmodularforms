#!/bin/bash

# Wrapper script for monitoring with scrolling support
# Usage: ./monitor-watch.sh

while true; do
    clear
    ./monitor-status-live.sh
    sleep 1
done
