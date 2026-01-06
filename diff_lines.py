#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("Usage: python diff_lines.py <file_A> <file_B>")
    sys.exit(1)

file_a = sys.argv[1]
file_b = sys.argv[2]

with open(file_b, 'r') as f:
    lines_b = set(line.rstrip('\n') for line in f)

with open(file_a, 'r') as f:
    for line in f:
        line_stripped = line.rstrip('\n')
        if line_stripped not in lines_b:
            print(line_stripped)
