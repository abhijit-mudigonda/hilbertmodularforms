#!/usr/bin/env python3
"""
Find level norms x in no-ike-2.txt such that for every m,
x*2^m is either in both no-ike-2.txt and 12-all-labels-2.txt or in neither.
"""

import re

def extract_level_norm(line):
    """Extract the level norm (integer part before decimal) from a label."""
    # Pattern: -2.0.1_{x}.{y}_...
    match = re.search(r'-2\.0\.1_(\d+)\.\d+_', line)
    if match:
        return int(match.group(1))
    return None

def read_level_norms(filename):
    """Read all level norms from a file."""
    norms = set()
    with open(filename, 'r') as f:
        for line in f:
            norm = extract_level_norm(line)
            if norm is not None:
                norms.add(norm)
    return norms

# Read level norms from both files
no_ike_norms = read_level_norms('/home/abhijitm/hmf/no-ike-2.txt')
all_labels_norms = read_level_norms('/home/abhijitm/hmf/12-all-labels-2.txt')

print(f"Level norms in no-ike-2.txt: {len(no_ike_norms)}")
print(f"Level norms in 12-all-labels-2.txt: {len(all_labels_norms)}")
print()

# For each x in no-ike, check if x*2^m is consistent across both files
consistent_norms = []
inconsistent_norms = []

for x in sorted(no_ike_norms):
    # Check all powers of 2 times x up to a reasonable limit
    # We'll check up to 2^15 (which gives us numbers up to ~30000 for small x)
    is_consistent = True
    violations = []
    
    for m in range(1, 16):
        x_times_2m = x * (2 ** m)
        
        # Check if it's in both or in neither
        in_no_ike = x_times_2m in no_ike_norms
        in_all_labels = x_times_2m in all_labels_norms
        
        # Consistent means: both True or both False
        if in_no_ike != in_all_labels:
            is_consistent = False
            violations.append((m, x_times_2m, in_no_ike, in_all_labels))
    
    if is_consistent:
        consistent_norms.append(x)
    else:
        inconsistent_norms.append((x, violations))

print("Level norms x in no-ike-2.txt where x*2^m is CONSISTENT (in both or neither):")
print("=" * 70)

if consistent_norms:
    for x in consistent_norms:
        print(f"x = {x}")
else:
    print("No such cases found.")

print()
print(f"Total consistent: {len(consistent_norms)}")
print()

print("Level norms x in no-ike-2.txt where x*2^m is INCONSISTENT:")
print("=" * 70)

if inconsistent_norms:
    for x, violations in inconsistent_norms:
        print(f"\nx = {x}:")
        for m, x_times_2m, in_no_ike, in_all_labels in violations:
            no_ike_status = "IN no-ike" if in_no_ike else "NOT in no-ike"
            all_labels_status = "IN all-labels" if in_all_labels else "NOT in all-labels"
            print(f"  m={m}: x*2^{m}={x_times_2m:5d} -> {no_ike_status:15s}, {all_labels_status}")
else:
    print("No such cases found.")

print()
print(f"Total inconsistent: {len(inconsistent_norms)}")
