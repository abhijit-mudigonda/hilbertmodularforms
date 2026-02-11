#!/usr/bin/env python3
"""
Find level norms x in no-ike-5.txt such that for every m,
both x*5^m and x*4^m are either in both no-ike-5.txt and 12-all-labels-5.txt or in neither.
"""

import re

def extract_level_norm(line):
    """Extract the level norm (integer part before decimal) from a label."""
    # Pattern: -5.0.1_{x}.{y}_...
    match = re.search(r'-5\.0\.1_(\d+)\.\d+_', line)
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
no_ike_norms = read_level_norms('/home/abhijitm/hmf/no-ike-5.txt')
all_labels_norms = read_level_norms('/home/abhijitm/hmf/12-all-labels-5.txt')

print(f"Level norms in no-ike-5.txt: {len(no_ike_norms)}")
print(f"Level norms in 12-all-labels-5.txt: {len(all_labels_norms)}")
print()

# For each x in no-ike, check if x*5^m and x*4^m are consistent across both files
consistent_norms = []
inconsistent_norms = []

for x in sorted(no_ike_norms):
    # Check all powers of 5 and 4 times x up to a reasonable limit
    is_consistent = True
    violations = []
    
    # Check powers of 5
    for m in range(1, 10):  # 5^9 is already quite large
        x_times_5m = x * (5 ** m)
        
        # Check if it's in both or in neither
        in_no_ike = x_times_5m in no_ike_norms
        in_all_labels = x_times_5m in all_labels_norms
        
        # Consistent means: both True or both False
        if in_no_ike != in_all_labels:
            is_consistent = False
            violations.append(('5', m, x_times_5m, in_no_ike, in_all_labels))
    
    # Check powers of 4
    for m in range(1, 12):  # 4^11 is similar magnitude to 5^9
        x_times_4m = x * (4 ** m)
        
        # Check if it's in both or in neither
        in_no_ike = x_times_4m in no_ike_norms
        in_all_labels = x_times_4m in all_labels_norms
        
        # Consistent means: both True or both False
        if in_no_ike != in_all_labels:
            is_consistent = False
            violations.append(('4', m, x_times_4m, in_no_ike, in_all_labels))
    
    if is_consistent:
        consistent_norms.append(x)
    else:
        inconsistent_norms.append((x, violations))

print("Level norms x in no-ike-5.txt where x*5^m and x*4^m are CONSISTENT (in both or neither):")
print("=" * 70)

if consistent_norms:
    for x in consistent_norms:
        print(f"x = {x}")
else:
    print("No such cases found.")

print()
print(f"Total consistent: {len(consistent_norms)}")
print()

print("Level norms x in no-ike-5.txt where x*5^m or x*4^m are INCONSISTENT:")
print("=" * 70)

if inconsistent_norms:
    for x, violations in inconsistent_norms:
        print(f"\nx = {x}:")
        for base, m, x_times_bm, in_no_ike, in_all_labels in violations:
            no_ike_status = "IN no-ike" if in_no_ike else "NOT in no-ike"
            all_labels_status = "IN all-labels" if in_all_labels else "NOT in all-labels"
            print(f"  base={base}, m={m}: x*{base}^{m}={x_times_bm:6d} -> {no_ike_status:15s}, {all_labels_status}")
else:
    print("No such cases found.")

print()
print(f"Total inconsistent: {len(inconsistent_norms)}")
