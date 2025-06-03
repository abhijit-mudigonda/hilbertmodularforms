import math

f = open("commands.sh", 'w')

max_norm = 50
chunk_size = 10
cores = 1

for chunk_idx in range(math.floor(max_norm / chunk_size)):
    chunk_min_norm = chunk_idx * chunk_size
    chunk_max_norm = chunk_min_norm + chunk_size - 1
    for core in range(cores):
        f.write(f"magma -b CORES:={cores} CORE:={core} MIN_NORM:={chunk_min_norm} MAX_NORM:={chunk_max_norm} single-target.m\n")
