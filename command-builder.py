f = open("commands.sh", 'w')
for N_norm in range(1, 1200):
    f.write(f"magma -b N_norm:={N_norm} search-12.m\n")
