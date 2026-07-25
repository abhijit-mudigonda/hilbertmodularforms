infile = open("ideal-labels-400.txt", "r")
outfile = open("precompute-commands.sh", "w")

for line in infile:
    label = line.strip()
    if label:
        outfile.write(f"magma -b N_label:={label} precompute-weight2-single.m\n")
