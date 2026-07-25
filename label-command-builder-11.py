infile = open("11-labels-1620-4000.txt", "r")
outfile = open("commands-11.sh", "w")

for line in infile:
    chi_label = line.strip()
    if chi_label:
        outfile.write(f"magma -b chi_label:={chi_label} label-search-11.m\n")
