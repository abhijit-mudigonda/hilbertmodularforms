outfile = open("label-commands.sh", 'w')
infile = open("12-labels.txt", 'r')

for line in infile.readlines():
    chi_label = line.strip()
    outfile.write(f"magma -b chi_label:={chi_label} label-search-12.m\n")



