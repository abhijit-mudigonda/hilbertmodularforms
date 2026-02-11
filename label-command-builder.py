outfile = open("no-ike-commands.sh", 'w')
infile = open("12-all-labels-2.txt", 'r')

for line in infile.readlines():
    chi_label = line.strip()
    outfile.write(f"magma -b chi_label:={chi_label} blah16.m\n")

