outfile = open("label-commands.sh", 'w')
infile = open("cm-12-labels-5.txt", 'r')

for line in infile.readlines():
    chi_label = line.strip()
    outfile.write(f"magma -b chi_label:={chi_label} tensor-induction-eigs.m\n")



