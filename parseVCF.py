# Import argv form the sys module
from sys import argv

# Read the command line arguments for the vcf file and gene position
vcf_file = argv[1]
start = int(argv[2])
end = int(argv[3])

outfile = open("parsed_vcf_file.vcf", "w")

# Open and read the vcf file
with open(vcf_file) as infile:
	for line in infile:
		# Only interested in the variant lines and not the comment lines at the start
		if not line.startswith("#"):  # Skips the comment lines
			# Split the line into a llist so the postion at index 1 can be checked
			split_line = line.split()
			# Identify the lines that are within the required gene and print them to the output file
			pos = split_line[1]
			if int(pos) >= start and int(pos) <= end:
				outfile.write(line)
outfile.close()	
