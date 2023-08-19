# Read the command line arguments for the vcf file and gene position
from sys import argv

vcf_file = argv[1]
start = int(argv[2])
end = int(argv[3])

# Open and read the vcf file
with open(vcf_file) as file_in:
	for line in file_in: 
		if not line.startswith("#"): #split the line into a list, variant position is index 1
			split_line = line.split()
			var_pos = split_line[1]
			#lines within the required gene
			if int(var_pos) >= start and int(var_pos) <= end: 
				#splitting annotation description on the pipe "|"
				ann_list = split_line[7].split("|")
				#using a for loop, stem of 15 as each transcript annotation has 15 elements in the list
				for i in range(0, len(ann_list)-1, 15): #final value ignored to not produce out of bounds error
					variant_type = ann_list[i+1]
					#is a missense, extract other relevant info
					if variant_type == "missense_variant":
						transcript = ann_list[i+6]
						nuc_change = ann_list[i+9]
						aa_change = ann_list[i+10]
						print("transcript"+ transcript)
						print("Chrom Pos:", var_pos, "Variant effect:", variant_type, "Transcript ID:", transcript, "Nuc change:", nuc_change, "AA change:", aa_change)

