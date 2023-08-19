# Import argv form the sys module
from sys import argv

# Read the command line arguments for the vcf file and gene position
vcf_file = argv[1]
start = int(argv[2])
end = int(argv[3])

# Open and read the vcf file
with open(vcf_file) as infile:
        for line in infile:
                # Only interested in the variant lines and not the comment lines at the start
                if not line.startswith("#"):  # Skips the comment lines
                        # Split the line into a list and the variant position is index 1
                        split_line = line.split()
                        pos = split_line[1]
                        # Identify the lines that are within the required gene
                        if int(pos) >= start and int(pos) <= end:
                                # Split the annotation description on the pipe '|'
                                ann_list = split_line[7].split("|")
                                # Using a for loop use a step of 15 as each trasnscript annotation is 15 elements in the list
                                for i in range(0, len(ann_list)-1, 15): # Length is ann_list - 1 so final value is ignored as it is not part of the annotation and would otherwise cause an index out of bounds error
                                        variant_type = ann_list[i+1]
                                        # If it is a missense variant extract the other relevant information and print it

                                        if variant_type == "splice_region_variant&intron_variant":
                                                transcript = ann_list[i+6]
                                                nuc_change = ann_list[i+9]
                                                aa_change = ann_list[i+10]
                                                print('transcript'+ transcript)
                                                print("Chrom Pos:", pos, "Variant effect:", variant_type, "Transcript ID:", transcript, "Nuc change:", nuc_change, "AA change:", aa_change)


