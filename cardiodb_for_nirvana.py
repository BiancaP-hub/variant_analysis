import csv

# I have a CardioDB CSV file with 9 columns : Gene, Nucleotide.Change, Protein.Change, Consequence, Ogml.class, Lmm.Class, Phenotype, Type and Location.GRCh37
# I need to create a new file .tsv in a format that is accepted by the tool "Nirvana" (https://github.com/Illumina/Nirvana.git)
# The new file is created using variant information from the CardioDB CSV file

def create_tsv_file():
    with open("cardioDBwithREF.csv", "r") as csv_file:
        with open("CardioDB.tsv", "w") as tsv_file:
            # Write the header of the new file
            write_header(tsv_file)
            # Write the data of the new file
            write_data(tsv_file, csv_file)

    # Eliminate duplicates
    eliminate_duplicates()

    # Sort the file by chromosome and position
    sort_file()
    
    print("File created successfully")

# Create the header of the new file, considering it has 7 rows and 7 columns
# use csv.writer to write the data, and use tab to separate the columns
def write_header(file):
    csv_writer = csv.writer(file, delimiter='\t', lineterminator='\n')
    # First row, first column is #title=CardioDB, the rest of the columns are empty
    csv_writer.writerow(["#title=CardioDB", "", "", "", "", "", ""])
    # Second row, first column is #assembly=GRCh37, the rest of the columns are empty
    csv_writer.writerow(["#assembly=GRCh37", "", "", "", "", "", ""])
    # Third row, first column is #matchVariantsBy=allele, the rest of the columns are empty
    csv_writer.writerow(["#matchVariantsBy=allele", "", "", "", "", "", ""])
    # Fourth row, columns are #CHROM, POS, REF, ALT, OmglClass, LmmClass, Phenotype
    csv_writer.writerow(["#CHROM", "POS", "REF", "ALT", "OmglClass", "LmmClass", "Phenotype"])
    # Fifth row, first colum is #categories, the rest of the columns are filled with "."
    csv_writer.writerow(["#categories", ".", ".", ".", ".", ".", "."])
    # Sixth row, first column is #descriptions, the rest of the columns are filled with "."
    csv_writer.writerow(["#descriptions", ".", ".", ".", ".", ".", "."])
    # Seventh row, first column is #type, columns 5,6,7 are filled with "string" and the rest of the columns are filled with "."
    csv_writer.writerow(["#type", ".", ".", ".", "string", "string", "string"])

# Write the data of the new file
# For each column in the new file (7 columns), get the data from the CardioDB CSV file
# use csv.writer to write the data, and use tab to separate the columns
# Here are the rules for each column:
# #CHROM: Get from LOCATION.GRCH37, POS: Get from LOCATION.GRCH37, REF: Equivalent to CORRECT_REF, ALT: Get from NUCLEOTIDE.CHANGE 
# OmglClass: Equivalent to OGML.CLASS, LmmClass: Equivalent to LMM.CLASS, Phenotype: Equivalent to PHENOTYPE
def write_data(file, csv_file):
    # Read CSV file using csv.reader, first row is header, delimiter is "," and consider quotes around fields
    reader = csv.DictReader(csv_file, delimiter=',')
    writer = csv.writer(file, delimiter='\t', lineterminator='\n')
    # For each row in the CardioDB CSV file, write the data in the new file
    for row in reader:
        # If the type of the variant is not "substitution", skip it
        if row["Type"] != "substitution":
            continue
        # Get the data from column LOCATION.GRCH37.: the characters before the ":" are the chromosome
        chrom = row["Location.GRCh37."].split(":")[0]
        # Get the data from LOCATION.GRCH37: the characters after the ":" are the position
        pos = row["Location.GRCh37."].split(":")[1]
        
        if ">" in row["Nucleotide.Change"]:
            # Get the data from correct_ref
            ref = row["correct_ref"]
            # Get the data from NUCLTEOTIDE.CHANGE: the character after the ">" is the alternative nucleotide
            alt = row["Nucleotide.Change"].split(">")[1][0]
        else:
            continue

        # Get the data from OGML.CLASS
        omgl_class = row["OMGL.class"]
        # Get the data from LMM.CLASS
        lmm_class = row["LMM.class"]
        # Get the data from PHENOTYPE
        phenotype = row["Phenotype"]
        # if either of omgl_class, lmm_class or phenotype are empty, replace them with "."
        if omgl_class == "":
            omgl_class = "."
        if lmm_class == "":
            lmm_class = "."
        if phenotype == "":
            phenotype = "."
        # Write the data in the new file, but do not skip a line between rows
        writer.writerow([chrom, pos, ref, alt, omgl_class, lmm_class, phenotype])

def sort_file():
    # Read the file and sort it by chromosome and position
    with open("CardioDB.tsv", "r") as file:
        # Read the file
        lines = file.readlines()
        # Sort lines that do not start with "#"
        lines[7:] = sorted(lines[7:], key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    # Write the sorted file
    with open("CardioDB.tsv", "w") as file:
        file.writelines(lines)

# Eliminate duplicates by checking if the variant is already in the file (same chromosome, position, reference and alternative nucleotides)
def eliminate_duplicates():
    # Read the file
    with open("CardioDB.tsv", "r") as file:
        # Read the file
        lines = file.readlines()
        # Create a list of variants
        variants = []
        # For each line in the file, check if the variant is already in the list of variants
        for line in lines:
            # If the line starts with "#", skip it
            if line.startswith("#"):
                continue
            
            # Get the chromosome, position, reference and alternative nucleotides
            chrom = line.split("\t")[0]
            pos = line.split("\t")[1]
            ref = line.split("\t")[2]
            alt = line.split("\t")[3]

            # If the variant is not in the list, add it
            if (chrom, pos, ref, alt) not in variants:
                variants.append((chrom, pos, ref, alt))
            # If the variant is already in the list, remove it
            else:
                lines.remove(line)
    
    # Write the file without duplicates
    with open("CardioDB.tsv", "w") as file:
        file.writelines(lines)

if __name__ == "__main__":
    create_tsv_file()

