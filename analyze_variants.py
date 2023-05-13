import json

original_file = "output_VCF_2.json"
output_file = "output_VCF_2_corrected.json"

# In my JSON file, I have lines that produces errors when I try to read them as JSON, so I will remove them
# More specifically, these are the lines that start with a number 0-9


def remove_errors(file, output_file):
    with open(file, "r") as input_file:
        lines = input_file.readlines()
        lines = [line.strip() for line in lines if not line.strip(
        ).startswith(tuple(map(str, range(10))))]

    with open(output_file, "w") as corrected_file:
        corrected_file.writelines(lines)


# Read content from JSON file as python dictionary
def read_json(file):
    with open(file, "r") as file:
        return json.load(file)

# Get list from 'positions' key in dictionary and print number of elements


def get_positions(dictionary):
    print("Number of variants: ", len(dictionary["positions"]))
    return dictionary["positions"]

# Get list from 'genes' key in dictionary


def get_genes(dictionary):
    return dictionary["genes"]

# Check how many variants are in clinvar and cardiodb and how many are in both


def count_variants_in_db(dictionary):
    cnt_clinvar = 0
    cnt_cardiodb = 0
    cnt_both = 0
    for position in dictionary:
        is_in_clinvar = False
        is_in_cardiodb = False
        is_in_both = False
        for variant in position["variants"]:
            if "clinvar" in variant:
                is_in_clinvar = True
            if "cardiodb" in variant:
                is_in_cardiodb = True
            if is_in_clinvar and is_in_cardiodb:
                is_in_both = True
        if is_in_clinvar:
            cnt_clinvar += 1
        if is_in_cardiodb:
            cnt_cardiodb += 1
        if is_in_both:
            cnt_both += 1
    print('Number of clinvar appearances: ', cnt_clinvar)
    print('Number of cardiodb appearances: ', cnt_cardiodb)
    print('Number of appearances in both clinvar and cardiodb: ', cnt_both)

# Check how many variants are classified as pathogenic (i.e. significance = "Pathogenic") and show them


def get_pathogenic_variants(dictionary):
    cnt = 0
    pathogenic_variants = []
    for position in dictionary:
        is_pathogenic = False
        for variant in position["variants"]:
            if "clinvar" in variant:
                for clinvar_variant in variant["clinvar"]:
                    if "significance" in clinvar_variant:
                        for significance in clinvar_variant["significance"]:
                            if significance == "pathogenic":
                                is_pathogenic = True
        if is_pathogenic:
            # Add to list
            pathogenic_variants.append(position)
            cnt += 1
    print('Number of pathogenic variants: ', cnt)
    return pathogenic_variants

# Find in which genes are the pathogenic variants:
# For each transcript in list at key transcripts for each variant in list at key variants of position in pathogenic_variants,
# get the value of hgnc key which corresponds to a gene and append it to a list where the key is chromosome and position of pathogenic variant


def get_pathogenic_genes(pathogenic_variants, genes):
    pathogenic_genes = {}
    for pathogenic_variant in pathogenic_variants:
        chromosome = pathogenic_variant["chromosome"]
        position = pathogenic_variant["position"]
        for variant in pathogenic_variant["variants"]:
            if "transcripts" in variant:
                for transcript in variant["transcripts"]:
                    if "hgnc" in transcript:
                        gene = transcript["hgnc"]
                        if chromosome not in pathogenic_genes:
                            pathogenic_genes[chromosome] = {}

                        if position not in pathogenic_genes[chromosome]:
                            pathogenic_genes[chromosome][position] = []

                        if gene not in pathogenic_genes[chromosome][position]:
                            pathogenic_genes[chromosome][position].append(gene)
    # Now, for each chromosome and position, keep only the genes that are found in the genes dictionary, at key name of each element of genes list at key genes
    genes_names = [gene["name"] for gene in genes]
    for chromosome in pathogenic_genes:
        for position in pathogenic_genes[chromosome]:
            # remove elements of pathogenic_genes[chromosome][position] that do not match any element in genes list at key name
            pathogenic_genes[chromosome][position] = [
                gene for gene in pathogenic_genes[chromosome][position] if gene in genes_names]
    print("Pathogenic genes: ", pathogenic_genes)

# I have a CardioDB.tsv file in this same folder that contains variants (chromosome and position) that are associated with cardiovascular diseases
# I want to manually check if any of the variants in my CardioDB file are also in my JSON file


def manual_check(dictionary):
    # skip the lines that begin with #
    cnt = 0
    with open("CardioDB.tsv", "r") as file:
        lines = file.readlines()
        for line in lines:
            if line[0] == "#":
                continue
            else:
                line_split = line.split("\t")
                chromosome = line_split[0]
                position = int(line_split[1])
                for variant in dictionary:
                    if chromosome == variant["chromosome"] and position == variant["position"]:
                        print("CardioDB variant found in JSON file: ",
                              chromosome, position)
                        cnt += 1
    print("Number of CardioDB variants found in JSON file: ", cnt)


def main():
    # remove_errors(original_file, output_file)
    dictionary = read_json(output_file)
    positions = get_positions(dictionary)
    genes = get_genes(dictionary)
    count_variants_in_db(positions)
    pathogenic_variants = get_pathogenic_variants(positions)
    get_pathogenic_genes(pathogenic_variants, genes)
    manual_check(positions)


if __name__ == "__main__":
    main()
