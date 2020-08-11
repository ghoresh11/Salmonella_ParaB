import argparse
import csv
import os
import networkx as nx


def get_conversion_table(gene_presence_absence_file, first, last):
    ''' read the gene presence absence CSV file in order to be able
    to convert between gene name in the pan-genome and gene id in a specific
    GFF file
    dict ID -> gene name
    genome -> ID of the first gene
    genome -> ID of the last gene
    '''
    print("Reading the gene presence absence file...")
    id_to_gene = {}

    with open(gene_presence_absence_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        for row in reader:
            gene = row["Gene"]
            if gene == first:
                first_per_genome = row  # this row contains the ID of the first gene for all genomes
            elif gene == last:
                last_per_genome = row  # this row contains the ID of the last gene for all genomes

            # either way, save the conversion per ID to the gene name
            for key in row:
                ids = row[key].split()
                for _id in ids:
                    id_to_gene[_id] = gene
    return id_to_gene, first_per_genome, last_per_genome


def add_gff_to_network(G, gff_file, first, last, id_to_gene):
    ''' read the gff file from first gene to last gene, and add them to the panaroo graph'''
    # swap first and last if need be

    first_count = int(first.split("_")[-1])
    last_count = int(last.split("_")[-1])
    if first_count > last_count:  # order is wrong, depends on GFF, flip them
        tmp = first
        first = last
        last = tmp

    # sometimes there's a duplication a gene
    first = first.split("\t")[0]  # if first, take the very first one = 0
    last = last.split("\t")[-1]  # if last, take the very last one = -1

    signature = []
    found_first = False
    found_last = False

    # variables for saving the flanking genes
    prev_gene = "?"
    next_gene = "?"
    prev_contig = "?"

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            if line.startswith("##FASTA"):
                break

            toks = line.strip().split("\t")
            curr_id = toks[-1].split(";")[0].split("=")[-1]

            if found_last:
                if toks[0] != contig:
                    next_gene = "contig break"
                    break
                if curr_id in id_to_gene:  # save the gene next to first
                    next_gene = id_to_gene[curr_id]
                    break
                continue

            if curr_id == first:
                signature.append(id_to_gene[curr_id])
                found_first = True
                contig = toks[0]
                if contig != prev_contig:
                    prev_gene = "contig break"
                continue

            if not found_first:  # don't care about the genes before you found the first one
                if curr_id in id_to_gene:  # save the gene previous to first
                    prev_gene = id_to_gene[curr_id]
                    prev_contig = toks[0]
                continue

            if toks[0] != contig:  # there's a contig break between the
                return(["different contigs"], prev_gene, "NA")
                break

            if curr_id not in id_to_gene:  # some ids don't end up in the roary output
                continue

            signature.append(id_to_gene[curr_id])

            # if we're at the last one, break
            if curr_id == last:
                found_last = True

    # add the edges to the graph, after knowing they're on the same contig
    for i in range(0, len(signature) - 1):
        prev_gene_name = signature[i]
        curr_gene_name = signature[i + 1]
        if G.has_edge(curr_gene_name, prev_gene_name):
            # if there is, add 1 to the weight
            G[curr_gene_name][prev_gene_name]["weight"] += 1
        else:
            G.add_edge(curr_gene_name, prev_gene_name, weight=1)
    return(signature, prev_gene, next_gene)


def run(options):
    id_to_gene, first_per_genome, last_per_genome = get_conversion_table(
        options.gene_presence_absence_file, options.first, options.last)
    gff_files = os.listdir(options.gff_dir)
    G = nx.Graph()
    out = open(options.first + "_" + options.last + "_sigs.csv", "w")
    out.write("Genome,Signature,Previous Gene, Next Gene\n")
    print("Reading GFF files...")
    for gff_file in gff_files:
        genome = gff_file.replace(".gff", "")
        # first and last genes must be present in the genome to be processed.
        if first_per_genome[genome] == "" or last_per_genome[genome] == "":
            out.write(genome + "," + "NA" + "\n")
            continue
        signature, prev_gene, next_gene = add_gff_to_network(G, os.path.join(options.gff_dir, gff_file), first_per_genome[genome], last_per_genome[genome], id_to_gene)
        if signature[0] == options.last:
            signature.reverse()
            tmp = prev_gene
            prev_gene = next_gene
            next_gene = tmp
        out.write(genome + "," + ":".join(signature) + "," + prev_gene + "," + next_gene + "\n")

    out.close()
    # save the network to a file
    nx.write_gml(G, options.first + "_" + options.last + ".gml")

    return


def get_options():
    parser = argparse.ArgumentParser(
        description='Create a panaroo like graph between two genes, from a Roary output')
    # input options
    parser.add_argument('--first', required=False,
                        type=str, default="group_543",
                        help='Name of the first gene')
    parser.add_argument('--last',
                        required=False,
                        type=str, default="trg_3",
                        help='Name of the last gene')
    parser.add_argument('--gene_presence_absence_file', help='csv output from Roary',
                        type=str,
                        default="/lustre/scratch118/infgen/team216/gh11/alyce/gene_presence_absence.csv")  # change when on farm
    parser.add_argument('--gff_dir', help='directory containing all the GFF files',
                        type=str,
                        default="/lustre/scratch118/infgen/team216/gh11/alyce/gffs/")  # change when on farm
    parser.add_argument('--out', help='Names of output files', default="first_second",
                        type=str,
                        required=False)
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
