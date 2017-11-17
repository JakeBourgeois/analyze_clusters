#! usr/bin/python

""" analyze_clusters is a python file written for Oggy's inversion clustering project. As input, it takes in a list
# of NCBI GenBank accession numbers and corresponding clustering data from Oggy's bash and R script analysis. It
# then downloads a flat gb file from entrez, parses that file to extract CDS information, and then lines up each
# cluster position with the nearby genes.
"""

import time

# I keep my classes and junk in cluster_tools.py
try:
    from cluster_tools import *
except ImportError:
    print("Clustering tools file not found. Please ensure cluster_tools.py is in the proper directory.")
    quit()


def align_clusters_to_genes(max_genes):

    # Define our folder and file names holding the data
    entrez_folder_name = 'Entrez Data'
    gene_folder_name = 'Gene Data'
    results_folder_name = 'Cluster Data'
    input_filename = 'accession_list.txt'
    translations_filename = '__cluster_gene_translations_fasta.txt'
    params_filename = '__result_parameters.txt'


    # Define our folder paths
    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    working_path = os.path.join(desktop_path, "Cluster Detection")
    results_path = os.path.join(working_path, results_folder_name)


    entrez_path = os.path.join(working_path, entrez_folder_name)
    gene_path = os.path.join(working_path, gene_folder_name)




    # If our folder paths do not exist, make them.
    paths = [desktop_path, entrez_path, gene_path, results_path]
    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)
    del paths

    # Define our input files
    accessions_input = os.path.join(working_path, input_filename)
    translations_output = os.path.join(results_path, translations_filename)

    # If the translation file exists, delete it
    if os.path.exists(translations_output):
        os.remove(translations_output)

    # Define running variables
    ntol = 100000  # number of nucleotides to look at around given cluster position
    max_genes = 6  # maximum number of genes to assign to a cluster
    accessions_list = list()  # holds our input accession numbers
    my_bugs = list()  # holds our bugs

    # Load the accession numbers
    with open(accessions_input, 'r') as f:
        for line in f:
            accessions_list.append(line.split('\n')[0])

    # For each accession number...
    for acc_num in accessions_list:

        print("Analyzing clusters for accession number", acc_num)

        # Define the file names
        acc_path = os.path.join(results_path, acc_num)
        entrez_file = os.path.join(entrez_path, acc_num+'.txt')
        gene_file = os.path.join(gene_path, acc_num+'.csv')
        cluster_file = os.path.join(acc_path, acc_num+'.csv')
        results_file = os.path.join(acc_path, acc_num+'.tsv')
        graph_path = os.path.join(acc_path, 'Gene Diagrams')

        if not os.path.exists(graph_path):
            os.makedirs(graph_path)

        # Get Entrez Data, if necessary
        get_entrez_data(acc_num, entrez_file)

        # Generate the gene list, if necessary
        find_genes(acc_num, entrez_file, gene_file)

        # Load the genes from the gene list onto a bug class
        my_bug = Bug(accession_num=acc_num)
        my_bug.load_genes_from_file(gene_file)

        # Scan for clusters
        match_clusters_to_genes(my_bug, cluster_file, results_file, translations_output, graph_path, ntol, max_genes)

    # Create the parameters file
    params_file = os.path.join(results_path, params_filename)
    with open(params_file, 'w') as f:
        f.write("Accession list: " + accessions_input + '\n')
        f.write("Nucleotide tolerance: " + str(ntol) + '\n')
        f.write("Maximum nearby genes: " + str(max_genes) + '\n')

    x = time.clock()
    print("Done! Operation took", x, "seconds.")
