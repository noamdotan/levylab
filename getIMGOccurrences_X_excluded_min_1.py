"""this script, given a csv file of potential PTs domains containing the column 'ID' (pfam id of domain's family)
and path to IMG database, outputs csv of all genomes and genes in IMG that contain potential PTS domains
according to a logic defined by the global parameters below
IMPORTANT NOTE: we need to know all domains types
"""
# todo add this functionality?
# IMPORTANT NOTE: domain of type REPEAT
# will be considered in the gene if it appears at least MIN_NUM_REPEATS consecutively

import os
import sys
import pandas as pd

TRAFFICKING = 'F'
REPEAT = 'R'
PRETOXIN = 'P'
TOXIN = 'X'

MIN_REPEATS = 1
MIN_NUM_DOMAINS = 1
MIN_NUM_DOMAINS_TYPES = 1
MIN_NUM_DOMAINS_PFAM = 1
MIN_NUM_DOMAINS_TIGRFAM = 1
MIN_NUM_DOMAINS_TYPES_PFAM = 1
MIN_NUM_DOMAINS_TYPES_TIGRFAM = 1
EXCLUDED_TYPES = set(list([TOXIN]))
ID_COL_NAME = 'ID'
TYPE_COL_NAME = 'Type'
CSV_ENCODEING = 'iso-8859-1'


def check_input():
    """
    :return: tuple of (potential_PTs_domains_csv_path, IMG_path, output_IMG_occurences_path)
    """
    if len(sys.argv) != 4:
        print(
            'Usage: <potential PTs domains csv path> <IMG path> <Output for IMG occurrences path>')
        sys.exit(1)
    if not os.path.exists(sys.argv[1]) or \
            not os.path.exists(sys.argv[2]):
        print('one or more of the paths is invalid')
        sys.exit(1)
    if os.path.exists(sys.argv[3]):
        if input('Be careful! output file already exists, are you sure you want to override? (any key/no) ') == 'no':
            sys.exit(1)

    return sys.argv[1], sys.argv[2], sys.argv[3]


def check_min_num_repeats(gene_domain_ids_ordered, common_ids):
    # getting all domains of type repeat in this gene
    repeats_ids = []
    for domain_id in common_ids:
        if all_domains.loc[all_domains[ID_COL_NAME] == domain_id][TYPE_COL_NAME].array[0] == REPEAT:
            repeats_ids.append(domain_id)
    # if one of repeat types domains is found min num repeats return True
    for domain_id in repeats_ids:
        if any([domain_id for i in range(MIN_REPEATS)] == gene_domain_ids_ordered[i:i + MIN_REPEATS]
               for i in range(len(gene_domain_ids_ordered) - MIN_REPEATS)):
            return True
    # if none were found, return False
    return False


def get_genomes_genes():
    """iterating over all genomes, for all genome iterating on all it's genes and for each gene
    checking if it's acceptable for us by the logic defined by above global variables
    returns data frame with columns GENOMEID and GENEID, CONTAINS with no duplicates"""

    pfam_domains_ids_set = set(all_domains[all_domains[ID_COL_NAME].str.startswith('PF')][ID_COL_NAME])
    tigrfam_domains_ids_set = set(all_domains[all_domains[ID_COL_NAME].str.startswith('TIGR')][ID_COL_NAME])
    genomes_genes = pd.DataFrame(columns=['GENOMEID',
                                          'GENEID',
                                          'CONTAINED_IDS',
                                          'NUM_CONTAINED_TYPES',
                                          'CONTAINED_TYPES'])
    # genomes = os.listdir(IMG_path)
    genomes = list([genome for genome in os.listdir(IMG_path) if os.path.isdir(os.path.join(IMG_path, genome))])
    # iterating over all genomes
    print('Begin itereting on ' + str(len(genomes)) + ' genomes')
    for index, genome in enumerate(genomes):
        print('checking genome ' + str(index) + ': ' + str(genome) + ', out of ' + str(len(genomes)))
        pfam_filename = genome + '.pfam.tab.txt'
        tigrfam_filename = genome + '.tigrfam.tab.txt'
        pfam_path = os.path.join(IMG_path, genome, pfam_filename)
        tigrfam_path = os.path.join(IMG_path, genome, tigrfam_filename)
        if not os.path.exists(pfam_path) or not os.path.exists(tigrfam_path):
            continue
        # finding genes in genome with at least one domain
        genome_genes_with_pfam = find_genes_in_genome_by_pfam(genome, pfam_domains_ids_set, pfam_path)
        genome_genes_with_tigrfam = find_genes_in_genome_by_tigrfam(genome, tigrfam_domains_ids_set, tigrfam_path)

        # merging pfam and tigrfam results
        genome_genes_merged = merge_pfam_tigrfam(genome_genes_with_pfam, genome_genes_with_tigrfam)

        # keeping only genes that have MIN_NUM_DOMAINS_TYPES and no EXCLUDED_TYPES
        # keeping only genes that have MIN_NUM_DOMAINS_TYPES
        genome_genes_merged.drop(genome_genes_merged[genome_genes_merged['NUM_CONTAINED_TYPES'] < MIN_NUM_DOMAINS_TYPES].index, inplace=True)
        # keeping only genes that have no EXCLUDED_TYPES
        indexes_to_drop = list()
        for idx, row in genome_genes_merged.iterrows():
            if len(EXCLUDED_TYPES.intersection(row['CONTAINED_TYPES'])) > 0:
                indexes_to_drop.append(idx)
        genome_genes_merged.drop(indexes_to_drop, inplace=True)
        
        # appending
        genomes_genes = genomes_genes.append(genome_genes_merged)
    # we deal with duplicates allready inside each genome and each genome genes should differ from all other genome genes
    #genomes_genes.drop_duplicates(subset=['GENOMEID',
    #                                      'GENEID',
    #                                      'CONTAINED_IDS',
    #                                      'NUM_CONTAINED_TYPES',
    #                                      'CONTAINED_TYPES'], inplace=True)
    return genomes_genes


def merge_pfam_tigrfam(genome_genes_with_pfam, genome_genes_with_tigrfam):
    genome_genes_merged = pd.DataFrame(columns=['GENOMEID',
                                                'GENEID',
                                                'CONTAINED_IDS',
                                                'NUM_CONTAINED_TYPES',
                                                'CONTAINED_TYPES'])
    for tigrfam_idx, tigrfam_row in genome_genes_with_tigrfam.iterrows():

        for pfam_idx, pfam_row in genome_genes_with_pfam.iterrows():

            if pfam_row['GENOMEID'] == tigrfam_row['GENOMEID'] and \
                    pfam_row['GENEID'] == tigrfam_row['GENEID']:
                merged_ids = pfam_row['CONTAINED_IDS'] | tigrfam_row['CONTAINED_IDS']
                merged_types = pfam_row['CONTAINED_TYPES'] | tigrfam_row['CONTAINED_TYPES']
                genome_genes_merged.loc[len(genome_genes_merged)] = [pfam_row['GENOMEID'],
                                                                     pfam_row['GENEID'],
                                                                     merged_ids,
                                                                     len(merged_types),
                                                                     merged_types]
    genome_genes_merged = genome_genes_merged.append(genome_genes_with_pfam, ignore_index=True)
    genome_genes_merged = genome_genes_merged.append(genome_genes_with_tigrfam, ignore_index=True)
    genome_genes_merged.drop_duplicates(subset=['GENOMEID', 'GENEID'], keep='first', inplace=True)

    return genome_genes_merged


def find_genes_in_genome_by_pfam(genome, pfam_domains_ids_set, pfam_path):
    genome_genes_with_pfam = pd.DataFrame(columns=['GENOMEID',
                                                   'GENEID',
                                                   'CONTAINED_IDS',
                                                   'NUM_CONTAINED_TYPES',
                                                   'CONTAINED_TYPES'])
    with open(pfam_path, 'r') as genes_pfam_annot:

        # iterating on all it's genes and for each gene checking if it's acceptable
        gene_pfam_ids = set()
        gene_pfam_ids_ordered = list()
        prev_gene_IMG_id = None
        genes_pfam_annot.readline()  # skip first headline line
        for row in genes_pfam_annot:
            cur_gene_IMG_id = row.split('\t')[0]

            pfam_id = row.split('\t')[8]
            if len(pfam_id) < 5:
                continue
            PF_id = pfam_id[:2].upper() + pfam_id[4:]

            # if we are looking at the same gene and only (maybe) different pfam domain id, update set
            if cur_gene_IMG_id == prev_gene_IMG_id:
                gene_pfam_ids.add(PF_id)
                gene_pfam_ids_ordered.append(PF_id)
            # else checking if prev_gene is acceptable
            else:
                intercept = pfam_domains_ids_set.intersection(gene_pfam_ids)
                # if contains at least MIN_NUM_DOMAINS num of known domains
                if len(intercept) >= MIN_NUM_DOMAINS_PFAM:

                    types = set(all_domains.loc[all_domains[ID_COL_NAME].isin(intercept)][TYPE_COL_NAME])
                    # checking contains at least MIN_NUM_DOMAINS_TYPES different types
                    if len(types) >= MIN_NUM_DOMAINS_TYPES_PFAM:
                        # checking that no excluded types are in gene
                        #if len(types.intersection(EXCLUDED_TYPES)) == 0:
                            # todo as asaf if to put or not this check, right now not:
                            #  check that if one of necessary types is Repeat than it appears enough times
                            # if len(types) == 2 and REPEAT in types:
                            #    if check_min_num_repeats(gene_pfam_ids_ordered, intercept):
                            #        genomes_genes.loc[len(genomes_genes)] = [genome, prev_gene_IMG_id, str(types)]
                            #        print('gene ' + str(prev_gene_IMG_id) + ' in genome ' + str(genome) + ' contains '
                            #              + str(intercept))
                            # else:
                            # adding location in IMG and types of domains in gene to df
                            genome_genes_with_pfam.loc[len(genome_genes_with_pfam)] = [genome,
                                                                                       prev_gene_IMG_id,
                                                                                       intercept,
                                                                                       len(types),
                                                                                       types]
                            # print('gene ' + str(prev_gene_IMG_id) + ' in genome ' + str(genome) + ' contains '
                            #      + str(intercept))
                #  resetting and updating with current row pfam id
                gene_pfam_ids.clear()
                gene_pfam_ids_ordered.clear()
                gene_pfam_ids.add(PF_id)
                gene_pfam_ids_ordered.append(PF_id)

            prev_gene_IMG_id = cur_gene_IMG_id
    return genome_genes_with_pfam


def find_genes_in_genome_by_tigrfam(genome, tigrfam_domains_ids_set, tigrfam_path):
    genome_genes_with_tigrfam = pd.DataFrame(columns=['GENOMEID',
                                                      'GENEID',
                                                      'CONTAINED_IDS',
                                                      'NUM_CONTAINED_TYPES',
                                                      'CONTAINED_TYPES'])

    with open(tigrfam_path, 'r') as genes_tigrfam_annot:

        # iterating on all it's genes and for each gene checking if it's acceptable
        gene_tigrfam_ids = set()
        gene_tigrfam_ids_ordered = list()
        prev_gene_IMG_id = None
        genes_tigrfam_annot.readline()  # skip first headline line
        for row in genes_tigrfam_annot:
            cur_gene_IMG_id = row.split('\t')[0]

            tigrfam_id = row.split('\t')[6]
            if len(tigrfam_id) < 9:
                continue
            TIGR_id = tigrfam_id[:4].upper() + tigrfam_id[4:]

            # if we are looking at the same gene and only (maybe) different pfam domain id, update set
            if cur_gene_IMG_id == prev_gene_IMG_id:
                gene_tigrfam_ids.add(TIGR_id)
                gene_tigrfam_ids_ordered.append(TIGR_id)
            # else checking if prev_gene is acceptable
            else:
                intercept = tigrfam_domains_ids_set.intersection(gene_tigrfam_ids)
                # if contains at least MIN_NUM_DOMAINS num of known domains
                if len(intercept) >= MIN_NUM_DOMAINS_TIGRFAM:

                    types = set(all_domains.loc[all_domains[ID_COL_NAME].isin(intercept)][TYPE_COL_NAME])
                    # checking contains at least MIN_NUM_DOMAINS_TYPES different types
                    if len(types) >= MIN_NUM_DOMAINS_TYPES_TIGRFAM:
                        # checking that no excluded types are in gene
                        #if len(types.intersection(EXCLUDED_TYPES)) == 0:
                            # todo as asaf if to put or not this check, right now not:
                            #  check that if one of necessary types is Repeat than it appears enough times
                            # if len(types) == 2 and REPEAT in types:
                            #    if check_min_num_repeats(gene_pfam_ids_ordered, intercept):
                            #        genomes_genes.loc[len(genomes_genes)] = [genome, prev_gene_IMG_id, str(types)]
                            #        print('gene ' + str(prev_gene_IMG_id) + ' in genome ' + str(genome) + ' contains '
                            #              + str(intercept))
                            # else:
                            # adding location in IMG and types of domains in gene to df
                            genome_genes_with_tigrfam.loc[len(genome_genes_with_tigrfam)] = [genome,
                                                                                             prev_gene_IMG_id,
                                                                                             intercept,
                                                                                             len(types),
                                                                                             types]
                            # print('gene ' + str(prev_gene_IMG_id) + ' in genome ' + str(genome) + ' contains '
                            #      + str(intercept))
                #  resetting and updating with current row pfam id
                gene_tigrfam_ids.clear()
                gene_tigrfam_ids_ordered.clear()
                gene_tigrfam_ids.add(TIGR_id)
                gene_tigrfam_ids_ordered.append(TIGR_id)

            prev_gene_IMG_id = cur_gene_IMG_id

    return genome_genes_with_tigrfam


if __name__ == '__main__':
    domains_path, IMG_path, out_path = check_input()
    print('domains path: ' + str(domains_path))
    print('IMG path: ' + str(IMG_path))
    print('output path: ' + str(out_path))
    print('Reading domains csv')
    all_domains = pd.read_csv(domains_path, encoding=CSV_ENCODEING)
    print('Searching for acceptable occurrences in IMG')
    output_df = get_genomes_genes()
    output_df.to_csv(out_path, index=False, encoding=CSV_ENCODEING)
