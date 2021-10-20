"""
this script gets the folder which conatins all fasta, msa and hmm of these clusters, and outputs to
output directory only faa msa hmm files which are considered taxonomically diversed enough according
to below defined parameters. in addition teh script outputs csv file of taxon both in in_dir and out_dir
"""

import os
import sys
import pandas as pd
from shutil import copyfile
import fastaparser

MIN_DIFFERENT_PHYLUM = 3
#MIN_DIFFERENT_CLASS = 2
#MIN_DIFFERENT_ORDER = 3
#MIN_DIFFERENT_FAMILY = 4
MIN_DIFFERENT_GENUS = 3





def check_input():
    if len(sys.argv) != 4:
        print(
            'Usage: <IMG metadata csv> <directory with all faa msa hmm of clusters> <output directory>')
        sys.exit(1)

    return sys.argv[1], sys.argv[2], sys.argv[3]


def getGenomesDF():
    df = pd.DataFrame(columns=['cluster', 'genomes',  'num_phylum', 'num_class', 'num_order', 'num_family', 'num_genus'])
    cluster = None
    for idx, file in enumerate(os.listdir(in_dir_path)):
        print(idx, end='\r', flush=True)
        if file.endswith('.msa.faa'):
            cluster = file[:file.index('.')]
            with open(os.path.join(in_dir_path, file), 'r') as fasta_file:
                reader = fastaparser.Reader(fasta_file, parse_method='quick')
                genomes = set()
                for seq in reader:
                    genome = seq.header[1:].split('_')[0]
                    genomes.add(genome)
                df.loc[len(df)] = [cluster, genomes, 0, 0, 0, 0, 0]

    return df


def getDiverseClusters():
    clusters_genomes_df = getGenomesDF()
    metadata_df = pd.read_csv(metadata_file_path, dtype=object)
    phyla = set()
    classes = set()
    orders = set()
    families = set()
    genuses = set()
    print(len(clusters_genomes_df))
    for idx1, row1 in clusters_genomes_df.iterrows():
        print(idx1, end='\r', flush=True)
        tmp = metadata_df[metadata_df['IMG Genome ID'].isin(row1['genomes'])]
        for idx2, row2 in tmp.iterrows():
            phyla.add(row2['Phylum'])
            classes.add(row2['Class'])
            orders.add(row2['Order'])
            families.add(row2['Family'])
            genuses.add(row2['Genus'])

        for tax_set in [phyla, classes, orders, families, genuses]:
            if 'Unclassified' in tax_set:
                tax_set.remove('Unclassified')
            if 'unclassified' in tax_set:
                tax_set.remove('unclassified')

        row1['num_phylum'] = len(phyla)
        row1['num_class'] = len(classes)
        row1['num_order'] = len(orders)
        row1['num_family'] = len(families)
        row1['num_genus'] = len(genuses)
        phyla.clear()
        classes.clear()
        orders.clear()
        families.clear()
        genuses.clear()

    clusters_genomes_df.to_csv(os.path.join(in_dir_path, 'taxon_counts.csv'))
    clusters_genomes_df = clusters_genomes_df[clusters_genomes_df['num_genus'] >= MIN_DIFFERENT_GENUS]
    clusters_genomes_df.reset_index(drop=True, inplace=True)
    clusters_genomes_df.to_csv(os.path.join(out_dir_path, 'taxon_counts.csv'))
    return clusters_genomes_df


if __name__ == '__main__':
    metadata_file_path, in_dir_path, out_dir_path = check_input()
    df = getDiverseClusters()
    for index, row in df.iterrows():
        for postfix in ['.faa', '.msa.faa', '.msa.faa.hmm']:
            copy_this = os.path.join(in_dir_path, row['cluster'] + postfix)
            create_this = os.path.join(out_dir_path, row['cluster'] + postfix)
            if os.path.exists(copy_this):
                copyfile(copy_this, create_this)
