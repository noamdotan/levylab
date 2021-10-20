"""
This script gets as input IMG metadata csv, IMG genomes to work with (all IMG genomes with Bacteria Domain), Ribosomal COGs list,
dir with COGs.faa for each genome from Alex's list, output dir,
and num of representitves for each family.
It creates the MSA needed for a phylo tree to be created after with fasttree2/VeryFastTree for instance.
The MSA created is from 1 genome per uniq family out of genomes to work with, each genome with the COGs chosen and
all gaps if some specific COG doesnt exists in it.
The representative for each family is the one with the highest number of ribosomal COGs from the list of 29(30?) COGs
"""


import pandas as pd
import os
import sys
from collections import defaultdict
import fastaparser
import numpy as np
from sklearn.cluster import KMeans
import seaborn as sbn
import matplotlib.pyplot as plt
import pickle
import functools

PROBLEMATIC_FAMILIES = {'Campylobacteraceae', 'Balneolaceae', 'Thermodesulfobacteriaceae', 'Thermodesulfobiaceae',
                        'Acholeplasmataceae', 'Anaeroplasmataceae', 'Entomoplasmataceae', 'Spiroplasmataceae',
                        'Mycoplasmataceae', 'Nitrospinaceae', 'Nitrospiraceae', 'unclassified'}


def check_input():
    if len(sys.argv) != 7:
        print(
            'Usage: <IMG metadata csv> <IMG genomes to work with(non metagenomes) from alex file> <universal COGS list> <dir with COGS.faa for each IMG genome (no metagenomes)> <output directory> <num representatives from each family>')
        sys.exit(1)

    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]


def create_fastas_from_genomes(res_df, COGS_dir_path, cogs_set):
    # create COGS_to_SEQS dict
    COG_to_SEQS = defaultdict(list)
    for genome in list(res_df['GENOME']):
        if os.path.exists(os.path.join(COGS_dir_path, str(genome)+'_ribosomal_cogs_no_dups_lowestEval.faa')):
            with open(os.path.join(COGS_dir_path, str(genome)+'_ribosomal_cogs_no_dups_lowestEval.faa'), 'r') as genome_faa:
                reader = fastaparser.Reader(genome_faa, parse_method='quick')
                for seq in reader:
                    seq_id = seq.header[1:].split('_')[0]
                    if seq_id in cogs_set:
                        COG_to_SEQS[seq_id].append((genome, seq.sequence))

    os.system('mkdir ' + os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs'))
    # create fasta and msa for each COG
    for cog_id in sorted(list(cogs_set)):
        # create faa
        with open(os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id+'.faa'), 'w') as cog_faa:
            for genome, seq in COG_to_SEQS[cog_id]:
                cog_faa.write('>'+str(genome)+'_'+str(cog_id) + '\n' +
                              seq + '\n')
        # create msa
        clustalo_path = 'C:\\Users\\noamdotan\\Desktop\\Study\\2021\\Project\\tools\\clustal-omega-1.2.2-win64\\clustalo'
        os.system(clustalo_path + ' -i ' + os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id+'.faa') + ' -o ' + os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id+'.msa'))

    # fill in non-existing COGS with gaps alignment
    genomes_set = set(res_df['GENOME'])
    for cog_id in sorted(list(cogs_set)):
        gaps = None
        if os.path.exists(os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id + '.msa')):
            with open(os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id + '.msa'), 'r') as cog_msa:
                reader = fastaparser.Reader(cog_msa)
                seq = next(reader)
                gaps = '-'*len(seq.sequence)
        if os.path.exists(os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id + '.msa')):
            with open(os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id + '.msa'), 'a') as cog_msa:
                cog_genomes = set([tup[0] for tup in COG_to_SEQS[cog_id]])
                missing_genomes = genomes_set.difference(cog_genomes)
                for genome in missing_genomes:
                    cog_msa.write('>'+str(genome)+'_'+str(cog_id) + '\n' +
                                  gaps + '\n')

    # concat for each genome it's COGS to a single MSA file
    os.system('mkdir ' + os.path.join(out_dir_path, 'msa_for_tree' + '_' + str(num_rprs) + '_rprs'))
    genomes_list = list(res_df['GENOME'])
    COG_GENOME_to_SEQ_dict = dict()
    for cog_id in cogs_set:
        with open(os.path.join(out_dir_path, 'COGS_faa_msa' + '_' + str(num_rprs) + '_rprs', cog_id+ '.msa'), 'r') as cog_faa:
            reader = fastaparser.Reader(cog_faa, parse_method='quick')
            for seq in reader:
                genome = seq.header[1:].split('_')[0]
                COG_GENOME_to_SEQ_dict[(genome, cog_id)] = seq.sequence

    with open(os.path.join(os.path.join(out_dir_path, 'msa_for_tree' + '_' + str(num_rprs) + '_rprs', 'msa_for_tree.faa')), 'w') as msa_for_tree_file:
        for genome in genomes_list:
            msa_for_tree_file.write('>' + str(genome) + '\n')
            for cog_id in sorted(list(cogs_set)):
                msa_for_tree_file.write(COG_GENOME_to_SEQ_dict[(genome, cog_id)])
            msa_for_tree_file.write('\n')


def MSA(clustalo_dir_path, seqs_to_align, output_file):
    clustalo = os.path.join(clustalo_dir_path, 'clustalo')
    input_arg = '-i ' + seqs_to_align
    output_arg = '-o ' + output_file
    os.system(clustalo + ' ' + input_arg + ' ' + output_arg)


def get_present_cogs_list_for_genome(genome, fasta_dir_path):
    cogs_list = list()
    fasta_file_path = os.path.join(fasta_dir_path, genome+'_ribosomal_cogs_no_dups_lowestEval.faa')
    if os.path.exists(fasta_file_path):
        with open(fasta_file_path, 'r') as genome_cogs_faa_file:
            parser = fastaparser.Reader(genome_cogs_faa_file, parse_method='quick')
            for seq in parser:
                cog_id = seq.header[1:].split('_')[0]
                cogs_list.append(cog_id)
    return cogs_list


if __name__ == "__main__":
    IMG_metadata_csv_path, genomes_to_work_with_list_path, COGS_list_path, COGS_dir_path, out_dir_path, num_rprs = check_input()
    num_rprs = int(num_rprs)
    # getting relevant metadata only of genome to work with
    with open(genomes_to_work_with_list_path) as genomes_list_file:
        genomes_to_work_with_set = set(line.strip() for line in genomes_list_file)
    with open(COGS_list_path) as COGS_list_file:
        COGS_ordered = list(line.strip() for line in COGS_list_file)

    cogs_set = set(COGS_ordered)

    # get meta to include only genomes to work with
    df_meta = pd.read_csv(IMG_metadata_csv_path)
    df_meta['IMG Genome ID'] = df_meta['IMG Genome ID'].astype(str).str.split(pat='.').str[0]
    df_meta = df_meta[['IMG Genome ID', 'Family', 'Phylum', 'Sequencing Status']]
    df_meta = df_meta[df_meta['IMG Genome ID'].isin(genomes_to_work_with_set)]

    # get meta to include only genomes with classified phylum
    df_meta = df_meta[df_meta['Phylum'] != 'unclassified']

    # ommit problematic families
    df_meta = df_meta[~df_meta['Family'].isin(PROBLEMATIC_FAMILIES)]

    # families in meta
    all_genomes_to_work_with_families_set = set(df_meta['Family'])

    # count how many cogs for each genome and sort
    df_meta['list_present_cogs'] = df_meta['IMG Genome ID'].apply(lambda genome: get_present_cogs_list_for_genome(
                                                                genome, COGS_dir_path))

    df_meta['num_cogs'] = df_meta['list_present_cogs'].str.len()
    # create tree with genomes representatives that have the most COGs out of the family.
    # choose Finished or Permanent draft if have, if not than Draft.
    res_df = pd.DataFrame(columns=['GENOME', 'FAMILY', 'COGS_LIST', 'COGS_NUM'])
    for family, df_fam in df_meta.groupby(by=['Family']):
        df_fam_sorted = df_fam.sort_values(by=['num_cogs'], ascending=False)
        #df_fam_sorted_no_draft = df_fam_sorted[df_fam_sorted['Sequencing Status'] != 'Draft']
        # if there are no Draft, take one with top num cogs
        #if(len(df_fam_sorted_no_draft) > 0):
        #    for i in range(min(num_rprs, len(df_fam_sorted_no_draft))):
        #        res_df.loc[len(res_df),] = [df_fam_sorted_no_draft.iloc[i]['IMG Genome ID'],
        #                                    df_fam_sorted_no_draft.iloc[i]['Family'],
        #                                    df_fam_sorted_no_draft.iloc[i]['list_present_cogs'],
        #                                    df_fam_sorted_no_draft.iloc[i]['num_cogs']]
        ## else, take top num_rprs or whatever we have
        #else:
        for i in range(min(num_rprs, len(df_fam_sorted))):
            res_df.loc[len(res_df), ] = [df_fam_sorted.iloc[i]['IMG Genome ID'],
                                        df_fam_sorted.iloc[i]['Family'],
                                        df_fam_sorted.iloc[i]['list_present_cogs'],
                                         df_fam_sorted.iloc[i]['num_cogs']]
    # output chosen representatives to be used in the R script

    res_df = res_df.sort_values(by=['COGS_NUM'], ascending=False)
    res_df.to_csv('C:\\Users\\noamdotan\\Desktop\\Study\\2021\\Project\\Toxins_Paper\\phylo_tree\\ribosomal_COGs_GENOMES_for_phylo_tree_plot_'+str(num_rprs)+'_rprs_each_family_WO_PROBLEMATIC_FAMILIES_WO_unclassified.txt',
                  columns=['GENOME'], index=False)
    create_fastas_from_genomes(res_df, COGS_dir_path, cogs_set)