""""
this script gets as input:
- path to hmmscan output file with pfam annotations of proteins which the cluster is Terminal(N if we look for
    trafficking and C if we look for toxin)
- initial list of known PTs domain which asaf gave me
- path to (non-exiting)csv with modularity score for each query with genome_gene from IMG as id (?) or maybe cluster# as id (?)
"""

import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np
from plotnine import *
from itertools import groupby

SCORE_COLUMNS = ['cluster',
                 'numOfFusedDomains', 'fusedDomains',
                 'numOfKnownPTsFusedDomains', 'knownPTsFusedDomains',
                 'PTStructureScore', 'architectures', 'max_pt_score', 'max_pt_score_architecture', 'avg_arcis_pt_score',
                 'numOfUniqueKnownPTsArchitectures', 'uniqueKnownPTsArchitectures', 'uniqueKnownPTsArchitectures_DOMAINS_NAMES']

TRAFFICKING = 'F'
REPEAT = 'R'
PRETOXIN = 'P'
TOXIN = 'X'
TERMINAL_LENGTH = 133
N_TERMINAL_ARG = 'N'
C_TERMINAL_ARG = 'C'
N_FROM = 30
C_TO = 30
EVAL_THRES = 0.001
COVER_PERCENTAGE_THRES = 0.8
ACC_THRES = 0.85

DOM_NAME_IDX = 0
DOM_ID_IDX = 1
DOM_LENGTH = 2
CLUSTER_GENOME_GENE_PROTEIN_IDX = 3
PROTEIN_LENGTH_IDX = 5
DOM_E_VALUE_IDX = 12
HMM_FROM_IDX = 15
HMM_TO_IDX = 16
PROTEIN_FROM_IDX = 17
PROTEIN_TO_IDX = 18
ACC_IDX = 21

COLOR_MAP = {'F': 'black', 'R': 'green', 'P': 'blue', 'X': 'red', 'PFXXXXX': 'white'}


def check_input():
    if len(sys.argv) != 5:
        print(
            'Usage: <path to hmmscan annotations file> < known PTs domains> <N or C Terminal to score> <Output file path>')
        sys.exit(1)
    if not os.path.exists(sys.argv[1]) or not os.path.exists(sys.argv[2]):
        print('one or more of the paths is invalid')
        sys.exit(1)
    if sys.argv[3] != N_TERMINAL_ARG and sys.argv[3] != C_TERMINAL_ARG:
        print('Terminal arg has to be N or C')
        sys.exit(1)
    if os.path.exists(sys.argv[4]):
        if input('Be careful! output file already exists, are you sure you want to override? (any key/no) ') == 'no':
            sys.exit(1)
    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]


def domain_align_is_ok(evalue, dom_len, hmm_from, hmm_to, acc, protein_from, protein_to, protein_len):
    if evalue <= EVAL_THRES and \
            ((hmm_to - hmm_from) / dom_len) >= COVER_PERCENTAGE_THRES and \
            acc >= ACC_THRES:
        if (terminal == N_TERMINAL_ARG and protein_from >= TERMINAL_LENGTH - N_FROM) or \
                (terminal == C_TERMINAL_ARG and protein_to <= protein_len - TERMINAL_LENGTH + C_TO):
            return True
    return False


def insert_domain_to_architecture(cluster, protein, domain, architectures_dict):
    if cluster not in architectures_dict.keys():
        # architectures_dict[cluster].append(domain)
        architectures_dict[cluster].append((protein, list()))
        architectures_dict[cluster][0][1].append(domain)
    else:
        for idx, arci in enumerate(architectures_dict[cluster]):
            if arci[0] == protein:
                for i, dom in enumerate(arci[1]):
                    if dom[1] > domain[2]:
                        arci[1].insert(i, domain)
                        return
                arci[1].append(domain)  # if it's downstream from all other domains
                return
        # if it's arci of a new protein we add arci of this protein and add this domain to the arci
        architectures_dict[cluster].append((protein, list()))
        architectures_dict[cluster][-1][1].append(domain)


def make_architectures():
    architectures_dict = defaultdict(list)  # key is cluster, value is list of architectures of domains

    with open(annotation_path, 'r') as annotation_file:
        for idx, line in enumerate(annotation_file):
            if line.startswith('#'):
                continue
            print('cur line: ' + str(idx))
            split = line.split()
            dom_name = split[DOM_NAME_IDX]
            dom_id = split[DOM_ID_IDX]
            if '.' in dom_id:
                dom_id = dom_id[:dom_id.index('.')]  # removing postfixes
            dom_length = int(split[DOM_LENGTH])
            cluster_genome_gene_protein = split[CLUSTER_GENOME_GENE_PROTEIN_IDX]
            cluster = cluster_genome_gene_protein.split('~')[0]
            genome_gene = cluster_genome_gene_protein.split('~')[1]
            if len(cluster_genome_gene_protein.split('~')) > 2:
                protein = cluster_genome_gene_protein.split('~')[2]
            else:
                protein = cluster_genome_gene_protein.split('~')[1]
            protein_length = int(split[PROTEIN_LENGTH_IDX])
            evalue = float(split[DOM_E_VALUE_IDX])
            hmm_from = int(split[HMM_FROM_IDX])
            hmm_to = int(split[HMM_TO_IDX])
            protein_from = int(split[PROTEIN_FROM_IDX])
            protein_to = int(split[PROTEIN_TO_IDX])
            acc = float(split[ACC_IDX])
            if domain_align_is_ok(evalue, dom_length, hmm_from, hmm_to, acc, protein_from, protein_to, protein_length):
                domain = (dom_id, protein_from, protein_to)
                insert_domain_to_architecture(cluster, protein, domain, architectures_dict)

    return architectures_dict


def score_modularity():
    known_pts_df = pd.read_csv(known_pt_csv_path, encoding='ISO-8859-1')
    known_pts_set = set(list(known_pts_df['ID']))
    df = pd.DataFrame(columns=SCORE_COLUMNS)
    for cluster in architectures.keys():
        fused_to_domains = set()
        known_pts_fused_to_domains = set()
        knownPTsArcis = set()
        knownPTsArcis_names = set()
        cluster_architectures = list()
        architectures_pt_scores = list()
        for (protein, arci) in architectures[cluster]:
            if len(arci) > 0:
                idx_by_terminal = 0 if terminal == 'N' else -1  # todo check if not overlapping cluster domain
                fused_to_domains.add(arci[idx_by_terminal][0])
                architectures_pt_scores.append(getPTStructureScore(arci, known_pts_df))
                cluster_architectures.append(arci)
                knownPTsArci_arr = [dom[0] for dom in arci if dom[0] in known_pts_set]
                if len(knownPTsArci_arr) > 0:
                    knownPTsArci_arr = [x[0] for x in groupby(knownPTsArci_arr)]  # remove consecutive duplicates
                    knownPTsArci_string = ''.join([(dom + '~') for dom in knownPTsArci_arr])
                    knownPTsArcis.add(knownPTsArci_string)
                    # getting domains names
                    knownPTsArci_arr_names = [known_pts_df.loc[known_pts_df.ID == dom_id, 'Name'].values[0] \
                                              for dom_id in knownPTsArci_arr]
                    knownPTsArci_names_string = ''.join([(dom + '~') for dom in knownPTsArci_arr_names])
                    knownPTsArcis_names.add(knownPTsArci_names_string)

        max_pt_score = max(architectures_pt_scores)
        max_pt_score_arcis = [arci for arci in cluster_architectures if
                              architectures_pt_scores[cluster_architectures.index(arci)] == max_pt_score]
        known_pts_fused_to_domains = fused_to_domains.intersection(known_pts_set)
        df.loc[len(df)] = [cluster,
                           len(fused_to_domains), fused_to_domains,
                           len(known_pts_fused_to_domains), known_pts_fused_to_domains,
                           architectures_pt_scores, cluster_architectures,
                           max_pt_score, max_pt_score_arcis,
                           sum(architectures_pt_scores)/len(architectures_pt_scores),
                           len(knownPTsArcis), knownPTsArcis, knownPTsArcis_names]

    return df


def getPTStructureScore(arci, known_pts_df):
    arci_ids_list = [dom[0] for dom in arci]
    arci_ids_set = set([dom[0] for dom in arci])
    score = 0

    if terminal == 'N':
        known_pts_ids_set = set(list(known_pts_df[known_pts_df['Type'] != TRAFFICKING]['ID']))
        known_in_arci = arci_ids_set.intersection(known_pts_ids_set)
        idxs_types = [(arci_ids_list.index(dom_id), known_pts_df[known_pts_df['ID'] == dom_id]['Type'].values[0])
                      for dom_id in known_in_arci]
        sorted_idxs_types = sorted(idxs_types, key=lambda x: x[0])
        sorted_types = [tup[1] for tup in sorted_idxs_types]
        sorted_types = [x[0] for x in groupby(sorted_types)]  # remove consecutive duplicates to easier scoring

        if sorted_types == ['R', 'P', 'X']:
            score += 9
        elif sorted_types in [['R', 'X'], ['P', 'X']]:
            score += 6
        elif sorted_types == ['X']:
            score += 3
        elif sorted_types == ['R', 'P']:
            score += 2
        elif sorted_types in [['R'], ['P']]:
            score += 1

    elif terminal == 'C':
        known_pts_ids_set = set(list(known_pts_df[known_pts_df['Type'] != TOXIN]['ID']))
        known_in_arci = arci_ids_set.intersection(known_pts_ids_set)

        idxs_types = [(arci_ids_list.index(dom_id), known_pts_df[known_pts_df['ID'] == dom_id]['Type'].values[0])
                      for dom_id in known_in_arci]
        sorted_idxs_types = sorted(idxs_types, key=lambda x: x[0])
        sorted_types = [tup[1] for tup in sorted_idxs_types]
        sorted_types = [x[0] for x in groupby(sorted_types)]  # remove consecutive duplicates to easier scoring
        if sorted_types == ['F', 'R', 'P']:
            score += 9
        elif sorted_types in [['F', 'R'], ['F', 'P']]:
            score += 6
        elif sorted_types == ['F']:
            score += 3
        elif sorted_types == ['R', 'P']:
            score += 2
        elif sorted_types in [['R'], ['P']]:
            score += 1

    return score


def plotArchitectures(df):
    known_pts_df = pd.read_csv(known_pt_csv_path)
    known_pts_ids_set = set(list(known_pts_df['ID']))
    for idx, row in df.iterrows():
        cluster = row['cluster']
        arcis_df = pd.DataFrame(columns=['ID', 'start', 'end', 'Type', 'arci_idx', 'y'])

        for idx1, arci in enumerate(row['architectures']):
            for domain in arci:
                type = 'PFXXXXX'
                if str(domain[0]) in known_pts_ids_set:
                    type = known_pts_df.loc[known_pts_df['ID'] == domain[0]]['Type'].values[0]
                arcis_df.loc[len(arcis_df)] = [domain[0], domain[1], domain[2], type, idx1, idx1*10 + 5]

        arcis_df['start'] = arcis_df['start'].astype(float)
        arcis_df['end'] = arcis_df['end'].astype(float)
        arcis_df['y'] = arcis_df['y'].astype(float)

        for idx,tup in enumerate(arcis_df.groupby(np.arange(len(arcis_df)) // 10)):
            max_x_coord = max(list(tup[1]['end'])) + 133
            plot = ggplot(data=tup[1]) + xlim(0, max_x_coord) + ggtitle(cluster)
            plot = plot +\
                   geom_segment(aes(x='start', y=1, xend='end', yend=1, color='Type', size=1)) +\
                   facet_wrap('arci_idx', nrow=10)
            plot.save(filename=output_file_path + '_' + str(cluster) + '_' + str(tup[0]) + '.png', format='png')


if __name__ == '__main__':
    annotation_path, known_pt_csv_path, terminal, output_file_path = check_input()
    architectures = make_architectures()
    df = score_modularity()
    df.to_csv(output_file_path, columns=['cluster', 'numOfFusedDomains', 'fusedDomains',
                                         'numOfKnownPTsFusedDomains', 'knownPTsFusedDomains',
                                         'numOfUniqueKnownPTsArchitectures', 'uniqueKnownPTsArchitectures',
                                         'uniqueKnownPTsArchitectures_DOMAINS_NAMES'])
    # plotArchitectures(df)
