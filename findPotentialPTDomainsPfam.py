"""this script, given csv of known PTs domains, and given the path to .txt file of pfamA anf a .txt file of all
architectures in pfam, uses them to find other potential PTs domains, by iterating over all known domains from the
csv, and for each domain, iterating over all architectures found in pfam that contain this domain and for each
architecture the script looks for potential PTs domains by examining it's position relative to the current known
domain being used, and it's annotation. i.e if we are using known domain which is a pre toxin then we'll look for (
currently with no consideration in distance from current known domain because i couldn't find the data for it) toxin
annotated domains found after, and multiple repetitive domains found before (up-stream, closer to N-terminal). the
output is all previously known domains and new potential domains according to asaf format """

import os
import sys
import pandas as pd
from itertools import groupby

ID_COL = 'ID'
TYPE_COL = 'Type'
NAME_COL = 'Name'
OCCURRENCES_COL = 'Occurrences_in_IMG_(July 2019)'
TRAFFICKING = 'F'
REPEAT = 'R'
PRETOXIN = 'P'
TOXIN = 'X'
TYPES = list([TRAFFICKING, REPEAT, PRETOXIN, TOXIN])
PT_KEYWORDS = ['polymorphic']
TRAFFICKING_KEYWORDS = [['bind', 'binding'],
                        ['trafficking', 'mediate', 'n-terminal', 'bind', 'binding']]
REPEAT_KEYWORDS = [['repeat'],
                   ['repeat', 'repeating', 'repetitive']]
PRETOXIN_KEYWORDS = [['pt', 'pre-toxin', 'pretoxin'],
                     ['pre-toxin', 'pretoxin']]
TOXIN_KEYWORDS = [['nuclease', 'dnase', 'rnase', 'deaminase', 'deaminases', 'tox', 'endonuclease'],  # name keywords
                  ['toxin', 'nuclease', 'dnase', 'rnase', 'deaminase', 'endonuclease']]  # description keywords
MIN_REPEATS = 3
MIN_TYPE_SCORE = 1.5
MAX_DOMAINS_PRETOXIN_FROM_TOXIN = 1


def check_input():
    """
    :return: tuple of (known_PTs_domains_csv_path, pfamA path, pfam_architectures_path,
     output_domains_data_path, output_IMG_occurences_path)
    """
    if len(sys.argv) != 7:
        print(
            'Usage: <Known PTs domains csv path> <pfamA path> <Architectures in pfam path> <IMG path>' +
            '<Output for domains data path> <Output for IMG occurrences path>')
        sys.exit(1)
    if not os.path.exists(sys.argv[1]) or \
            not os.path.exists(sys.argv[2]) or \
            not os.path.exists(sys.argv[3]) or \
            not os.path.exists(sys.argv[4]):
        print('one or more of the paths is invalid')
        sys.exit(1)
    if os.path.exists(sys.argv[5]) or os.path.exists(sys.argv[6]):
        if input('Be careful! output file already exists, are you sure you want to override? (any key/no) ') == 'no':
            sys.exit(1)

    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]


def get_types_score(name, description):
    scores = list([0, 0, 0, 0])  # trafficking, repeat, pretoxin, toxin
    for i, keywords in enumerate([TRAFFICKING_KEYWORDS, REPEAT_KEYWORDS, PRETOXIN_KEYWORDS, TOXIN_KEYWORDS]):
        # i decided that it is more significant if a keyword is in name than in description
        scores[i] = len(set(keywords[0]).intersection(set(name))) * 1.5 + \
                    len(set(keywords[1]).intersection(set(description))) * 0.5
    #  if description contains the 2 words - 'polymorphic toxin\s' add to score
    if ('polymorphic' in description and 'toxin' in description and
        description.index('polymorphic') + 1 == description.index('toxin')) or \
            ('polymorphic' in description and 'toxins' in description and
             description.index('polymorphic') + 1 == description.index('toxins')):

        for j in range(len(scores)):
            scores[j] += 1.5
    return scores


def get_domain_type(domain_id, get_scores=False):
    """returns domain type according to pfamA description of domain family. returns nan if type could not be decided"""
    scores = [0, 0, 0, 0]  # TRAFFICKING, REPEAT, PRETOXIN, TOXIN
    pfamA = open(pfamA_path, 'r')
    for row in pfamA:
        split = row.split('\t')
        if domain_id == split[0]:
            name = split[3].lower().split()
            description = split[7].lower().split()
            scores = get_types_score(name, description)
            break
    pfamA.close()
    if max(scores) >= MIN_TYPE_SCORE:
        # if there are several types with maximum score, we dont know type of domain
        if len([i for i, score in enumerate(scores) if score == max(scores)]) == 1:
            return TYPES[scores.index(max(scores))] if not get_scores else TYPES[scores.index(max(scores))], scores
    return float('nan') if not get_scores else float('nan'), scores


def get_num_repeats(repeat_domain_id, arc_domains_ids_ordered):
    """returns maximum number of repeats of repeat_domain_id and start and end index"""
    grouped = [(k, sum(1 for i in g)) for k, g in groupby(arc_domains_ids_ordered)]
    max_r = 0
    start = 0
    max_r_start = 0
    for d_id, n in grouped:
        if d_id == repeat_domain_id and n > max_r:
            max_r = n
            max_r_start = start
        start += n
    return max_r, max_r_start, max_r_start + max_r - 1


def check_if_potential_PT_domain(cur_domain_id, cur_domain_type, domain_id, domain_id_idx, arc_domains_ids_ordered):
    """given current known domain being used id and type, domain pfam id and an architecture it is found in, returns
    if domain_id is potential PT domain and it's type (trafficking, repeat, pre-toxin, toxin) as a tuple in the format
    of (type, isPTDomain). This is done according to domains positions and common structure of polymorphic toxins.
    if domain type could not be decided, returns (nan, False)"""
    domain_type, domains_type_scores = get_domain_type(domain_id, True)
    is_potential_PT_domain = False

    # if we managed to infer potential domain type follow this logic based on type and position
    if pd.notna(domain_type):

        # if cur_domain is trafficking and N-terminal we accept if potential domain is repeat and appears downstream and
        # MIN_REPEATS consecutive repeats and starts after N-terminal and ends before C-terminal,
        # or if domain is toxin and C-terminal
        if cur_domain_type == TRAFFICKING and arc_domains_ids_ordered[0] == cur_domain_id and \
                domain_type != cur_domain_type:
            if domain_type == REPEAT:
                num_max_repeats, start, end = get_num_repeats(domain_id, arc_domains_ids_ordered)
                if num_max_repeats >= MIN_REPEATS and \
                        start > 0 and end < len(arc_domains_ids_ordered) - 1:
                    is_potential_PT_domain = True

            #if domain_type == REPEAT and 0 < domain_id_idx < len(arc_domains_ids_ordered) - 1:
            #    is_potential_PT_domain = True

            #elif domain_type == PRETOXIN and domain_id in \
            #        arc_domains_ids_ordered[-min(MAX_DOMAINS_PRETOXIN_FROM_TOXIN, len(arc_domains_ids_ordered) - 1):]\
            #        and get_domain_type(arc_domains_ids_ordered[-1]) == TOXIN:
            #    is_potential_PT_domain = True

            elif domain_type == TOXIN and arc_domains_ids_ordered[-1] == domain_id:
                is_potential_PT_domain = True

        # if cur_domain is repeat we accept if potential domain is trafficking and appears upstream and N-terminal
        # or is toxin and c-terminal downstream
        elif cur_domain_type == REPEAT and domain_type != cur_domain_type:

            if domain_type == TRAFFICKING and arc_domains_ids_ordered[0] == domain_id:
                is_potential_PT_domain = True

            #elif domain_type == PRETOXIN:
            #    if arc_domains_ids_ordered.index(cur_domain_id) < domain_id_idx:
            #        is_potential_PT_domain = True

            elif domain_type == TOXIN and arc_domains_ids_ordered[-1] == domain_id:
                is_potential_PT_domain = True

        # if cur_domain is pre-toxin we accept if potential domain appears downstream and it is of type toxin
        # and domains that appear upstream and are repeats or trafficking
        elif cur_domain_type == PRETOXIN and domain_type != cur_domain_type:
            # if domain_type == REPEAT:
            #    num_max_repeats, start, end = get_num_repeats(domain_id, arc_domains_ids_ordered)
            #    if num_max_repeats >= MIN_REPEATS and end < arc_domains_ids_ordered.index(cur_domain_id):
            #        is_potential_PT_domain = True
            if domain_type == TOXIN and arc_domains_ids_ordered[-1] == domain_id:
                if cur_domain_id in arc_domains_ids_ordered[-min(MAX_DOMAINS_PRETOXIN_FROM_TOXIN + 1, len(arc_domains_ids_ordered)):]:
                    is_potential_PT_domain = True

            elif domain_type == TRAFFICKING and arc_domains_ids_ordered[0] == domain_id:
                is_potential_PT_domain = True

        # if cur_domain is toxin we accept if potential domain is N-terminal Trafficking or close enough Pretoxin
        # or repeat that appears at least MIN_REPEAT consecutive repeats and
        # starts after at the N-terminal and ends before the C-terminal
        elif cur_domain_type == TOXIN and arc_domains_ids_ordered[-1] == cur_domain_id and\
                domain_type != cur_domain_type:
            if domain_type == REPEAT:
                num_max_repeats, start, end = get_num_repeats(domain_id, arc_domains_ids_ordered)
                if num_max_repeats >= MIN_REPEATS and start > 0 and end < len(arc_domains_ids_ordered) - 1:
                    is_potential_PT_domain = True

            elif domain_type == TRAFFICKING and arc_domains_ids_ordered[0] == domain_id:
                is_potential_PT_domain = True

            #elif domain_type == REPEAT and 0 < domain_id_idx < len(arc_domains_ids_ordered) - 1:
            #    is_potential_PT_domain = True

            elif domain_type == PRETOXIN and domain_id in \
                arc_domains_ids_ordered[-min(MAX_DOMAINS_PRETOXIN_FROM_TOXIN + 1, len(arc_domains_ids_ordered)):]:
                is_potential_PT_domain = True

    # if we couldn't infer potential domain type follow this logic based mostly on position
    else:
        # if potential domain is N-terminal with a chance to be trafficking,
        # and cur_domain is C-terminal toxin we accept as trafficking
        if arc_domains_ids_ordered[0] == domain_id and domains_type_scores[0] >= 1 \
                and cur_domain_type == TOXIN and arc_domains_ids_ordered[-1] == cur_domain_id:
            domain_type, is_potential_PT_domain = TRAFFICKING, True

        # if potential domain is C-terminal with chance to be toxin and cur_domain is upstream adjacent Pre-toxin or
        # Trafficking N-terminal we accept as toxin
        if arc_domains_ids_ordered[-1] == domain_id and domains_type_scores[3] >= 1:
            if (cur_domain_type == TRAFFICKING and arc_domains_ids_ordered[0] == cur_domain_id) or \
               (cur_domain_type == PRETOXIN and arc_domains_ids_ordered[len(arc_domains_ids_ordered) - 1] == cur_domain_id):
                domain_type, is_potential_PT_domain = TOXIN, True

    return domain_type, is_potential_PT_domain


def get_architectures_containing_domain_id(domain_id, domain_type, architectures_path):
    # todo right now if repeat accept an architecture even if theres only one appearance of domain id
    #  ask asaf if to add functionality of accepting architecture for repeats only if appears MIN_REPEATS
    #  this functionality is commented below right now
    arcs = list()
    with open(architectures_path, 'r') as all_pfam_architectures:
        for row in all_pfam_architectures:
            good_row = row[:-1]  # removing new line char
            good_row_split = list(good_row.split('\t')[4].split(' '))
            if domain_id in good_row_split:
                # if domain_type == REPEAT:
                # if any([domain_id for i in range(MIN_REPEATS)] == good_row_split[i:i+MIN_REPEATS]
                #       for i in range(len(good_row_split) - MIN_REPEATS)):
                arcs.append(good_row)
                # else:
                #    arcs.append(good_row)
    return arcs


def split_architecture(architecture):
    """splitting architecture to it's inorder domains names, inorder domains ids
    and the protein with this architecture"""
    split = architecture.split('\t')
    return split[1].split('~'), split[2].split(' '), split[4].split(' ')


def search_new_domains_by_cur_domain_id_type(all_ids_set, cur_domain_id, cur_domain_type, df_columns,
                                             architectures_path):
    """returns pandas data frame containing all found potential PT domains in pfam
     according to a known domain id and type (allowing duplicates by id)"""
    new_domains = pd.DataFrame(columns=df_columns)
    architectures = get_architectures_containing_domain_id(cur_domain_id, cur_domain_type, architectures_path)
    for architecture in architectures:

        print(architectures.index(architecture), ' out of ', len(architectures))

        arc_domains_names, arc_protein_id, arc_domains_ids = split_architecture(architecture)
        for potential_domain_id_idx, potential_domain_id in enumerate(arc_domains_ids):
            if potential_domain_id != cur_domain_id and \
                    potential_domain_id not in all_ids_set and \
                    potential_domain_id not in all_potentially_found_domains_ids_used:  # checking redundancy
                potential_domain_type, is_PT_domain = check_if_potential_PT_domain(cur_domain_id, cur_domain_type,
                                                                                   potential_domain_id,
                                                                                   potential_domain_id_idx,
                                                                                   arc_domains_ids)
                # if decided to be PT domain
                # adding domain id, type and architecture and cur_domain_id the decision based on
                if is_PT_domain:
                    new_domains.loc[len(new_domains)] = [arc_domains_names[potential_domain_id_idx],
                                                         potential_domain_id,
                                                         float('nan'),
                                                         float('nan'),
                                                         potential_domain_type,
                                                         architecture,
                                                         cur_domain_id]

    new_domains.drop_duplicates(subset=[ID_COL], inplace=True)
    new_domains.reset_index(inplace=True, drop=True)
    return new_domains


def search_new_domains():
    """return df of previously known and new found potential PT domains according to script description logic"""
    # todo decide if to allow duplicates or not, right now not allowing
    # loading known domains
    df = pd.read_csv(known_PTs_domains_csv_path)

    # initializing with known domains
    all_domains = df
    all_domains_ids_set = set(all_domains[ID_COL])  # for search stop condition

    # exhaustive search until no new unique domains ids are found
    prev_num_uniq_domain_ids = len(all_domains_ids_set) - 1
    prev_domains_ids = set()
    i = 0
    while len(all_domains_ids_set) > prev_num_uniq_domain_ids:
        i += 1
        new_domains = pd.DataFrame(columns=df.columns)  # creating empty df for new domains to be found
        cur_iter_prev_domains_ids = set()  # for not searching again if duplicates appear
        for index in all_domains.index:
            print(i, index, all_domains[ID_COL][index], all_domains[NAME_COL][index])
            cur_domain_name, cur_domain_id, cur_domain_type = all_domains[NAME_COL][index], \
                                                              all_domains[ID_COL][index], \
                                                              all_domains[TYPE_COL][index]
            # if undefined type or id or we searched this cur_domain_id in a prev iter, skip this one
            if cur_domain_id in all_potentially_found_domains_ids_used or \
                    cur_domain_id in prev_domains_ids or \
                    cur_domain_id in cur_iter_prev_domains_ids or \
                    pd.isna(cur_domain_id) or \
                    cur_domain_type not in TYPES:
                continue
            all_potentially_found_domains_ids_used.add(cur_domain_id)
            cur_iter_prev_domains_ids.add(cur_domain_id)
            tmp = search_new_domains_by_cur_domain_id_type(all_domains_ids_set, cur_domain_id, cur_domain_type,
                                                           list(df.columns), pfam_architectures_path)
            new_domains = new_domains.append(tmp, ignore_index=True)  # adding all new domains using cur_domain
            print('from ' + str(cur_domain_id) + ' got ' + str(len(tmp)) + ' new domains')
        prev_domains_ids = all_domains_ids_set
        all_domains = all_domains.append(new_domains, ignore_index=True)  # adding all new domains from this iteration
        all_domains_ids_set = set(all_domains[ID_COL])  # for search stop condition
        # removing duplicates by ID only
        all_domains.drop_duplicates(subset=[ID_COL], inplace=True)
        all_domains.reset_index(inplace=True, drop=True)
        prev_num_uniq_domain_ids = len(prev_domains_ids)

    # removing duplicates by ID only
    all_domains.drop_duplicates(subset=[ID_COL], inplace=True)
    all_domains.reset_index(inplace=True, drop=True)
    return all_domains


def get_pfamA_name_by_domain_id(domain_id):
    pfamA = open(pfamA_path, 'r')
    name = None
    for row in pfamA:
        split = row.split('\t')
        if domain_id == split[0]:
            name = split[3]
            break
    pfamA.close()
    return name


def get_IMG_occur_by_domain_id(domain_id):
    """"returns tuple (num_occurrences, df of (IMG genome, IMG gene))"""
    num_occur = 0
    genomes_genes = pd.DataFrame(columns=['GENOMEID', 'GENEID'])
    genomes = os.listdir(IMG_path)

    for genome in genomes:
        filename = genome + '.pfam.tab.txt'
        path = os.path.join(IMG_path, genome, filename)

        if not os.path.exists(path):
            continue
        genes_pfam = open(path, 'r')

        for row in genes_pfam:
            pfam_id = row.split('\t')[8]
            if len(pfam_id) < 5:
                continue
            PF_id = pfam_id[:2].upper() + pfam_id[4:]
            if domain_id == PF_id:
                num_occur += 1
                genomes_genes.loc[len(genomes_genes)] = [genome, row.split('\t')[0]]

        genes_pfam.close()

    # num_occur = float('nan') if num_occur == 0 else num_occur
    return num_occur, genomes_genes


def put_names_and_IMG_occurrences_rm_dup(new_domains):
    """put info about name and IMG occur and remove duplicates in pfamA ID
    in addition returns all locations of occurrences of each of the domains in new_domains in IMG"""

    df = new_domains
    all_genomes_genes = pd.DataFrame(columns=['GENOMEID', 'GENEID'])

    # remove duplicates accroding to pfam ID and domain obtained from
    df.drop_duplicates(subset=[ID_COL], inplace=True)
    df.reset_index(inplace=True, drop=True)

    for index in df.index:
        if pd.isna(df[NAME_COL][index]):
            df[NAME_COL][index] = get_pfamA_name_by_domain_id(df[ID_COL][index])
        if pd.isna(df[OCCURRENCES_COL][index]):
            df[OCCURRENCES_COL][index], tmp = get_IMG_occur_by_domain_id(df[ID_COL][index])
            all_genomes_genes = all_genomes_genes.append(tmp)
    return df, all_genomes_genes


if __name__ == "__main__":
    known_PTs_domains_csv_path, pfamA_path, pfam_architectures_path, IMG_path, output_file_path, IMG_occur_path \
        = check_input()
    all_potentially_found_domains_ids_used = set()

    known_and_new_domains = search_new_domains()
    #  todo right now not calculating IMG occurrences automatically
    # known_and_new_domains, IMG_occurrences_locations = put_names_and_IMG_occurrences_rm_dup(known_and_new_domains)
    known_and_new_domains.to_csv(output_file_path, index=False, encoding='utf-8')
    # todo right now not writing to csv locations because it's not in the desired logic. desired logic implemented in script
    #  'getIMGOccurrences'
    # IMG_occurrences_locations.to_csv(IMG_occur_path, index=False, encoding='utf-8')
