"""
this script gets:
- fasta file of single gene domain #todo what should be the format of description
- diamond hits file of the domain containing sequences
- hmmscan output file of the above diamond hits - hmmscan NEEDS to be against PFAM + PredPTs
- fasta file of all unannotated genes (proteins) containing the above domain??
- maybe IMG metadata / taxon counts of domain??

Than the script creates plots of the domain:
- architectures
- conserved sequences (analysis by MSA - global?local?)
- phylogenetic tree?

"""
#todo MAKE SURE an ARCHITECTURE is PLOTTED ONLY IF WE RECOGNIZED THE PredPT in it with HMMSCAN

from collections import defaultdict
from itertools import groupby
import sys
import os
import pandas as pd
from plotnine import *


PRED_KEYWORDS = ['pred', 'cluster', 'nterm', 'predptraf']

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
COVER_PERCENTAGE_THRES = 0.5
#ACC_THRES = #0.85
FRACTION_CUTOFF = 0.2
OVERLAP_THRES_AA = 50


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

TXT_ARCIS_FILENAME = 'text_printed_architectures.txt'
PLOT_ARCIS_FILENAME = 'architectures.png'
JSON_FILENAME = 'architecture.json'
COLOR_MAP = {'F': '#5e3c99',
             'R': '#fdb863',
             'P': '#b2abd2',
             'X': '#e66101',
             'PredPT': '#e66101',
             'Other': '#fddbc7'}
REPEATS_COLORS = {'PF05860': '#e5f5e0',
                  'PF13332': '#c7e9c0',
                  'PF03527': '#a1d99b',
                  'PF04830': '#74c476',
                  'PF03752': '#41ab5d',
                  'PF05593': '#238b45',
                  'PF05594': '#006d2c',
                  'PF13517': '#00441b'}

SCALE = 2.5
FIXED_IN_BETWEEN_LENGTH = 50


def check_input():
    if len(sys.argv) != 6:
        print(
            'Usage: <path to hmmscan annotations file> < known PTs domains> <N or C Terminal> <Output dir for architectures and jsons> <is only known pt domains (T,F)>')
        sys.exit(1)
    if not os.path.exists(sys.argv[1]) or not os.path.exists(sys.argv[2]):
        print('one or more of the paths is invalid')
        sys.exit(1)
    if sys.argv[3] != N_TERMINAL_ARG and sys.argv[3] != C_TERMINAL_ARG:
        print('Terminal arg has to be N or C')
        sys.exit(1)
    #if os.path.exists(sys.argv[4]):
    #   if input('Be careful! output file already exists, are you sure you want to override? (any key/no) ') == 'no':
    #        sys.exit(1)
    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]


def calc_overlap(from_1, to_1, from_2, to_2):
    return max(0, min(to_1, to_2) - max(from_1, from_2))


def domain_align_is_ok(evalue, dom_len, hmm_from, hmm_to, acc, protein_from, protein_to, protein_len):
    if evalue <= EVAL_THRES and \
            ((hmm_to - hmm_from) / dom_len) >= COVER_PERCENTAGE_THRES: # and \
            #acc >= ACC_THRES:
        #if (terminal == N_TERMINAL_ARG and protein_from >= TERMINAL_LENGTH - N_FROM) or \
        #        (terminal == C_TERMINAL_ARG and protein_to <= protein_len - TERMINAL_LENGTH + C_TO):
            return True
    return False


def insert_domain_to_architecture(PredPT, protein, protein_length, domain, architectures_dict):
    if PredPT not in architectures_dict.keys():
        architectures_dict[PredPT].append((protein, protein_length, list()))
        architectures_dict[PredPT][0][2].append(domain)
    else:
        for idx, arci in enumerate(architectures_dict[PredPT]):
            if arci[0] == protein:
                for i, dom in enumerate(arci[2]):
                    if calc_overlap(domain[3], domain[4], dom[3], dom[4]) > OVERLAP_THRES_AA:  # first check if there is overlap case
                        if dom[-1] < domain[-1]:  # if current eval smaller, do nothing
                            return
                        else:  # if new eval smaller, replace current with it
                            arci[2][i] = domain
                            return
                    if dom[3] > domain[4] or dom[3] > domain[3]:  # if not overlappig check regular conditions
                        arci[2].insert(i, domain)
                        return
                arci[2].append(domain)  # if it's downstream from all other domains
                return
        # if it's arci of a new protein we add arci of this protein and add this domain to the arci
        architectures_dict[PredPT].append((protein, protein_length, list()))
        architectures_dict[PredPT][-1][2].append(domain)


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
            PredPT_genome_gene_protein = split[CLUSTER_GENOME_GENE_PROTEIN_IDX]
            PredPT = PredPT_genome_gene_protein.split('~')[0]
            if len(PredPT_genome_gene_protein.split('~')) > 1:
                genome_gene = PredPT_genome_gene_protein.split('~')[1]
            if len(PredPT_genome_gene_protein.split('~')) > 2:
                protein = PredPT_genome_gene_protein.split('~')[2]
            else:
                protein = PredPT_genome_gene_protein.split('~')[1]
            protein_length = int(split[PROTEIN_LENGTH_IDX])
            evalue = float(split[DOM_E_VALUE_IDX])
            hmm_from = int(split[HMM_FROM_IDX])
            hmm_to = int(split[HMM_TO_IDX])
            protein_from = int(split[PROTEIN_FROM_IDX])
            protein_to = int(split[PROTEIN_TO_IDX])
            acc = float(split[ACC_IDX])
            if domain_align_is_ok(evalue, dom_length, hmm_from, hmm_to, acc, protein_from, protein_to, protein_length):
                domain = (dom_id, dom_name, dom_length, protein_from, protein_to, hmm_from, hmm_to, evalue)
                insert_domain_to_architecture(PredPT, protein, protein_length, domain, architectures_dict)

    return architectures_dict


def print_architectures(architectures):
    # output to txt file
    with open(os.path.join(output_dir_path, TXT_ARCIS_FILENAME), 'w') as out_file:
        for predPT in architectures.keys():
            for (protein, protein_len, arci) in architectures[predPT]:
                protein_len = int(protein_len / 5)
                arci_str_names = ['-' for i in range(protein_len)]
                arci_str_ids = ['-' for i in range(protein_len)]
                predPTFlag = False
                for (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval) in arci:
                    for keyword in PRED_KEYWORDS:
                        if keyword in domain_name.lower():
                            predPTFlag = True
                            break
                    protein_from = int(protein_from / 5)
                    protein_to = int(protein_to / 5)
                    arci_str_names[protein_from-1] = '('
                    arci_str_names[protein_to-1] = ')'
                    arci_str_ids[protein_from - 1] = '('
                    arci_str_ids[protein_to - 1] = ')'

                    mid = int((protein_from + protein_to) / 2 - len(domain_name) / 2)
                    for i, letter in enumerate(list(domain_name)):
                        arci_str_names[mid + i + 1] = letter
                    for i, letter in enumerate(list(domain_id)):
                        arci_str_ids[mid + i] = letter

                if predPTFlag:
                    out_file.write(''.join(arci_str_names) + '\n')
                    out_file.write(''.join(arci_str_ids) + '\n\n')


def plot_architectures(architectures):

    for predPT in architectures.keys():
        # to add only one of each architecture structur
        added_arcis_ids = set()
        # creating data frame to plot easily
        df = pd.DataFrame(data={}, columns=['PredPT',
                                            'protein_idx', 'protein_name', 'protein_len',
                                            'dom_id', 'dom_name', 'domain_length', 'protein_from', 'protein_to',
                                            'arci_its_from_names', 'arci_its_from_ids', 'arci_its_from_count'])
        for idx, (protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count) in enumerate(architectures[predPT]):
            predPTFlag = False
            if arci_str_ids in added_arcis_ids:
                continue
            else:
                added_arcis_ids.add(arci_str_ids)
            for (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval) in arci:
                for keyword in PRED_KEYWORDS:
                    if keyword in domain_name.lower():
                        predPTFlag = True
                        break
                df.loc[len(df)] = [predPT,
                                    idx, protein, protein_len,
                                    domain_id, domain_name, domain_length, protein_from, protein_to,
                                    arci_str_names, arci_str_ids, 'X'+str(arci_count)]
            # for not including an architecture if for some reason hmmscan didnt detect the prePT in it
            # (although it should've..)
            if not predPTFlag:
                df.drop(df.index[df['protein_name'] == protein], inplace=True)
                df.reset_index(drop=True, inplace=True)

        # sorting accroding to protein_len
        df.sort_values(['protein_len', 'protein_name', 'protein_from'], inplace=True)
        # sorting according to abundance of architecture
        df['arci_its_from_count_int'] = df['arci_its_from_count'].str[1:].astype(int)
        df.sort_values(['arci_its_from_count_int'], ascending=False, inplace=True)
        df['protein_name'] = pd.Categorical(df['protein_name'], categories=reversed(df['protein_name'].unique()), ordered=True)
        df.reset_index(drop=True)
        df['text_x_pos'] = (df['protein_from'] + df['protein_to']) / 2
        df['protein_len'], df['protein_from'], df['protein_to'], df['text_x_pos'] = df['protein_len'].astype(int), \
                                                                                    df['protein_from'].astype(int), \
                                                                                    df['protein_to'].astype(int), \
                                                                                    df['text_x_pos'].astype(int)
        # plotting
        y_axis_value = 'protein_name'

        p = ggplot(df, aes(x='protein_from', xend='protein_to', y=y_axis_value, color='dom_name')) +\
            geom_segment(aes(x=0, y=y_axis_value, xend='protein_len', yend=y_axis_value), size=1) +\
            geom_segment(aes(x='protein_from', y=y_axis_value, xend='protein_to', yend=y_axis_value), size=10) +\
            geom_label(aes(label='dom_name', x='text_x_pos'), size=7) +\
            geom_label(aes(label='arci_its_from_count', x='protein_len', size=9)) +\
            theme(panel_grid_minor=element_blank(),
                  panel_grid_major=element_blank(),
                  legend_position='none') + \
            ylab(ylab='representative_protein_with_this_architecture') + \
            xlab(xlab='protein_aa_position') + \
            ggtitle(title='Architectures of ' + str(predPT))
        print(p)
        p.save(filename=os.path.join(output_dir_path, str(predPT) + '_' + PLOT_ARCIS_FILENAME),
               format='png', height=40, width=20, limitsize=False)


def count_unique_architectures(arcis):
    arcis_with_counts = defaultdict(list)  # with additional cells in list - str_name_arci, str_id_arci, count
    for predPT in arcis.keys():
        unique_arcis_names = defaultdict(int)
        unique_arcis_ids = defaultdict(int)
        for (protein, protein_len, arci) in arcis[predPT]:
            arci_str_names = (''.join([dom_name+'~' for (dom_id, dom_name, dom_len, p_from, p_to, hmm_from, hmm_to, eval) in arci]))[:-1]
            arci_str_ids = (''.join([dom_id + '~' for (dom_id, dom_name, dom_len, p_from, p_to, hmm_from, hmm_to, eval) in arci]))[:-1]
            unique_arcis_names[arci_str_names] += 1
            unique_arcis_ids[arci_str_ids] += 1

        for (protein, protein_len, arci) in arcis[predPT]:
            arci_str_names = (''.join([dom_name + '~' for (dom_id, dom_name, dom_len, p_from, p_to, hmm_from, hmm_to, eval) in arci]))[:-1]
            arci_str_ids = (''.join([dom_id + '~' for (dom_id, dom_name, dom_len, p_from, p_to, hmm_from, hmm_to, eval) in arci]))[:-1]
            arcis_with_counts[predPT].append(
                (protein, protein_len, arci, arci_str_names, arci_str_ids, unique_arcis_names[arci_str_names]))
    return arcis_with_counts


def get_domain_color(domain_id, domain_name, known_pts_df):
    type = known_pts_df[known_pts_df['ID']==domain_id]['Type'].values
    if len(type) == 0:
        type ='Other'
    else:
        type = type[0]
    if 'PredPT' in domain_name:
        type = 'X'
    if 'PredPTraf' in domain_name:
        type='F'
    color = COLOR_MAP[type]
    #if type == 'R': # getting the right shade of green if it's repeat domain
    #    color = REPEATS_COLORS[domain_id]
    return color


def get_style(hmm_from, hmm_to, domain_length, isStart):
    style = 'curved'
    if isStart:
        if (hmm_from / domain_length) >= FRACTION_CUTOFF:
            style = 'jagged'
    else:
        if (domain_length - hmm_to) / domain_length >= FRACTION_CUTOFF:
            style = 'jagged'
    return style


def build_pfam_json_from_architecture(predPT, protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count):
    known_pts_df = pd.read_csv(known_pt_csv_path, encoding='iso-8859-1')
    with open(os.path.join(output_dir_path, 'jsons', str(predPT) + '_' + str(protein) + '_' + JSON_FILENAME), 'w') as json:
        json.write('{ \n')
        json.write('  "length" : ' + '"'+str(protein_len)+'",  \n')
        json.write('  "regions" : [  \n')
        for i, (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval) in enumerate(arci):
            color = get_domain_color(domain_id, domain_name, known_pts_df)
            start_style = get_style(hmm_from, hmm_to, domain_length, True)
            end_style = get_style(hmm_from, hmm_to, domain_length, False)
            # start of domain json
            json.write('    {  \n')
            json.write('      "type" : "pfama",  \n')
            json.write('      "text" : "'+str('RHS_r' if domain_name == 'RHS_repeat' else domain_name)+'",  \n')
            json.write('      "colour" : "' + str(color) + '",  \n')
            json.write('      "display" : "true", \n')
            json.write('      "startStyle" : "' + str(start_style) + '", \n')
            json.write('      "endStyle" : "' + str(end_style) + '", \n')
            json.write('      "start" : "' + str(protein_from) + '", \n')
            json.write('      "end" : "' + str(protein_to) + '",  \n')
            json.write('      "aliEnd" : "' + str(protein_to) + '",  \n')
            json.write('      "aliStart" : "' + str(protein_from) + '",\n')
            # end of domain json
            json.write('    }')
            if i < len(arci_str_names) - 1:
                json.write(', \n')
            else:
                json.write(' \n')
        json.write('  ] \n')
        json.write('}')


def build_pfam_jsons(architectures_with_counts):
    for predPT in architectures_with_counts.keys():
        created = set()
        for (protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count) in architectures_with_counts[predPT]:
            if arci_str_ids not in created:
                created.add(arci_str_ids)
                build_pfam_json_from_architecture(predPT, protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count)


def print_architectures_simple(architectures_with_counts):  #todo add argument of immunity downstream and print as well
    for predPT in architectures_with_counts.keys():
        with open(os.path.join(output_dir_path, 'txt_architectures', str(predPT) + '_' + TXT_ARCIS_FILENAME), 'w') as txt:
            printed = set()
            for (protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count) in architectures_with_counts[predPT]:
                if arci_str_ids not in printed:
                    printed.add(arci_str_ids)
                    split = arci_str_names.split('~')
                    split = [(i, sum(1 for i in group)) for i, group in groupby(split)]
                    txt.write(protein + '\n')
                    for i, (dom, dups) in enumerate(split):
                        txt.write(dom)
                        txt.write(' X '+str(dups) if dups > 1 else '')
                        txt.write('~' if i < len(split)-1 else '')
                    txt.write('\t'+str(arci_count)+' occurrences in db')
                    txt.write('\n\n\n')


def scale_up_alignments_to_fit_text(architectures_with_counts):
    new_dict = dict()
    for predPT in architectures_with_counts.keys():
        new_dict[predPT] = list()
        for (protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count) in architectures_with_counts[predPT]:
            new_arci = list()
            for (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval) in arci:
                domain = (domain_id, domain_name, domain_length*SCALE, protein_from*SCALE, protein_to*SCALE, hmm_from*SCALE, hmm_to*SCALE, eval)
                new_arci.append(domain)
            new_dict[predPT].append((protein, protein_len*SCALE, new_arci, arci_str_names, arci_str_ids, arci_count))
    return new_dict


def shrink_architectures(architectures_with_counts):
    new_dict = dict()
    for predPT in architectures_with_counts.keys():
        new_dict[predPT] = list()
        for (protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count) in architectures_with_counts[predPT]:
            new_arci = list()
            prev_dom_end = 0
            diff = 0
            total_shrink = 0
            prev_shrink = 0
            for idx, (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval) in\
                    enumerate(arci):
                diff = protein_from - prev_dom_end
                shrinked_protein_from = protein_from - diff + FIXED_IN_BETWEEN_LENGTH - total_shrink\
                    if (diff > FIXED_IN_BETWEEN_LENGTH) else protein_from - total_shrink
                shrinked_protein_to = protein_to - diff + FIXED_IN_BETWEEN_LENGTH - total_shrink\
                    if (diff > FIXED_IN_BETWEEN_LENGTH)  else protein_to - total_shrink
                domain = (domain_id, domain_name, domain_length, shrinked_protein_from, shrinked_protein_to,
                          hmm_from, hmm_to, eval)
                new_arci.append(domain)
                prev_dom_end = protein_to
                prev_shrink = (diff - FIXED_IN_BETWEEN_LENGTH) \
                    if (diff > FIXED_IN_BETWEEN_LENGTH) else 0
                total_shrink += prev_shrink

            new_dict[predPT].append((protein, protein_len - total_shrink, new_arci, arci_str_names, arci_str_ids, arci_count))
    return new_dict


def clean_no_pt_domains(architectures_with_counts):
    known_pts_ids_set = set(pd.read_csv(known_pt_csv_path, encoding='iso-8859-1')['ID'])
    new_dict = dict()
    for predPT in architectures_with_counts.keys():
        new_dict[predPT] = list()
        for (protein, protein_len, arci, arci_str_names, arci_str_ids, arci_count) in architectures_with_counts[predPT]:
            new_arci = list()
            for (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval) in arci:
                domain = (domain_id, domain_name, domain_length, protein_from, protein_to, hmm_from, hmm_to, eval)
                if domain_id in known_pts_ids_set\
                        or any([keyword.lower() in domain_name.lower() for keyword in PRED_KEYWORDS]):
                    new_arci.append(domain)
            new_dict[predPT].append((protein, protein_len, new_arci, arci_str_names, arci_str_ids, arci_count))
    return new_dict


if __name__ == '__main__':
    annotation_path, known_pt_csv_path, terminal, output_dir_path, is_only_known = check_input()
    is_only_known = True if is_only_known == 'T' else False
    #known_pts_df = pd.read_csv(known_pt_csv_path, encoding='iso-8859-1')
    #repeat_domains = list(known_pts_df[known_pts_df['Type'] == 'R'])
    os.system('mkdir '+os.path.join(output_dir_path, 'jsons'))
    os.system('mkdir '+os.path.join(output_dir_path, 'txt_architectures'))
    architectures = make_architectures()
    architectures_with_counts = count_unique_architectures(architectures)

    architectures_with_counts_onlyPTdoms = clean_no_pt_domains(architectures_with_counts)

    arcis_to_plot = architectures_with_counts_onlyPTdoms if is_only_known else architectures_with_counts

    plot_architectures(arcis_to_plot)

    build_pfam_jsons(shrink_architectures(scale_up_alignments_to_fit_text(arcis_to_plot)))

    print_architectures_simple(arcis_to_plot)
