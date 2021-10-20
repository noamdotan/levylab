"""
for each cluster hits, writes the names of clusters with no hits to novel_clusters.txt
input_dir = dir with hmm_msas, outpud_dir = dir with hits from pfam
"""

import os
import sys
import time
import subprocess

COVER_PERCENTAGE_THRESHOLD = 0.4
E_VALUE_THRESHOLD = 0.001
NO_HITS_MESSAGE = 'No hits detected that satisfy reporting thresholds'
COMBINED_HMMS_FILE_NAME = 'all_clusters_profile_hmms.hmm'
NOVEL_CLUSTERS_FILE = 'novel_clusters.txt'
LINE_NUMBER_START_OF_HITS = 4
LINE_NUMBER_FROM_END_LAST_LINE_OF_HITS = 11
LINE_NOT_HITS_INIT = '#'


def check_input():
    if len(sys.argv) != 3:
        print(
            'Usage: <input directory path> <Output dir for results>')
        sys.exit(1)
    if not os.path.exists(sys.argv[1]) or \
            not os.path.isdir(sys.argv[2]):
        print('one or more of the paths is invalid')
        sys.exit(1)


def write_novel_clusters(input_dir_path, output_dir_path):
    """for each cluster hits, writes the names of clusters with no hits to novel_clusters.txt"""
    clusters_with_hits = set()
    all_clusters = set()
    with open(os.path.join(output_dir_path, NOVEL_CLUSTERS_FILE), 'w') as novel_clust_file:
        with open(os.path.join(output_dir_path, COMBINED_HMMS_FILE_NAME + '.hits'), 'r') as out_hits_file:
            for line in out_hits_file:
                if line.startswith(LINE_NOT_HITS_INIT):
                    continue
                else:
                    clusters_with_hits.add(line.split()[2][:-4])

        for file in os.listdir(input_dir_path):
            #if file != COMBINED_HMMS_FILE_NAME:
            if file.endswith('.msa.faa'):
                all_clusters.add(file[:file.index('.')])

        clusters_with_no_hits = all_clusters.difference(clusters_with_hits)
        for cluster in clusters_with_no_hits:
            novel_clust_file.write(cluster + '\n')



# def write_novel_clusters(output_dir_path):
#    """for each .out file for a cluster with no hits or only weak hits, outputs it's name to a novel_clusters.txt"""
#    with open(os.path.join(output_dir_path, 'novel_clusters.txt'), 'w') as novel_clust_file:
#        for out_file in [file for file in os.listdir(output_dir_path) if file.endswith('.hits')]:
#            with open(os.path.join(output_dir_path, out_file), 'r') as out_hits_file:
#                for line in out_hits_file:
#                    if re.search(NO_HITS_MESSAGE, line):
#                        novel_clust_file.write(out_file[:out_file.index('.')])
#                        break



if __name__ == '__main__':
    check_input()
    in_dir, out_dir = sys.argv[1], sys.argv[2]
    print('in_dir_path: ' + in_dir)
    print('out_dir_path: ' + out_dir)
    time.sleep(5)
    write_novel_clusters(in_dir, out_dir)
