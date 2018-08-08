#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song
# songwz03@gmail.com

# get_purity.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# get_purity.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import glob
import random
import shutil
import argparse
import numpy as np
from time import sleep
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


usage = '''

python3 get_purity.py -r reference_folder -q query.fasta -x fasta -n 100000 -p 90 -bbmap path/to/bbmap.sh
python3 ~/PycharmProjects/PacBio_script/get_purity.py -r reference_genomes -q query.fas -n 10000 -x fasta -bbmap ~/Softwares/bbmap/bbmap.sh 

'''


def get_total_length(seq_file):
    total_length = 0
    for each_seq in SeqIO.parse(seq_file, 'fasta'):
        total_length += len(each_seq.seq)
    return total_length


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def simulate_reads(pwd_ref_file, read_number, read_length, insert_size, pwd_output_folder):

    # get the id and length list of reference sequences
    seq_id_list = []
    seq_length_list = []
    for each_seq in SeqIO.parse(pwd_ref_file, 'fasta'):
        seq_id_list.append(each_seq.id)
        seq_length_list.append(len(each_seq.seq))

    # get the number of reads need to be simulated from each sequence
    read_number_list = []
    n = 1
    for each_length in seq_length_list:
        if n < len(seq_length_list):
            current_reads_number = round(read_number*each_length/sum(seq_length_list))
            n += 1
            read_number_list.append(current_reads_number)
        elif n == len(seq_length_list):
            current_reads_number = read_number - sum(read_number_list)
            read_number_list.append(current_reads_number)

    # prepare output file name
    path, file_name = os.path.split(pwd_ref_file)
    genome_name, ext = os.path.splitext(file_name)
    output_r1 = '%s/%s_R1.fasta' % (pwd_output_folder, genome_name)
    output_r2 = '%s/%s_R2.fasta' % (pwd_output_folder, genome_name)
    output_r12 = '%s/%s_R12.fasta' % (pwd_output_folder, genome_name)

    # create output reads file
    if split == True:
        output_r1_handle = open(output_r1, 'w')
        output_r2_handle = open(output_r2, 'w')
    else:
        output_combined_handle = open(output_r12, 'w')

    # simulate reads
    fragment_length = 2 * read_length + insert_size
    m = 0
    total_simulated = 1
    for each_seq2 in SeqIO.parse(pwd_ref_file, 'fasta'):
        reads_num_to_simulate = read_number_list[m]
        current_seq = str(each_seq2.seq)
        current_seq_len = seq_length_list[m]
        simulated = 1
        while simulated <= reads_num_to_simulate:
            random_position = random.randint(1, current_seq_len) # random_position
            current_fragment = ''

            # if simulated reads located at the middle
            if (random_position + fragment_length) <= current_seq_len:
                current_fragment = current_seq[random_position - 1: random_position + fragment_length - 1]

            # if simulated reads located at the end
            elif (random_position + fragment_length) > current_seq_len:
                if circle == True:
                    seq_part_1_seq = current_seq[random_position - 1:]
                    seq_part_2_seq = current_seq[:fragment_length - current_seq_len + random_position - 1]
                    current_fragment = seq_part_1_seq + seq_part_2_seq

            # get forward and reverse reads
            current_fragment_r1 = current_fragment[:read_length]
            current_fragment_r2 = current_fragment[-read_length:]
            current_fragment_r2_reverse_complement = str(Seq(current_fragment_r2, generic_dna).reverse_complement())
            current_read_r1_id = '%s_||_%s_r1' % (genome_name, total_simulated)
            current_read_r2_id = '%s_||_%s_r2' % (genome_name, total_simulated)

            # write out sequence
            if current_fragment != '':
                if split == True:
                    export_dna_record(current_fragment_r1, current_read_r1_id, '', output_r1_handle)
                    export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_r2_handle)
                else:
                    export_dna_record(current_fragment_r1, current_read_r1_id, '', output_combined_handle)
                    export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_combined_handle)

                simulated += 1
                total_simulated += 1
        m += 1

    # close file
    if split == True:
        output_r1_handle.close()
        output_r2_handle.close()
    else:
        output_combined_handle.close()


def get_table(list_of_ref_genomes, ctg_id_list, min_mapped_reads, in_file, out_table, out_unplotted_seqs):
    genome_1 = list_of_ref_genomes[0]
    genome_2 = list_of_ref_genomes[1]
    mapped_reads_num_dict = {}
    for each in open(in_file):
        each_split = each.strip().split('\t')
        ref_name = each_split[0]
        source_genome = each_split[1]
        mapped_reads_num = int(each_split[2])
        key = '%s___%s' % (ref_name, source_genome)
        mapped_reads_num_dict[key] = mapped_reads_num

    out_table_handle = open(out_table, 'w')
    out_unplotted_seqs_handle = open(out_unplotted_seqs, 'w')
    out_unplotted_seqs_handle.write('Sequence_ID\tMapped_reads_from_%s\tMapped_reads_from_%s\n' % (genome_1, genome_2))

    for each_ctg in ctg_id_list:
        key_genome_1 = '%s___%s' % (each_ctg, genome_1)
        key_genome_2 = '%s___%s' % (each_ctg, genome_2)
        genome_1_mapped_reads = 0
        genome_2_mapped_reads = 0

        if key_genome_1 in mapped_reads_num_dict:
            genome_1_mapped_reads = mapped_reads_num_dict[key_genome_1]
        if key_genome_2 in mapped_reads_num_dict:
            genome_2_mapped_reads = mapped_reads_num_dict[key_genome_2]

        if (genome_1_mapped_reads + genome_2_mapped_reads) < min_mapped_reads:
            out_unplotted_seqs_handle.write('%s\t%s\t%s\n' % (each_ctg, genome_1_mapped_reads, genome_2_mapped_reads))
        else:
            out_table_handle.write('%s\t%s\t%s\n' % (each_ctg, genome_1_mapped_reads, genome_2_mapped_reads))

    out_table_handle.close()
    out_unplotted_seqs_handle.close()


def get_purity_and_assignment(purity_cutoff, list_of_ref_genomes, seq_id_to_length_dict, file_in, file_out):
    genome_1 = list_of_ref_genomes[0]
    genome_2 = list_of_ref_genomes[1]
    file_out_handle = open(file_out, 'w')
    file_out_handle.write('Sequence_ID\tAssignment\tPurity(%)\tLength(Mbp)\n')
    for each in open(file_in):
        each_split = each.strip().split('\t')
        query_id = each_split[0]
        mapped_reads_genome_1 = int(each_split[1])
        mapped_reads_genome_2 = int(each_split[2])
        total_mapped_reads = mapped_reads_genome_1 + mapped_reads_genome_2
        mapped_reads_genome_1_percent = mapped_reads_genome_1 / total_mapped_reads * 100
        mapped_reads_genome_1_percent = float("{0:.2f}".format(mapped_reads_genome_1_percent))

        # get query assignment and purity for output
        query_assignment = ''
        purity_for_output = ''
        if mapped_reads_genome_1_percent >= purity_cutoff:
            query_assignment = genome_1
            purity_for_output = mapped_reads_genome_1_percent
        elif mapped_reads_genome_1_percent <= (100 - purity_cutoff):
            query_assignment = genome_2
            purity_for_output = 100 - mapped_reads_genome_1_percent
        else:
            query_assignment = 'Ambiguous'
            purity_for_output = '%s/%s' % (
            mapped_reads_genome_1_percent, float("{0:.2f}".format(100 - mapped_reads_genome_1_percent)))

        # write out
        seq_length_bp = seq_id_to_length_dict[query_id]
        seq_length_Mbp = float("{0:.2f}".format(seq_length_bp/(1024*1024)))
        file_out_handle.write('%s\t%s\t%s\t%s\n' % (query_id, query_assignment, purity_for_output, seq_length_Mbp))
    file_out_handle.close()


def get_overall_purity(query_file, pwd_out_txt_file, pwd_output_folder):
    # get sum of product
    sum_of_product = 0
    total_length_of_plotted_queries = 0
    for each_assignment in open(pwd_out_txt_file):
        if not each_assignment.startswith('Sequence_ID'):
            each_assignment_split = each_assignment.strip().split('\t')
            ref_assignment = each_assignment_split[1]
            purity_value = each_assignment_split[2]
            seq_len = float(each_assignment_split[3])
            if ref_assignment != 'Ambiguous':
                product = float(purity_value) * float(seq_len)
                sum_of_product += product
            total_length_of_plotted_queries += seq_len

    # check whether total_length_of_plotted_queries is 0
    if total_length_of_plotted_queries == 0:
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' No query sequence got enough mapped reads, please specify a higher "-n" or a lower "-m"')
        os.system('rm -r %s' % pwd_output_folder)
        exit()

    # get total length of query sequences in Mbp
    query_seq_total_length_bp = get_total_length(query_file)
    query_seq_total_length_Mbp = float("{0:.2f}".format(query_seq_total_length_bp / (1024 * 1024)))

    # calculate overall purity
    purity_overall = float("{0:.2f}".format(sum_of_product / query_seq_total_length_Mbp))
    purity_plotted = float("{0:.2f}".format(sum_of_product / total_length_of_plotted_queries))

    return purity_overall, purity_plotted


####################################################### Arguments ######################################################

parser = argparse.ArgumentParser()

parser.add_argument('-r',
                    required=True,
                    help='folder holds reference sequences')

parser.add_argument('-x',
                    required=True,
                    help='the extension of reference files')

parser.add_argument('-q',
                    required=True,
                    help='query sequences')

parser.add_argument('-n',
                    required=True,
                    type=int,
                    help='the number of reads to simulate')

parser.add_argument('-p',
                    required=False,
                    type=int,
                    default=85,
                    help='purity cut-off for reference assignment, default: 85')

parser.add_argument('-l',
                    required=False,
                    type=int,
                    default=250,
                    help='reads length, default: 250')

parser.add_argument('-i',
                    required=False,
                    type=int,
                    default=500,
                    help='reads insertion size, default: 500')

parser.add_argument('-m',
                    required=False,
                    type=int,
                    default=100,
                    help='the minimum number (NONZERO) of mapped reads required for purity calculation, default: 100')

parser.add_argument('-bbmap',
                    required=False,
                    default='bbmap.sh',
                    help='path to BBMAP executable file, default: bbmap.sh')

parser.add_argument('-circle',
                    action="store_true",
                    required=False,
                    help='specify to simulate reads crossing the breaking point for circularized reference sequences')

parser.add_argument('-split',
                    action="store_true",
                    required=False,
                    help='specify to export simulated forward and reverse reads into separate files')

parser.add_argument('-keep_temp',
                    action="store_true",
                    required=False,
                    help='specify to keep temporary files')

parser.add_argument('-tuning',
                    action="store_true",
                    required=False,
                    help='tuning mode')

parser.add_argument('-quiet',
                    action="store_true",
                    required=False,
                    help='specify to suppress reporting information')

args = vars(parser.parse_args())
reference_folder = args['r']
ref_extension = args['x']
query_file = args['q']
read_number = args['n']
percentage_cutoff = args['p']
read_length = args['l']
insert_size = args['i']
min_mapped_reads = args['m']
keep_temp = args['keep_temp']
silent = args['quiet']
pwd_bbmap = args['bbmap']
circle = args['circle']
split = args['split']
tuning_mode = args['tuning']

if reference_folder[-1] == '/':
    recipients_folder = reference_folder[:-1]

# check whether reference folder exist
if not os.path.isdir(reference_folder):
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Program exited: %s not found' % (reference_folder))
    exit()

# check whether query file exist
if not os.path.isfile(query_file):
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Program exited: %s not found' % (query_file))
    exit()

############################################### Define file name and location ##########################################

wd = os.getcwd()
query_base_name, query_ext = os.path.splitext(query_file)

temp_folder =                       '%s_temp'                           % query_base_name
image_file =                        '%s.purity.png'                     % query_base_name
out_txt_file =                      '%s.purity.txt'                     % query_base_name
out_unplotted_seqs =                '%s.purity.unplotted.txt'           % query_base_name
bbmap_report =                      '%s.bbmap_report.txt'               % (query_base_name)
pwd_output_folder =                 '%s/%s'                             % (wd, temp_folder)
pwd_bbmap_report =                  '%s/%s'                             % (wd, bbmap_report)
pwd_out_txt_file =                  '%s/%s'                             % (wd, out_txt_file)
pwd_image_file =                    '%s/%s'                             % (wd, image_file)
pwd_out_unplotted_seqs =            '%s/%s'                             % (wd, out_unplotted_seqs)
pwd_combined_reads =                '%s/%s'                             % (pwd_output_folder, 'combined_simulated_reads.fasta')
pwd_sam_file =                      '%s/%s.sam'                         % (pwd_output_folder, query_base_name)
tmp_out_file =                      '%s/%s_tmp.txt'                     % (pwd_output_folder, query_base_name)
tmp_out_file_sorted =               '%s/%s_sorted.txt'                  % (pwd_output_folder, query_base_name)
tmp_out_file_sorted_count =         '%s/%s_sorted_count.txt'            % (pwd_output_folder, query_base_name)
tmp_out_file_sorted_count_table =   '%s/%s_sorted_count_table.txt'      % (pwd_output_folder, query_base_name)

########################################################################################################################

# get reference file list
ref_file_rx = '%s/*.%s' % (reference_folder, ref_extension)
reference_file_list = [os.path.basename(file_name) for file_name in glob.glob(ref_file_rx)]
reference_file_list = sorted(reference_file_list)


# get the list of reference file without extension
reference_file_list_no_ext = []
for each_reference_file in reference_file_list:
    each_reference_file_base_name, each_reference_file_ext = os.path.splitext(each_reference_file)
    reference_file_list_no_ext.append(each_reference_file_base_name)


# get the total length of sequences in each reference file
reference_file_total_length_list = []
for each_ref in reference_file_list:
    pwd_ref_file = '%s/%s' % (reference_folder, each_ref)
    reference_file_total_length_list.append(get_total_length(pwd_ref_file))


# get the number of reads need to be simulated from each reference file (for equal depth)
ref_read_number_list = []
n = 1
for each_total_length in reference_file_total_length_list:
    if n < len(reference_file_total_length_list):
        current_reads_number = round(read_number*each_total_length/sum(reference_file_total_length_list))
        ref_read_number_list.append(current_reads_number)
    elif n == len(reference_file_total_length_list):
        current_reads_number = read_number - sum(ref_read_number_list)
        ref_read_number_list.append(current_reads_number)
    n += 1


# create temporary folder
if tuning_mode != True:
    if os.path.isdir(pwd_output_folder):
        shutil.rmtree(pwd_output_folder)
        os.makedirs(pwd_output_folder)
    else:
        os.makedirs(pwd_output_folder)


# simulate reads from each reference file
if tuning_mode != True:
    n = 0
    for each_ref2 in reference_file_list:
        if silent != True:
            print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Simulating %s reads from %s' % (ref_read_number_list[n], each_ref2))
        pwd_ref2_file = '%s/%s' % (reference_folder, each_ref2)
        simulate_reads(pwd_ref2_file, ref_read_number_list[n], read_length, insert_size, pwd_output_folder)
        n += 1


# combine simulated reads
if tuning_mode != True:
    if silent != True:
        sleep(1)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Combine all simulated reads together')
    os.system('cat %s/*.fasta > %s/combined_simulated_reads.fasta' % (pwd_output_folder, pwd_output_folder))


# run BBMAP
if tuning_mode != True:
    if silent != True:
        sleep(1)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running BBMAP')
    current_wd = os.getcwd()
    os.chdir(pwd_output_folder)
    bbmap_cmd = '%s in=%s out=%s ref=../%s ambiguous=toss &> %s' % (pwd_bbmap, pwd_combined_reads, pwd_sam_file, query_file, pwd_bbmap_report)
    os.system(bbmap_cmd)
    os.chdir(current_wd)


# get the ID list and length of reference sequences
query_ctg_id_list = []
seq_id_length_dict = {}
for each_seq in SeqIO.parse(query_file, 'fasta'):
    query_ctg_id_list.append(each_seq.id)
    seq_id_length_dict[each_seq.id] = len(each_seq.seq)


# parse mapping results
if silent != True:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Parsing mapping results')
out_handle = open(tmp_out_file, 'w')
for each_aln in open(pwd_sam_file):
    if not each_aln.startswith('@'):
        each_aln_split = each_aln.strip().split('\t')
        query_source = each_aln_split[0].split('_||_')[0]
        ref_name = each_aln_split[2].split(' ')[0]
        if ref_name != '*':
            out_handle.write('%s\t%s\n' % (ref_name, query_source))
out_handle.close()


# get the number of reads mapped to each query sequences from the two references
if silent != True:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Calculating purity')

os.system('cat %s | sort > %s' % (tmp_out_file, tmp_out_file_sorted))
out_file_sorted_count_handle = open(tmp_out_file_sorted_count, 'w')
default_line = ''
current_count = 0
for each_line in open(tmp_out_file_sorted):
    each_line_strip = each_line.strip()
    each_line_split = each_line.strip().split('\t')
    ref = each_line_split[0]
    query = each_line_split[1]
    if default_line == '':
        default_line = each_line_strip
        current_count += 1
    elif default_line == each_line_strip:
        current_count += 1
    elif default_line != each_line_strip:
        out_file_sorted_count_handle.write('%s\t%s\n' % (default_line, current_count))
        default_line = each_line_strip
        current_count = 1
out_file_sorted_count_handle.write('%s\t%s\n' % (default_line, current_count))
out_file_sorted_count_handle.close()


# transfer the number of mapped reads to table
get_table(reference_file_list_no_ext, query_ctg_id_list, min_mapped_reads, tmp_out_file_sorted_count, tmp_out_file_sorted_count_table, pwd_out_unplotted_seqs)


# write out purity and reference assignment for query sequences
get_purity_and_assignment(percentage_cutoff, reference_file_list_no_ext, seq_id_length_dict, tmp_out_file_sorted_count_table, pwd_out_txt_file)


# prepare data for plotting
if silent != True:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Preparing plot')
# read in purity and sequence length for plot
purity_list_genome_1 = []
purity_list_genome_2 = []
purity_list_Ambiguous = []
seq_length_list_genome_1 = []
seq_length_list_genome_2 = []
seq_length_list_Ambiguous = []
for each_assignment in open(pwd_out_txt_file):
    if not each_assignment.startswith('Sequence_ID'):
        each_assignment_split = each_assignment.strip().split('\t')
        ref_assignment = each_assignment_split[1]
        purity_value = each_assignment_split[2]
        seq_len = float(each_assignment_split[3])

        if ref_assignment == reference_file_list_no_ext[0]:
            purity_list_genome_1.append(100 - float(purity_value))
            seq_length_list_genome_1.append(seq_len)

        if ref_assignment == reference_file_list_no_ext[1]:
            purity_list_genome_2.append(float(purity_value))
            seq_length_list_genome_2.append(seq_len)

        if ref_assignment == 'Ambiguous':
            purity_value = float(purity_value.split('/')[0])
            purity_list_Ambiguous.append(100 - purity_value)
            seq_length_list_Ambiguous.append(seq_len)

# calculate purity for all and only plotted sequences
purity_overall, purity_plotted = get_overall_purity(query_file, pwd_out_txt_file, pwd_output_folder)

# plot
plt_2_10 = plt.scatter(purity_list_genome_1, seq_length_list_genome_1, c='g', lw=0, s=10)
plt_BS107 = plt.scatter(purity_list_genome_2, seq_length_list_genome_2, c='b', lw=0, s=10)
plt_Ambiguous = plt.scatter(purity_list_Ambiguous, seq_length_list_Ambiguous, c='r', lw=0, s=10)
plt.title('Purity: ' + str(purity_overall) + '% (overall), ' + str(purity_plotted) + '% (only plotted)', fontsize=10)
#plt.title('Overall purity: ' + str(purity_overall) + '%', fontsize=10)
plt.ylabel('Sequence length (Mbp)', fontsize=10)
plt.axis([-5, 105, 0, None])  # set the range of X-axis

# add separating lines
plt.axvline(x=0, c='black', linestyle='dashed', dashes=(5, 10), linewidth=0.3)
plt.axvline(x=100, c='black', linestyle='dashed', dashes=(5, 10), linewidth=0.3)
plt.axvline(x=(100 - percentage_cutoff), c='black', linewidth=0.3)
plt.axvline(x=percentage_cutoff, c='black', linewidth=0.5)

# name of X-axis labels
my_xticks = ['100%',
             reference_file_list_no_ext[0],
             str(percentage_cutoff) + '%',
             'Ambiguous',
             str(percentage_cutoff) + '%',
             reference_file_list_no_ext[1],
             '100%']

# position and format of X-axis label
plt.xticks(np.array([0,
                     (100 - percentage_cutoff)/2,
                     (100 - percentage_cutoff),
                     50,
                     percentage_cutoff,
                     (100 - (100 - percentage_cutoff)/2),
                     100]),
           my_xticks,
           rotation=315,
           fontsize=10,
           horizontalalignment='left')

plt.yticks(fontsize=10)  # format of Y-axis label
plt.tick_params(axis='x', bottom='off', top='off')  # hide X-axis ticks
plt.tight_layout()
plt.savefig(pwd_image_file, dpi=300)
plt.close()


if silent != True:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Purity of query sequences exported to: %s' % (out_txt_file))

if keep_temp != True:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Deleting temporary files')
    os.system('rm -r %s' % pwd_output_folder)

if silent != True:
    sleep(1)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' All done!')
