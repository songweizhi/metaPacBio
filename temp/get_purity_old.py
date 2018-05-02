#!/usr/bin/env python3

import os
import argparse
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

usage = '''

module load python/3.5.2
module load bbmap/35.82
python3 /srv/scratch/z5039045/Scripts/get_purity.py -ref cns_h_ctg.fasta -reads simulated_reads.fasta -c 75

'''

parser = argparse.ArgumentParser()

parser.add_argument('-ref',
                    required=True,
                    help='reference sequences')

parser.add_argument('-reads',
                    required=True,
                    help='reads file')

parser.add_argument('-c',
                    required=False,
                    default=75,
                    type=int,
                    help='percentage_cutoff')

args = vars(parser.parse_args())
ref_file = args['ref']
reads_file = args['reads']
percentage_cutoff = int(args['c'])
basename, ext = os.path.splitext(ref_file)
sam_file = '%s.sam' % basename
image_file = '%s.png' % basename
out_txt_file = '%s.txt' % basename

# tmp files
tmp_out_file = '%s_tmp.txt' % basename
tmp_out_file_sorted = '%s_sorted.txt' % basename
tmp_out_file_sorted_count = '%s_sorted_count.txt' % basename
tmp_out_file_sorted_count_table = '%s_sorted_count_table.txt' % basename


# run BBMAP
os.system('/Users/songweizhi/Softwares/bbmap/bbmap.sh in=%s out=%s ref=%s ambiguous=toss' % (reads_file, sam_file, ref_file))


# get length of reference sequences
seq_id_length_dict = {}
for each_seq in SeqIO.parse(ref_file, 'fasta'):
    seq_id_length_dict[each_seq.id] = len(each_seq.seq)

# parse sam file
out_handle = open(tmp_out_file, 'w')
for each_aln in open(sam_file):
    if not each_aln.startswith('@'):
        each_aln_split = each_aln.strip().split('\t')
        query_name = each_aln_split[0]
        query_source = query_name.split('_')[0]
        ref_name = each_aln_split[2].split(' ')[0]
        query_seq = each_aln_split[9]
        if ref_name != '*':
            out_handle.write('%s\t%s\n' % (ref_name, query_source))
out_handle.close()


os .system('cat %s | sort > %s' % (tmp_out_file, tmp_out_file_sorted))
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


out_file_sorted_count_table_handle = open(tmp_out_file_sorted_count_table, 'w')
current_ref = ''
current_out = []
for each in open(tmp_out_file_sorted_count):
    each_split = each.strip().split('\t')
    ref_name = each_split[0]
    source_genome = each_split[1]
    mapped_reads_num = str(each_split[2])
    if current_ref == '':
        current_ref = ref_name
        current_out = [source_genome, mapped_reads_num]
    elif current_ref == ref_name:
        current_out.append(source_genome)
        current_out.append(mapped_reads_num)
    elif current_ref != ref_name:
        if len(current_out) == 4:
            out_file_sorted_count_table_handle.write('%s\t%s\t%s\n' % (current_ref, current_out[1], current_out[3]))
        elif len(current_out) == 2:
            if current_out[0] == '2.10':
                out_file_sorted_count_table_handle.write('%s\t%s\t%s\n' % (current_ref, current_out[1], '0'))
            elif current_out[0] == 'BS107':
                out_file_sorted_count_table_handle.write('%s\t%s\t%s\n' % (current_ref, '0', current_out[1]))
        current_ref = ref_name
        current_out = []
        current_out = [source_genome, mapped_reads_num]


if len(current_out) == 4:
    out_file_sorted_count_table_handle.write('%s\t%s\t%s\n' % (current_ref, current_out[1], current_out[3]))
elif len(current_out) == 2:
    if current_out[0] == '2.10':
        out_file_sorted_count_table_handle.write('%s\t%s\t%s\n' % (current_ref, current_out[1], '0'))
    elif current_out[0] == 'BS107':
        out_file_sorted_count_table_handle.write('%s\t%s\t%s\n' % (current_ref, '0', current_out[1]))
out_file_sorted_count_table_handle.close()


cotig_id_list_2_10 = []
cotig_id_list_BS107 = []
cotig_id_list_Ambiguous = []
purity_list_2_10 = []
purity_list_BS107 = []
purity_list_Ambiguous = []
length_list_2_10 = []
length_list_BS107 = []
length_list_Ambiguous = []
out_txt_file_handle = open(out_txt_file, 'w')
out_txt_file_handle.write('Ctg_ID\tAssignment\tPurity(%)\tLength(bp)\n')

for each_ref in open(tmp_out_file_sorted_count_table):
    each_ref_split = each_ref.strip().split('\t')
    each_ref_name = each_ref_split[0]
    ref_seq_length = seq_id_length_dict[each_ref_name]
    reads_from_2_10 = int(each_ref_split[1])
    reads_from_BS107 = int(each_ref_split[2])
    reads_from_2_10_percent = round(reads_from_2_10 *100/(reads_from_2_10 + reads_from_BS107), 2)

    if reads_from_2_10_percent >= percentage_cutoff:
        cotig_id_list_2_10.append(each_ref_name)
        purity_list_2_10.append(reads_from_2_10_percent)
        length_list_2_10.append(ref_seq_length)
        out_txt_file_handle.write('%s\t%s\t%s\t%s\n' % (each_ref_name, '2.10', reads_from_2_10_percent,ref_seq_length))

    elif reads_from_2_10_percent <= (100 - percentage_cutoff):
        cotig_id_list_BS107.append(each_ref_name)
        purity_list_BS107.append(reads_from_2_10_percent)
        length_list_BS107.append(ref_seq_length)
        out_txt_file_handle.write('%s\t%s\t%s\t%s\n' % (each_ref_name, 'BS107', 100 - reads_from_2_10_percent, ref_seq_length))

    elif (reads_from_2_10_percent < percentage_cutoff) and (reads_from_2_10_percent > (100 - percentage_cutoff)):
        cotig_id_list_Ambiguous.append(each_ref_name)
        purity_list_Ambiguous.append(reads_from_2_10_percent)
        length_list_Ambiguous.append(ref_seq_length)
        out_txt_file_handle.write('%s\t%s\t%s(from 2.10)\t%s\n' % (each_ref_name, 'Ambiguous', reads_from_2_10_percent, ref_seq_length))
out_txt_file_handle.close()


# calculate overall purity
purity = (sum(length_list_2_10) + sum(length_list_BS107))/(sum(length_list_2_10) + sum(length_list_BS107) + sum(length_list_Ambiguous))
purity = round(purity, 3)

plt_2_10 = plt.scatter(purity_list_2_10,length_list_2_10, c='g', lw = 0)
plt_BS107 = plt.scatter(purity_list_BS107,length_list_BS107, c='b', lw = 0)
plt_Ambiguous = plt.scatter(purity_list_Ambiguous,length_list_Ambiguous, c='r', lw = 0)
plt.title('Overall Purity: ' + str(round(purity*100, 1)) + '%')
plt.ylabel('Contig length (bp)')
plt.axis([0, 100, 0, None])
plt.axvline(x=(100 - percentage_cutoff), c='black')
plt.axvline(x=percentage_cutoff, c='black')
my_xticks = ['100%', 'BS107', str(percentage_cutoff) + '%', 'Ambiguous', str(percentage_cutoff) + '%', '2.10WT', '100%']
plt.xticks(np.array([0, (100 - percentage_cutoff)/2, (100 - percentage_cutoff), 50, percentage_cutoff, (100 - (100 - percentage_cutoff)/2), 100]), my_xticks)
plt.savefig(image_file, dpi = 300)
plt.close()

# remove tmp files
os.system('rm %s' % tmp_out_file)
os.system('rm %s' % tmp_out_file_sorted)
os.system('rm %s' % tmp_out_file_sorted_count)
os.system('rm %s' % tmp_out_file_sorted_count_table)
