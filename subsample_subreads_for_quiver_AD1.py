import os
import glob
from Bio import SeqIO


os.chdir('/Users/songweizhi/Desktop/subreads_subset')
subreads_files = '/Users/songweizhi/Desktop/subreads_subset/all/*.fasta'


# get 210 subreads list
subreads_list_AD1 = []
for subread_AD1 in SeqIO.parse('AD1_subreads_subset.fasta', 'fasta'):
    if subread_AD1.id not in subreads_list_AD1:
        subreads_list_AD1.append(subread_AD1.id)

file_list = [os.path.basename(file_name) for file_name in glob.glob(subreads_files)]

for each_file in file_list:
    out_AD1 = open('AD1/%s' % each_file, 'w')
    for subread in SeqIO.parse('all/%s' % each_file, 'fasta'):
        if subread.id in subreads_list_AD1:
            out_AD1.write('>' + subread.description + '\n')
            out_AD1.write(str(subread.seq) + '\n')
    out_AD1.close()
