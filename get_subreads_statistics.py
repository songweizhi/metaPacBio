
import os
import glob
from Bio import SeqIO

os.chdir('/Users/songweizhi/Desktop/subreads_subset')


subreads_files = '/Users/songweizhi/Desktop/subreads_subset/PacBio_genomes/*.fas'
file_list = [os.path.basename(file_name) for file_name in glob.glob(subreads_files)]

total_len = 0
total_num = 0
for each_file in file_list:
    for subread in SeqIO.parse('PacBio_genomes/%s' % each_file, 'fasta'):
        total_len += len(subread.seq)
        total_num += 1
average_len = total_len/total_num

print(total_len)
print(total_len/(1024*1024))
print(total_num)
print(average_len)






