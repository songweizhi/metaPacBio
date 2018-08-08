import os
from Bio import SeqIO


wd = '/Users/songweizhi/Desktop/test_recircularize/LSS9'

os.chdir(wd)

out_handle = open('LSS9_new.fasta', 'w')
for each in SeqIO.parse('LSS9.fasta', 'fasta'):
    SeqIO.write(each, out_handle, 'fasta')

out_handle.close()

