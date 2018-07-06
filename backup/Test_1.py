import os
from Bio import SeqIO

os.chdir('/Users/songweizhi/Desktop/555')

for each in SeqIO.parse('/Users/songweizhi/Desktop/555/SON1053.SP16554.hcq.qv20.fas', 'fasta'):

    file_name = each.id.split('_')[0] + '.fasta'

    file_handle = open(file_name, 'w')

    print(each.id)
    SeqIO.write(each, file_handle, 'fasta')
    file_handle.close()





