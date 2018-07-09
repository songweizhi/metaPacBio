import os
from Bio import SeqIO

os.chdir('/Users/songweizhi/Desktop/888')

for each in SeqIO.parse('/Users/songweizhi/Desktop/888/AD91A_canu.contigs.circularized.quiver1.recircularized.fasta', 'fasta'):

    file_name = each.id.split('_')[0] + '.fasta'

    file_handle = open(file_name, 'w')

    print(each.id)
    SeqIO.write(each, file_handle, 'fasta')
    file_handle.close()





