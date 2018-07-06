
from Bio import SeqIO

output = open('/Users/songweizhi/Desktop/111/hcq44.fasta', 'w')
for each in SeqIO.parse('/Users/songweizhi/Desktop/SON1053.SP16554.hcq.qv20.fas', 'fasta'):
    print(each.id)
    if 'hcq44_' in each.id:
        SeqIO.write(each, output, 'fasta')
output.close()






