from Bio import SeqIO


seq_file = '/Users/songweizhi/Desktop/uuu/LSS9.faa'


total_num = 0
hypothetical_protein_num = 0
rRNA_num = 0

for each in SeqIO.parse(seq_file, 'fasta'):
    description = each.description

    #print(description)

    if 'hypothetical protein' in description:
        hypothetical_protein_num += 1

    total_num += 1


print('total_num: %s' % total_num)
print('hypothetical_protein_num: %s' % hypothetical_protein_num)











