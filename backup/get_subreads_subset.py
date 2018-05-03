
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


os.chdir('/Users/songweizhi/Desktop/subreads_subset')


# get 210 subreads list
subreads_list_210 = []
for subread_210 in SeqIO.parse('210_subreads_renamed.fasta', 'fasta'):
    subread_210_id_modified = '/'.join(subread_210.id.split('/')[:-1])
    if subread_210.id not in subreads_list_210:
        subreads_list_210.append(subread_210.id)

subreads_list_BS107 = []
for subread_BS107 in SeqIO.parse('BS107_subreads_renamed.fasta', 'fasta'):
    subread_BS107_id_modified = '/'.join(subread_BS107.id.split('/')[:-1])
    if subread_BS107.id not in subreads_list_BS107:
        subreads_list_BS107.append(subread_BS107.id)


print(len(subreads_list_210))
print(len(subreads_list_BS107))





subreads_files = '/Users/songweizhi/Desktop/subreads_subset/all/*.fasta'
file_list = [os.path.basename(file_name) for file_name in glob.glob(subreads_files)]
#print(file_list)
#print(len(file_list))


m = 0
n = 0
for each_file in file_list:
    out_210 = open('210/%s' % each_file, 'w')
    out_BS107 = open('BS107/%s' % each_file, 'w')

    for subread in SeqIO.parse('all/%s' % each_file, 'fasta'):
        #print(subread.id)
        #print(subread.description)

        if subread.id in subreads_list_210:
            out_210.write('>' + subread.description + '\n')
            out_210.write(str(subread.seq) + '\n')
            m += 1

        if subread.id in subreads_list_BS107:
            out_BS107.write('>' + subread.description + '\n')
            out_BS107.write(str(subread.seq) + '\n')
            n += 1

    out_210.close()
    out_BS107.close()

print(m)
print(n)


# for subread in SeqIO.parse('all/m151114_230846_42272_c100936552550000001823207905251627_s1_p0.1.subreads.fasta', 'fasta'):
#     if subread.id in subreads_list_210:
#         print(subread.id)


