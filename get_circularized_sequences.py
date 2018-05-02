import os
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


os.chdir('/Users/songweizhi/Desktop/test_recircularize')
# output = open('210_and_BS107_canu.contigs.circularized.quiver1.fasta', 'w')
# for each_seq in SeqIO.parse('210_and_BS107_canu.contigs.circularized.trimmed.quiver1.fasta', 'fasta'):
#     print('%s\t%s' % (each_seq.id, len(each_seq.seq)))
#     print(each_seq.seq)




for each in open('blast_results_BS107_quiver2.txt'):
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    qstart = each_split[6]
    qend = each_split[7]
    sstart = each_split[8]
    send = each_split[9]
    align_length = int(each_split[3])
    iden = float(each_split[2])
    if (iden >= 99) and (align_length >= 10000):
        print(each.strip())



