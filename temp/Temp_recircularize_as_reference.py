
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

G210_c = []
G210_p1 = []
G210_p2 = []
G210_p3 = []

BS107_c = []
BS107_pPGA1_65 = []
BS107_pPGA1_78 = []
BS107_pPGA1_262 = []


for each in open('break_point.txt'):
    #print(each)
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    qstart = each_split[6]
    qend = each_split[7]
    sstart = each_split[8]
    send = each_split[9]
    align_length = int(each_split[3])
    iden = float(each_split[2])
    if subject == '2.10_Egal_chromosome':
        G210_c.append(qstart)
        G210_c.append(qend)

    elif subject == '2.10_Egal_plasmid1':
        G210_p1.append(qstart)
        G210_p1.append(qend)

    elif subject == '2.10_Egal_plasmid2':
        G210_p2.append(qstart)
        G210_p2.append(qend)

    elif subject == '2.10_Egal_plasmid3':
        G210_p3.append(qstart)
        G210_p3.append(qend)

    elif subject == 'DSM17395_chromosome':
        BS107_c.append(qstart)
        BS107_c.append(qend)

    elif subject == 'DSM17395_pPGA1_65':
        BS107_pPGA1_65.append(qstart)
        BS107_pPGA1_65.append(qend)

    elif subject == 'DSM17395_pPGA1_78':
        BS107_pPGA1_78.append(qstart)
        BS107_pPGA1_78.append(qend)

    elif subject == 'DSM17395_pPGA1_262':
        BS107_pPGA1_262.append(qstart)
        BS107_pPGA1_262.append(qend)


G210_c = sorted(G210_c)
G210_p1 = sorted(G210_p1)
G210_p2 = sorted(G210_p2)
G210_p3 = sorted(G210_p3)
BS107_c = sorted(BS107_c)
BS107_pPGA1_65 = sorted(BS107_pPGA1_65)
BS107_pPGA1_78 = sorted(BS107_pPGA1_78)
BS107_pPGA1_262 = sorted(BS107_pPGA1_262)


print(G210_c)
print(G210_p1)
print(G210_p2)
print(G210_p3)
print(BS107_c)
print(BS107_pPGA1_65)
print(BS107_pPGA1_78)
print(BS107_pPGA1_262)








