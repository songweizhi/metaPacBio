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

out = open('/Users/songweizhi/Desktop/test/BS107_gap_sequence.fasta', 'w')
out_2000 = open('/Users/songweizhi/Desktop/test/BS107_gap_sequence_2000.fasta', 'w')
for each in SeqIO.parse('/Users/songweizhi/Desktop/test/DSM17395_chromosome.fasta', 'fasta'):
    print(each.id)
    print(each.seq)
    print(len(each.seq))
    gap_seq = str(each.seq)[1900540:1939441]
    gap_seq_2000 = str(each.seq)[1898540:1941441]
    print(gap_seq)
    print(len(gap_seq))
    print(len(gap_seq_2000))
    export_dna_record(gap_seq, 'BS107_gap_sequence', '', out)
    export_dna_record(gap_seq_2000, 'BS107_gap_sequence_2000', '', out_2000)

out.close()
out_2000.close()
