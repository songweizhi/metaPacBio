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


os.chdir('/Users/songweizhi/Desktop/test_recircularize_hcq7')
overlap_length = 3000

output_file = open('210_DSM17395_combined_ref_with_overlap.fasta', 'w')
for each in SeqIO.parse('combined_10.consensus_unitig7.fasta', 'fasta'):
    start_seq = each.seq[:(overlap_length)]
    seq_with_overlap = str(each.seq) + str(start_seq)
    export_dna_record(seq_with_overlap, each.id, '', output_file)
output_file.close()
