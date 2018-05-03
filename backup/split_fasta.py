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


for each in SeqIO.parse('/Users/songweizhi/Desktop/test_recircularize/BS107/BS107_canu.contigs.circularized.trimmed.quiver1.fasta', 'fasta'):

    output_handle = open('/Users/songweizhi/Desktop/test_recircularize/BS107/%s.fasta' % each.id, 'w')
    export_dna_record(str(each.seq), each.description, '', output_handle)

    output_handle.close()
