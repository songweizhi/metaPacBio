

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


os.chdir('/Users/songweizhi/Desktop/rename_reads')

subreads_210_renamed = open('210_subreads_renamed.fasta', 'w')
for each210 in open('210_subreads.fasta'):
    if each210.startswith('>'):
        each210_new = '/'.join(each210.strip().split('/')[:-1])
        subreads_210_renamed.write(each210_new + '\n')
    else:
        subreads_210_renamed.write(each210)
subreads_210_renamed.close()



