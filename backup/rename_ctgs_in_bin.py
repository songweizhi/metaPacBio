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


os.chdir('/Users/songweizhi/Desktop')

subreads_files = '/Users/songweizhi/Desktop/qualified_bins_noNs_lt2000/*.fa'
file_list = [os.path.basename(file_name) for file_name in glob.glob(subreads_files)]

for file in file_list:
    file_name_new = '%s%s.fasta' % (file.split('.')[0], file.split('.')[1])

    n = 1
    renamed_file = open('qualified_bins_noNs_lt2000_renamed/%s' % file_name_new, 'w')
    for each_seq in SeqIO.parse('qualified_bins_noNs_lt2000/%s' % file, 'fasta'):
        each_seq_new_id = '%s%s_ctg%s' % (file.split('.')[0], file.split('.')[1], n)
        export_dna_record(str(each_seq.seq), each_seq_new_id, '', renamed_file)
        n += 1
    renamed_file.close()
