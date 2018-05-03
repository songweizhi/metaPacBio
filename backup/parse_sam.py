
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


os.chdir('/Users/songweizhi/Desktop/check_sam')

for each in open('BS107_tig00000042_canu.sam'):
    if not each.startswith('@'):
        #print(each)
        each_split = each.strip().split('\t')
        qname = each_split[0]
        segment_seq = each_split[9]
        pos = int(each_split[3])
        mapQ = int(each_split[4])
        Tlen = int(each_split[8])
        #print('%s\t%s' % (pos, Tlen))
        #print('%s\t%s' % (qname, len(segment_seq)))
        if (pos <= 18230) and ((pos + len(segment_seq)) >= 19350):
            print(each.strip())



