import os
from Bio import SeqIO


os.chdir('/Users/songweizhi/Desktop/test_recircularize')
end_len = 30000
contig = 'hcq44.fasta'


for each in SeqIO.parse(contig, 'fasta'):

   each_id = each.id
   each_len = len(each.seq)
   file_start = '%s_%s_%s.fasta' % (each_id, 1, end_len)
   file_end = '%s_%s_%s.fasta' % (each_id, (each_len-end_len), each_len)

   file_start_handle = open(file_start, 'w')
   file_end_handle = open(file_end, 'w')

   seq_start = str(each.seq)[:end_len]
   seq_end = str(each.seq)[each_len-end_len:]

   file_start_handle.write('>%s_%s_%s\n' % (each_id, 1, end_len))
   file_start_handle.write('%s\n' % seq_start)
   file_start_handle.close()

   file_end_handle.write('>%s_%s_%s\n' % (each_id, (each_len-end_len), each_len))
   file_end_handle.write('%s\n' % seq_end)
   file_end_handle.close()

   # run blast
   blast_output = '%s_%s_blast_outfmt6.txt' % (each_id, end_len)
   blastn_cmd = 'blastn -query %s -subject %s -outfmt 6 -out %s' % (file_start, file_end, blast_output)
   os.system(blastn_cmd)

ds = 0


