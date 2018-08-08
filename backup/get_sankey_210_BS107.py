import os


# input
#os.chdir('/Users/songweizhi/Desktop')
os.chdir('/Users/songweizhi/Desktop/2018-08-07_Sankey_210_BS107')


#blast_result = 'blast_results_AD91A_vs_refs.txt'
blast_result = 'blast_results.txt'

align_len_cutoff = 5000


# define tmp file name
tmp1 = 'output_tmp1.txt'
tmp2 = 'output_tmp2.txt'
sankey_input = 'sankey_input.txt'


# for the nine
out_tmp_handle = open(tmp1, 'w')
for each in open(blast_result):
    print(each)
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    identity = float(each_split[2])
    alignment_len = int(each_split[3])
    for_out = '%s,%s,%s\n' % (subject, query, alignment_len)

    if (identity >= 99) and (alignment_len >= align_len_cutoff):
        out_tmp_handle.write(for_out)

out_tmp_handle.close()


# sort results
os.system('cat %s | sort > %s' % (tmp1, tmp2))



temp_2 = open(sankey_input, 'w')
temp_2.write('Reference,Contig,Length\n')
current_ref_query = ''
total_len = 0
for each2 in open(tmp2):
    each2_split = each2.strip().split(',')
    ref = each2_split[0]
    query = each2_split[1]
    length = int(each2_split[2])
    ref_query = '%s___%s' % (ref, query)
    if current_ref_query == '':
        current_ref_query = ref_query
        total_len += length
    elif current_ref_query == ref_query:
        total_len += length
    elif current_ref_query != ref_query:
        temp_2.write('%s,%s,%s\n' % (current_ref_query.split('___')[0], current_ref_query.split('___')[1], total_len))
        current_ref_query = ref_query
        total_len = length
temp_2.write('%s,%s,%s\n' % (current_ref_query.split('___')[0], current_ref_query.split('___')[1], total_len))
temp_2.close()


# # remove tmp
# os.remove(tmp1)
# os.remove(tmp2)


# call R
print('Rscript ~/PycharmProjects/metaPacBio/get_sankey_plot.R -f %s -x 600 -y 600' % sankey_input)

os.system('Rscript ~/PycharmProjects/metaPacBio/get_sankey_plot.R -f %s -x 600 -y 600' % sankey_input)

