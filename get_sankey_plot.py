import os


# input
os.chdir('/Users/songweizhi/Desktop')
blast_result = 'blast_results_AD91A_vs_refs.txt'
align_len_cutoff = 5000


# define tmp file name
tmp1 = 'output_tmp1.txt'
tmp2 = 'output_tmp2.txt'
sankey_input = 'sankey_input.txt'


# for the nine
out_tmp_handle = open(tmp1, 'w')
for each in open(blast_result):
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    query_genome = query.split('_')[0]
    subject_genome = subject.split('_')[0]
    identity = float(each_split[2])
    alignment_len = int(each_split[3])
    for_out = '%s,%s,%s\n' % (subject_genome, query_genome, alignment_len)
    if query_genome != 'hcq7':
        if (identity >= 99) and (alignment_len >= align_len_cutoff):
            out_tmp_handle.write(for_out)
    if subject_genome == 'PM05':
        if (identity >= 80) and (alignment_len >= align_len_cutoff):
            out_tmp_handle.write(for_out)
out_tmp_handle.close()


# sort results
os.system('cat %s | sort > %s' % (tmp1, tmp2))


# get all match pairs
all_match_pairs = []
for each_hit in open(tmp2):
    each_hit_split = each_hit.strip().split(',')
    match_pair_id = '%s___%s' % (each_hit_split[0], each_hit_split[1])
    if match_pair_id not in all_match_pairs:
        all_match_pairs.append(match_pair_id)


# initial length dict
length_dict = {}
for each_pair in all_match_pairs:
    length_dict[each_pair] = 0


# get total length
for each_pair2 in open(tmp2):
    each_hit_split2 = each_pair2.strip().split(',')
    match_pair_id2 = '%s___%s' % (each_hit_split2[0], each_hit_split2[1])
    length = int(each_hit_split2[2])
    length_dict[match_pair_id2] += length


# set plot order
order_list = ['2.10', 'BS107', 'PM05', 'LSS9', 'AD1', 'AD10', 'AU392', 'BL5', 'BL110', 'D2']


# write out length
sankey_input_handle = open(sankey_input, 'w')
sankey_input_handle.write('Reference,PacBio,Total_alignment_length\n')
for each_ref in order_list:
    for each_key in length_dict:
        each_key_l = each_key.split('___')[0]
        each_key_r = each_key.split('___')[1]
        if each_key_l == each_ref:
            for_out = '%s,%s,%s\n' % (each_key_l, each_key_r, float("{0:.2f}".format(length_dict[each_key] / (1024 * 1024))))
            sankey_input_handle.write(for_out)
sankey_input_handle.close()


# remove tmp
os.remove(tmp1)
os.remove(tmp2)


# call R
os.system('Rscript ~/PycharmProjects/metaPacBio/get_sankey_plot.R -f %s -x 600 -y 1300' % sankey_input)

