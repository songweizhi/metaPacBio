
strain_name = 'BS107'

blast_results = '/Users/songweizhi/Desktop/get_final_correlation/%s_blast_results.txt' % strain_name
blast_results_filtered = '/Users/songweizhi/Desktop/get_final_correlation/%s_blast_results_filtered.txt' % strain_name
blast_results_filtered_sum = '/Users/songweizhi/Desktop/get_final_correlation/%s_blast_results_sum.txt' % strain_name


temp_1 = open(blast_results_filtered, 'w')
for each in open(blast_results):
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    identity = float(each_split[2])
    align_len = int(each_split[3])
    query_len = int(each_split[12])
    subject_len = int(each_split[13])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))
    if (identity >= 99) and (align_len >= 1000):
        temp_1.write('%s\t%s\t%s\n' % (subject, query, align_len))
temp_1.close()


temp_2 = open(blast_results_filtered_sum, 'w')
temp_2.write('Reference,Contig,Length\n')
current_ref_query = ''
total_len = 0
for each2 in open(blast_results_filtered):
    each2_split = each2.strip().split('\t')
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




