
genome_list = ['210', 'AD1', 'AD10', 'AU392', 'BL110', 'BL5', 'BS107', 'D2', 'LSS9']


for each in genome_list:
    print('prokka --force --cpus 6 --prefix %s --locustag %s --strain %s --outdir %s %s.fasta &' % (each, each, each, each, each))