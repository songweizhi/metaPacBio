import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '02:59:00'
email = 'wythe1987@163.com'
modules_needed = ['']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_gzip_subreads'

###########################################################################################

os.chdir(wd)

# create outputs folder
if not os.path.exists(outputs_folder):
    os.makedirs(outputs_folder)
else:
    shutil.rmtree(outputs_folder)
    os.makedirs(outputs_folder)

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

################################################################

for each in open('files_for_gzip.txt'):
    each = each.strip()
    wd_on_katana = '/'.join(each.split('/')[:-1])
    file = each.split('/')[-1]
    output_handle = open('%s/qsub_gzip_%s.sh' % (outputs_folder, file), 'w')
    output_handle.write(header)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('gzip %s\n' % file)
    output_handle.close()




