import os
import glob
import shutil


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 20
walltime_needed = '02:59:00'
email = 'wythe1987@163.com'
modules_needed = ['python/3.5.2']

wd = '/Users/songweizhi/Desktop'
wd_on_katana = '/srv/scratch/z5039045/PacBio/2018-07-04_assess_depth'


#outputs_folder = 'qsub_subsample_10x'
#outputs_folder = 'qsub_subsample_25x'
#outputs_folder = 'qsub_subsample_50x'
outputs_folder = 'qsub_subsample_100x'

in_file_folder = '/srv/scratch/z5039045/PacBio/2018-07-04_assess_depth/subreads_all'
#out_file_folder = '/srv/scratch/z5039045/PacBio/2018-07-04_assess_depth/subreads_10x'
#out_file_folder = '/srv/scratch/z5039045/PacBio/2018-07-04_assess_depth/subreads_25x'
#out_file_folder = '/srv/scratch/z5039045/PacBio/2018-07-04_assess_depth/subreads_50x'
out_file_folder = '/srv/scratch/z5039045/PacBio/2018-07-04_assess_depth/subreads_100x'

#percent = 0.0595  # 10x
#percent = 0.149  # 25x
#percent = 0.298  # 50x
percent = 0.595  # 100x


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

for each in open('/Users/songweizhi/Desktop/subreads_file.txt'):
    print(each)
    each = each.strip()
    pwd_in = '%s/%s' % (in_file_folder, each)
    pwd_out = '%s/%s' % (out_file_folder, each)
    cmd = 'python3 /srv/scratch/z5039045/Scripts/subsample_unpaired_reads.py -i %s -p %s -o %s\n' % (pwd_in, percent, pwd_out)
    output_handle = open('%s/qsub_subsample_%s.sh' % (outputs_folder, each), 'w')
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write(cmd)
    output_handle.close()

