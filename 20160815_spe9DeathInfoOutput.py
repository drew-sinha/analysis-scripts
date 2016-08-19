import analyzeHealth
import pickle
import csv


import collections

import numpy as np
import annotation_file

# Modified to NEVER LAID EGGS
# 15B 003
# 13C 68


# Load annotation files
expt_dirs = [
    '/mnt/iscopearray/Zhang_William/2016.02.16 spe-9 Run 9/',      #0
    '/mnt/iscopearray/Zhang_William/2016.02.20 spe-9 Run 10A/',    #1
    '/mnt/iscopearray/Zhang_William/2016.02.20 spe-9 Run 10B/',    #2
    '/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11A/',    #3
    '/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11B/',    #4
    '/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11C/',    #5
    '/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11D/',    #6
    '/mnt/iscopearray/Zhang_William/2016.02.29 spe-9 Run 12A/',    #7
    '/mnt/iscopearray/Zhang_William/2016.02.29 spe-9 Run 12B/',    #8
    '/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13A/',    #9
    '/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13B/',    #10
    '/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13C/',    #11
    '/mnt/iscopearray/Zhang_William/2016.03.14 spe-9 Run 14/',     #12
    '/mnt/iscopearray/Zhang_William/2016.03.25 spe-9 Run 15A/',    #13
    '/mnt/iscopearray/Zhang_William/2016.03.25 spe-9 Run 15B/',    #14
    '/mnt/iscopearray/Zhang_William/2016.03.31 spe-9 Run 16/']     #15
md_list = [
    {'W':'/mnt/iscopearray/Zhang_William/2016.02.16 spe-9 Run 9/experiment_metadata.json','D':'/mnt/bulkhelium/Sinha_Drew/20160216_spe9Acquisition/experiment_metadata.json'},   #0
    {'':''},    #1
    {'':''},    #2
    {'':''},    #3
    {'':''},    #4
    {'':''},    #5
    {'':''},    #6
    {'W':'/mnt/iscopearray/Zhang_William/2016.02.29 spe-9 Run 12A/experiment_metadata.json', 'D':'/mnt/bulkhelium/Sinha_Drew/20160229_spe9Acquisition_DevVarA/experiment_metadata.json'},    #7
    {'W':'/mnt/iscopearray/Zhang_William/2016.02.29 spe-9 Run 12B/experiment_metadata.json','D':'/mnt/bulkhelium/Sinha_Drew/20160229_spe9Acquisition_DevVarB/experiment_metadata.json'},    #8
    {'':''},    #9
    {'W':'/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13B/experiment_metadata.json', 'D':'/mnt/bulkhelium/Sinha_Drew/20160304_spe9Acquisition_DevVarB/experiment_metadata.json'},    #10
    {'W':'/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13C/experiment_metadata.json', 'D':'/mnt/bulkhelium/Sinha_Drew/20160304_spe9Acquisition_DevVarC/experiment_metadata.json'},    #11
    {'':''},    #12
    {'':''},    #13
    {'':''},    #14
    {'':''},    #15
]

ann_timestamped_data = annotation_file.compile_expt_timestamped_data(expt_dirs,md_list)
ann_raw_data = annotation_file.compile_expt_raw_data(expt_dirs)
ann_timepointed_data = annotation_file.compile_expt_timestamped_data(expt_dirs,md_list, as_timepoints=True)

def save_dict_as_tsv(my_dict, output_file):
    with open(output_file,'w') as output_fp:
        output_writer = csv.writer(output_fp,delimiter='\t')
        output_writer.writerow([a_tag for a_tag in my_dict.keys()])
        for worm_data in zip(*(my_dict.values())):
            output_writer.writerow(worm_data)


# Load df
strain='spe-9'
with open('/mnt/bulkdata/wzhang/human_dir/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
    strain_df = pickle.load(my_file)['adult_df']
print(strain+": n = {}".format(len(strain_df.worms)))

# Build new data structure
# Worm name DONE
# Lifespan (real + reinterpolated?) DONE, DONE
# Frame of death (real) & timepoint: DONE, NOT DONE
# Dying health DONE
# Which cohort? Won't do now.


abb_data = collections.OrderedDict()
abb_data['Worm'] = ann_timestamped_data['Worm_FullName']
abb_data['AnnotatedDeathFrame'] = [death_value[1:] if death_value[0].isalpha() else death_value for death_value in ann_raw_data['Death']]   # Specify to use WZ's file
abb_data['AnnotatedDeathTimepoint'] = ann_timepointed_data['Death']
abb_data['AnnotatedLifespan'] = (ann_timestamped_data['Death']-ann_timestamped_data['Hatch'])/(3600*24) # Days
abb_data['InterpolatedLifespan'] = analyzeHealth.selectData.get_lifespans(strain_df)/24
abb_data['Dying Health'] = []
health_idx=-4
for i,worm_data in enumerate(strain_df.mloc(measures=['ghost_age','health'])):
    last_alive = np.nonzero((~np.isnan(worm_data[0])) & (~np.isnan(worm_data[1])))[0][-1]
    abb_data['Dying Health'].append(strain_df.display_variables(worm_data[1][last_alive],'health')[0])

save_dict_as_tsv(abb_data, '/media/Data/Work/ZPLab/temp_Emily.tsv')
