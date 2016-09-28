import annotation_file
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.stats

import plotting_tools

expt_dirs = [
    '/mnt/scopearray/Sinha_Drew/20160823_spe9-HT115GFP/',
    '/mnt/scopearray/Sinha_Drew/20160823_N2-HT115pos1g/'
]

annotation_fps =[
    '/mnt/scopearray/Sinha_Drew/20160823_spe9-HT115GFP/2016.08.23 spe-9+HT115-GFP - Annotations.tsv',
    '/mnt/scopearray/Sinha_Drew/20160823_N2-HT115pos1g/2016.08.23 N2+HT115-pos1_RNAi - Annotations.tsv'
]
make_labels = False
out_dir = '/media/Data/Work/ZPLab/Analysis/pos-1Diagnostics/for_presentation/'
if make_labels and 'labeled' not in out_dir: 
    out_dir = out_dir+'labeled/'
    if not os.path.isdir(out_dir): os.mkdir(out_dir)

my_ann_files = [annotation_file.AnnotationFile(a_ann_fp) for a_ann_fp in annotation_fps]
my_timestamped_data = [ann_file.data_as_timestamps_simple(expt_dir+os.path.sep+'experiment_metadata.json',
    restricted_list=np.nonzero([((hatch_str is not '') and ('NOT DEAD' in note_str) and ('DELAYED' not in note_str)) for hatch_str,note_str in zip(ann_file.data['Hatch'],ann_file.data['Notes'])])[0]) 
    for ann_file,expt_dir in zip(my_ann_files,expt_dirs)]

# Calculate relevant spans in the two populations
pop_L14_durations = [(timestamped_data['L4 ecdysis']-timestamped_data['L1 ecdysis'])/3600 for timestamped_data in my_timestamped_data]
pop_L1_durations = [(timestamped_data['L1 ecdysis']-timestamped_data['Hatch'])[timestamped_data['Hatch']>0]/3600 for timestamped_data in my_timestamped_data]

# Plot distribution of times spent in relevant portions of development
L14_min = min([min(L14_dur) for L14_dur in pop_L14_durations])
L14_max = max([max(L14_dur) for L14_dur in pop_L14_durations])
L1_min = min([min(L1_dur) for L1_dur in pop_L1_durations])
L1_max = max([max(L1_dur) for L1_dur in pop_L1_durations])


L14_bins = np.arange(start=L14_min-1,stop=L14_max+1+1)
L1_bins = np.arange(start=L1_min-1,stop=L1_max+1+1)


fig_h, ax_h = plt.subplots(2,1)
plot_data = [ax_h[0].hist(L14_dur, L14_bins) for L14_dur in pop_L14_durations]
if make_labels: ax_h[0].legend([pop_plot[2][0] for pop_plot in plot_data],['spe-9/HT115-GFP','N2/HT115-pos1'])
ax_h[0].set_title('Duration of L2-Adulthood; (N={},{} (spe-9/N2))'.format(*[len(L14_dur) for L14_dur in pop_L14_durations]))
ranksum_data = scipy.stats.mannwhitneyu(pop_L14_durations[0],pop_L14_durations[1])
ttest_data = scipy.stats.ttest_ind(pop_L14_durations[0],pop_L14_durations[1])
print('L2-A duration Median {:.2f}/{:.2f}'.format(np.median(pop_L14_durations[0]),np.median(pop_L14_durations[1])))
print('L2-A duration Mean {:.2f}/{:.2f}'.format(np.mean(pop_L14_durations[0]),np.mean(pop_L14_durations[1])))
print('L2-A duration Std {:.2f}/{:.2f}'.format(np.std(pop_L14_durations[0]),np.std(pop_L14_durations[1])))
print('L2-A duration Rank-sum p = {}'.format(ranksum_data[1]*2))
print('L2-A duration t-test p = {}'.format(ttest_data[1]))

plot_data = [ax_h[1].hist(L1_dur, L1_bins) for L1_dur in pop_L1_durations]
if make_labels: ax_h[1].legend([pop_plot[2][0] for pop_plot in plot_data],['spe-9/HT115-GFP','N2/HT115-pos1'])
ax_h[1].set_title('Duration of L1; (N={},{} (spe-9/N2))'.format(*[len(L1_dur) for L1_dur in pop_L1_durations]))
ranksum_data = scipy.stats.mannwhitneyu(pop_L1_durations[0],pop_L1_durations[1])
ttest_data = scipy.stats.ttest_ind(pop_L1_durations[0],pop_L1_durations[1])
print('L1 duration Median {:.2f}/{:.2f}'.format(np.median(pop_L1_durations[0]),np.median(pop_L1_durations[1])))
print('L1 duration Mean {:.2f}/{:.2f}'.format(np.mean(pop_L1_durations[0]),np.mean(pop_L1_durations[1])))
print('L1 duration Std {:.2f}/{:.2f}'.format(np.std(pop_L1_durations[0]),np.std(pop_L1_durations[1])))
print('L1 duration Rank-sum p = {}'.format(ranksum_data[1]*2))
print('L1 duration t-test p = {}'.format(ttest_data[1]))

[plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h]
if len(out_dir) is not 0:
    fig_h.savefig(out_dir+'pos-1_timing.png')
    
