import numpy as np
import matplotlib.pyplot as plt
import os

import annotation_file

def quick_plot_dev(ann_fps, expt_mds,bad_worm_kws=[]):
    my_ann_files = [annotation_file.AnnotationFile(ann_fp) for ann_fp in ann_fps]  
    
    timestamped_data = {}
    [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    for [expt_md_fp, ann_file] in zip(expt_mds, my_ann_files):
        ann_file_data = ann_file.data_as_timestamps_simple(expt_md_fp)
        for expt_key in timestamped_data.keys():
            timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
    print(timestamped_data)
    
    viable_worm = (timestamped_data['Hatch']!=-1) \
        & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    hatched_on_corral = timestamped_data['Hatch'] != 0
    
    L1_durations = (timestamped_data['L1 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    L2_durations = (timestamped_data['L2 ecdysis']-timestamped_data['L1 ecdysis'])[viable_worm]/3600
    L3_durations = (timestamped_data['L3 ecdysis']-timestamped_data['L2 ecdysis'])[viable_worm]/3600
    L4_durations = (timestamped_data['L4 ecdysis']-timestamped_data['L3 ecdysis'])[viable_worm]/3600
    larval_durations = (timestamped_data['L4 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    print(L1_durations)
    
    plt.show()
    plt.ion()
    
    plt.gcf().clf()
    fig_h,ax_h = plt.subplots(5,1,sharex=True)

    ax_h[0].hist(L1_durations)
    ax_h[0].set_xlabel('Duration (hr)')
    ax_h[0].set_ylabel('Frequency')
    ax_h[0].set_title('L1 (n={}) - Mean: {}, Std:{}'.format(np.size(L1_durations),np.mean(L1_durations),np.std(L1_durations)))

    ax_h[1].hist(L2_durations)
    ax_h[1].set_title('L2 (n={}) - Mean: {}, Std:{}'.format(np.size(L2_durations),np.mean(L2_durations),np.std(L2_durations)))

    ax_h[2].hist(L3_durations)
    ax_h[2].set_title('L3 - Mean: {}, Std:{}'.format(np.mean(L3_durations),np.std(L3_durations)))

    ax_h[3].hist(L4_durations)
    ax_h[3].set_title('L4 - Mean: {}, Std:{}'.format(np.mean(L4_durations),np.std(L4_durations)))

    ax_h[4].hist(larval_durations)
    ax_h[4].set_xlabel('Time to maturity (hr)')
    ax_h[4].set_title('Maturity (n={}) - Mean: {}, Std:{}'.format(np.size(larval_durations),np.mean(larval_durations),np.std(larval_durations)))

def quick_plot_lifespan(ann_fps, expt_mds, bad_worm_kws=[],debug_mode=False):
    my_ann_files = [annotation_file.AnnotationFile(ann_fp) for ann_fp in ann_fps]  
    
    timestamped_data = {}
    [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    for [expt_md_fp, ann_file] in zip(expt_mds, my_ann_files):
        if debug_mode: print(expt_md_fp)
        ann_file_data = ann_file.data_as_timestamps_simple(expt_md_fp)
        for expt_key in timestamped_data.keys():
            timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
    
    viable_worm = (timestamped_data['Hatch']!=-1) \
        & (timestamped_data['Death']!=-1) \
        & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    print(timestamped_data['Worm'][viable_worm])
    
    lifespan = (timestamped_data['Death']-timestamped_data['Hatch'])[viable_worm]/(3600*24) # Days
    sorted_ls = sorted(lifespan)
    prop_alive = 1 - (np.arange(start=0,stop=np.size(lifespan)))/np.size(lifespan)
    
    plt.show()
    plt.ion()
    
    plt.gcf().clf()
    fig_h, ax_h = plt.subplots(2,1,sharex=True)
    ax_h[0].plot(np.append([0],sorted_ls),np.append([1],prop_alive))
    ax_h[0].set_xlabel('Time since expt. start (d)')
    ax_h[0].set_ylabel('Proportion alive')
    ax_h[0].set_title('Survival curve - n = {}'.format(np.size(lifespan)))
    
    ax_h[1].hist(lifespan)
    ax_h[1].set_xlabel('Time to death (d)')
    ax_h[1].set_ylabel('Frequency')
    ax_h[1].set_title('Mean+/-STD: {:.2f}+/-{:.2f}d\nMedian:{}d'.format(np.mean(lifespan),np.std(lifespan),np.median(lifespan)))
    
def plot_lifespan(ann_fps, expt_mds, annotation_prefix_list = [], bad_worm_kws=[],debug_mode=False):
    if annotation_prefix_list == []:
        my_ann_files = [annotation_file.AnnotationFile(ann_fp) for ann_fp in ann_fps]
    else:
        my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix=prefix) for (ann_fp,prefix) in zip(ann_fps,annotation_prefix_list)]
    
    timestamped_data = {}
    [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    for [expt_md_fp, ann_file] in zip(expt_mds, my_ann_files):
        if debug_mode: print(expt_md_fp)
        if type(expt_mds[0]) == str:    # Simple - one md file
            ann_file_data = ann_file.data_as_timestamps_simple(expt_md_fp)
        if type(expt_mds[0]) == dict:   # Need to link multiple md files
            ann_file_data = ann_file.data_as_timestamps(expt_md_fp)
            
        for expt_key in timestamped_data.keys():
            timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
    
    viable_worm = (timestamped_data['Hatch']!=-1) \
        & (timestamped_data['Death']!=-1) \
        & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    print(timestamped_data['Worm'][viable_worm])
    
    lifespan = (timestamped_data['Death']-timestamped_data['Hatch'])[viable_worm]/(3600*24) # Days
    sorted_ls = sorted(lifespan)
    prop_alive = 1 - (np.arange(start=0,stop=np.size(lifespan)))/np.size(lifespan)
    
    plt.show()
    plt.ion()
    
    plt.gcf().clf()
    fig_h, ax_h = plt.subplots(2,1,sharex=True)
    ax_h[0].plot(np.append([0],sorted_ls),np.append([1],prop_alive))
    ax_h[0].set_xlabel('Time since expt. start (d)')
    ax_h[0].set_ylabel('Proportion alive')
    ax_h[0].set_title('Survival curve - n = {}'.format(np.size(lifespan)))
    
    ax_h[1].hist(lifespan)
    ax_h[1].set_xlabel('Time to death (d)')
    ax_h[1].set_ylabel('Frequency')
    ax_h[1].set_title('Mean+/-STD: {:.2f}+/-{:.2f}d\nMedian:{}d'.format(np.mean(lifespan),np.std(lifespan),np.median(lifespan)))


def clean_plot(my_plot,make_labels=False,suppress_ticklabels=False):
    my_plot.spines['right'].set_visible(False)
    my_plot.spines['top'].set_visible(False)
    
    my_plot.tick_params(axis='both',which='both', top='off', bottom='off', left='off', right='off')
    
    #if not suppress_ticklabels:
    full_xticks = my_plot.get_xticks().copy()
    full_xticklabels = my_plot.get_xticklabels().copy()
    if full_xticks[0]*full_xticks[-1] < 0:
        my_plot.set_xticks([full_xticks[0], 0, full_xticks[-1]])
    else:
        my_plot.set_xticks([full_xticks[0], full_xticks[-1]])
    
    full_yticks = my_plot.get_yticks().copy()
    full_yticklabels = my_plot.get_yticklabels().copy()
    if full_yticks[0]*full_yticks[-1] < 0:
        my_plot.set_yticks([full_yticks[0], 0, full_yticks[-1]])
    else:
        my_plot.set_yticks([full_yticks[0], full_yticks[-1]])
    #else:
    if suppress_ticklabels:
        my_plot.set_xticks([])
        my_plot.set_yticks([])
    
    if not make_labels:
        my_plot.set_xlabel('')
        my_plot.set_ylabel('')
        my_plot.set_title('')

def force_same_plot_attributes(my_axes, *args):
    for attr in args:
        if attr=='xlim':
            min_x = min([ax.get_xlim()[0] for ax in my_axes])
            max_x = max([ax.get_xlim()[1] for ax in my_axes])
            [ax.set_xlim([min_x,max_x]) for ax in my_axes]
        elif attr=='ylim':
            min_y = min([ax.get_ylim()[0] for ax in my_axes])
            max_y = max([ax.get_ylim()[1] for ax in my_axes])
            [ax.set_ylim([min_y,max_y]) for ax in my_axes]

# Take a list and unfold it all the way
def flatten_list(my_list, to_level=-1, this_level=0):
    import collections
    flat_list = []
    
    for my_item in my_list:
        if (not isinstance(my_item, collections.Iterable)) or (type(my_item) is str) or ((to_level is not -1) and (this_level is to_level)):
            flat_list.append(my_item)
        else:
            flat_list.extend(flatten_list(my_item,to_level,this_level+1))
    return flat_list
