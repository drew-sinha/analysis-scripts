import freeimage
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import csv

import annotation_file
import process_data

import pickle

def find_char(string,ch):
    return [idx for (idx, letter) in enumerate(string) if letter == ch]

if __name__ == "__main__":
    # Pull out timing data from annotation files
    ann_fps = ['/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12A/2016.02.29 spe-9 12A Worm Annotations - WormNotes.tsv',
        '/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12B/2016.02.29 spe-9 12B Worm Annotations - WormNotes.tsv',
        '/mnt/scopearray/ZhangWillie/2016.03.04 spe-9 Run 13B/2016.03.04 spe-9 13B Worm Annotations - WormNotes.tsv',
        '/mnt/scopearray/ZhangWillie/2016.03.04 spe-9 Run 13C/2016.03.04 spe-9 13C Worm Annotations - WormNotes.tsv']

    expt_mds = [
        {'D':'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/experiment_metadata.json',
        'W':'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12A/experiment_metadata.json'},
        {'D':'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB/experiment_metadata.json',
        'W':'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12B/experiment_metadata.json'},
        {'D':'/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarB/experiment_metadata.json',
        'W':'/mnt/scopearray/ZhangWillie/2016.03.04 spe-9 Run 13B/experiment_metadata.json'},
        {'D':'/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarC/experiment_metadata.json',
        'W':'/mnt/scopearray/ZhangWillie/2016.03.04 spe-9 Run 13C/experiment_metadata.json'}
    ]
    expt_paths = ['/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/',
        '/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB/',
        '/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarB/',
        '/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarC/']
    bad_worm_kws = ['FERTILITY', 'Nh', 'NO HATCH', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS']
    
    # Grab lifespans
    # Pull annotations and convert them to timestamps
    #Good State
    #my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix='D') for ann_fp in ann_fps]
    #timestamped_data = [[] for idx in range(len(expt_paths))]
    #for expt_num, [expt_md_fps, ann_file] in enumerate(zip(expt_mds, my_ann_files)):
    #    timestamped_data[expt_num] = ann_file.data_as_timestamps(expt_md_fps)    
    
    my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix='D') for ann_fp in ann_fps]
    timestamped_data = [[] for idx in range(len(expt_paths))]
    for expt_num, [expt_md_fps, ann_file,ann_fp] in enumerate(zip(expt_mds, my_ann_files,ann_fps)):
        timestamped_data[expt_num] = ann_file.data_as_timestamps(expt_md_fps,ann_fp[find_char(ann_fp,os.path.sep)[-1]+1:find_char(ann_fp,' ')[-4]]+' ')

    compiled_expts = [process_data.Experiment(expt_path) for expt_path in expt_paths]
    devtiming_allexpts = [[] for idx in range(len(expt_paths))]
    totaldev_allexpts = [[] for idx in range(len(expt_paths))]
    goodworms_allexpts = []

    for expt_num, [expt, expt_events] in enumerate(zip(compiled_expts, timestamped_data)):
        acq_worms = expt.get_acquired_wells()
        first_worm_num = int(expt_events['Worm'][0][-2:])   # Add in adjustment for the first index of the worms not being at 0
        viable_worm = (expt_events['Hatch']!=-1) \
            & (expt_events['Death']!=-1) \
            & ((expt_events['L1 ecdysis']-expt_events['Hatch'])/3600 <= 15) \
            & np.array([not any([kw in note for kw in bad_worm_kws]) for note in expt_events['Notes']])
        hatched_on_corral = timestamped_data[expt_num]['Hatch'] != 0
        worm_was_acquired = [str(worm_num+first_worm_num).zfill(len(expt.get_acquired_wells()[0])) in expt.get_acquired_wells() for worm_num in np.arange(len(expt_events['Hatch']))]

        print(np.count_nonzero(viable_worm & hatched_on_corral & worm_was_acquired))
        
        # Calculate dev timing for worms
        #lifespan_allexpts[expt_num] = (expt_events['Death']-expt_events['Hatch'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        
        # Define durations used to look at variability in each stage (combining all individuals that were viable, i.e. hatched and developed and acquired); durations in hours
        L1_durations = (expt_events['L1 ecdysis']-expt_events['Hatch'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        L2_durations = (expt_events['L2 ecdysis']-expt_events['L1 ecdysis'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        L3_durations = (expt_events['L3 ecdysis']-expt_events['L2 ecdysis'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        L4_durations = (expt_events['L4 ecdysis']-expt_events['L3 ecdysis'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
    
        # Do analysis on only those guys that we could see hatch on the corrals; durations in COLUMNS
        devtiming_allexpts[expt_num] = np.array([L1_durations,L2_durations,L3_durations,L4_durations]).T
        totaldev_allexpts[expt_num] = (expt_events['L4 ecdysis']-expt_events['Hatch'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        goodworms_allexpts.append([worm_name.replace('/','') for worm_name in expt_events['Worm_FullName'][viable_worm & hatched_on_corral & worm_was_acquired]])
    
    #lifespans = np.concatenate(lifespan_allexpts,axis=0)
    devtiming = np.concatenate(devtiming_allexpts,axis=0)
    totaldev = np.concatenate(totaldev_allexpts,axis=0)
    goodworms = [worm for good_group in goodworms_allexpts for worm in good_group]
    print(devtiming_allexpts)
    print(devtiming)
    
    bob = input('stop here')
    
    ##################################
    # Pull out egg data from Willie data
    
    # Compile all worm parameters into single data obj idx'd by worm (not experiment)
    worm_health = []
    longterm_dir = '/media/Data/Work/ZPLab/WormImages/spe9_Longterm_Data/semiprocessed_cheese_product/'
    for worm_fullname in goodworms:
        # Load the appropriate file of parameters
        with open(longterm_dir+os.path.sep+worm_fullname+'.tsv','r') as worm_file:
            worm_filereader = csv.reader(worm_file,delimiter='\t')
            num_vals = sum(1 for row in worm_filereader)-1  # First row is header
        
        with open(longterm_dir+os.path.sep+worm_fullname+'.tsv','r') as worm_file:
            worm_filereader = csv.reader(worm_file,delimiter='\t')
            tags = next(worm_filereader)   
            tags[0] = 'time'    # No first column header
            
            
            worm_data = {tag:np.zeros((num_vals)) for tag in tags}
            for row_num,row in enumerate(worm_filereader):
                for item_num, item in enumerate(row):
                    worm_data[tags[item_num]][row_num] = float(item) if item is not '' else -1.0
            
            # Finally add it to accumulating list
            worm_health.append(worm_data)
    
    totaleggs = np.array([worm['cumulative_eggs'].max() for worm in worm_health])
    eggs_firstday = np.zeros_like(totaleggs)
    for worm_idx, worm in enumerate(worm_health):
        first_egg_laid = np.where((worm['egg_age']<0) & (worm['egg_age'] != -1))[0][-1]
        end_of_day1 = np.where((worm['egg_age']<24) & (worm['egg_age'] != -1))[0][-1]
        eggs_firstday[worm_idx] = worm['cumulative_eggs'][end_of_day1]-worm['cumulative_eggs'][first_egg_laid]
        
    ##################################
    fig_h, ax_h = plt.subplots(1,5)
    for stg in range(5):
        ax_h[stg].set_ylabel('Total eggs laid')
        if stg < 4: 
            ax_h[stg].set_title('Eggs laid vs. L{} timing\n Spearman r:{:.2f} p={:.2f}'.format(stg+1,*scipy.stats.spearmanr(devtiming[:,stg],totaleggs)))
            ax_h[stg].set_xlabel('Total time spent in stage')
            ax_h[stg].scatter(devtiming[:,stg],totaleggs)
        else:
            ax_h[stg].set_title('Eggs laid vs. total time spent in development\n Spearman r:{:.2f} p={:.2f}'.format(*scipy.stats.spearmanr(totaldev,totaleggs)))
            ax_h[stg].set_xlabel('Total developmental time (hr)')
            ax_h[stg].scatter(totaldev,totaleggs)
    #fig_h.savefig('./20160506_DevTimingvsEgglaying.png')
    
    fig_h, ax_h = plt.subplots(1,5)
    for stg in range(5):
        ax_h[stg].set_ylabel('Eggs on first day')
        if stg < 4: 
            ax_h[stg].set_title('Eggs 1st day vs. L{} timing\n Spearman r:{:.2f} p={:.2f}'.format(stg+1,*scipy.stats.spearmanr(devtiming[:,stg],eggs_firstday)))
            ax_h[stg].set_xlabel('Total time spent in stage')
            ax_h[stg].scatter(devtiming[:,stg],eggs_firstday)
        else:
            ax_h[stg].set_title('Eggs 1st day vs. total time spent in development\n Spearman r:{:.2f} p={:.2f}'.format(*scipy.stats.spearmanr(totaldev,eggs_firstday)))
            ax_h[stg].set_xlabel('Total developmental time (hr)')
            ax_h[stg].scatter(totaldev,eggs_firstday)

    ##################################
    # Pull out size data from masks
    top_dir = '/media/Data/Work/ZPLab/WormImages/spe9Acquisition_DevImages_hmasks/'
    
    all_worm_sizes = np.empty((0,5))
    all_worm_sizechgs = np.empty((0,4))
    all_worm_devtime = np.empty((0,4))    # time in each larval stage
    all_worm_rates = np.empty((0,4))
    
    for expt_dir in [subdir for subdir in sorted(os.listdir(top_dir)) if os.path.isdir(top_dir+os.path.sep+subdir)]:
        for worm_dir in sorted(os.listdir(top_dir+os.path.sep+expt_dir)):
            
            #worm_sizes = np.zeros((5,1))
            #for stg in range(1,
            
            #Note: Skip first file
            #for img_num, img_fn in enumerate([img_file in sorted(os.listdir(top_dir+os.path.sep+expt_dir+os.path.sep+worm_dir+os.path.sep)) if 'hmask' in img_file][1:]):
            #    worm_image = freeimage.read(top_dir+os.path.sep+expt_dir+os.path.sep+worm_dir+os.path.sep+img_fn)
            worm_mask_files = [img_file for img_file in sorted(os.listdir(top_dir+os.path.sep+expt_dir+os.path.sep+worm_dir+os.path.sep)) if 'hmask' in img_file][1:]
            worm_images = [freeimage.read(top_dir+os.path.sep+expt_dir+os.path.sep+worm_dir+os.path.sep+img_fn) for \
                img_fn in worm_mask_files]
            worm_image_times = [process_data.extract_datetime_fromstr(mask_file[0:15]) for mask_file in worm_mask_files]    # len=5
            time_diffs = np.array([(expt_time-worm_image_times[0]).total_seconds()/3600 for expt_time in worm_image_times]) #len=5 with first item =0
            all_worm_sizes = np.append(all_worm_sizes, [[np.count_nonzero(w_img) for w_img in worm_images]],axis=0)
            all_worm_sizechgs = np.append(all_worm_sizechgs, [np.diff(all_worm_sizes[-1,:])], axis=0)
            all_worm_devtime = np.append(all_worm_devtime, [np.diff(time_diffs)],axis=0)
            all_worm_rates = np.append(all_worm_rates, [np.diff(all_worm_sizes[-1,:])/all_worm_devtime[-1,:]],axis=0)
            
            # Censor worms with too long taken to develop
            all_worm_sizes = all_worm_sizes[all_worm_devtime[:,0]<=15,:]
            all_worm_sizechgs = all_worm_sizechgs[all_worm_devtime[:,0]<=15,:]
            all_worm_rates = all_worm_rates[all_worm_devtime[:,0]<=15,:]
            all_worm_devtime = all_worm_devtime[all_worm_devtime[:,0]<=15,:]
    #print(all_worm_sizechgs)
    
    sizes_zscored = (all_worm_sizes-np.mean(all_worm_sizes,axis=0))/np.std(all_worm_sizes,axis=0)
    rates_zscored = (all_worm_rates-np.mean(all_worm_rates,axis=0))/np.std(all_worm_rates,axis=0)
    sizechgs_zscored = (all_worm_sizechgs-np.mean(all_worm_sizechgs,axis=0))/np.std(all_worm_sizechgs,axis=0)
    devtime_zscored = (all_worm_devtime-np.mean(all_worm_devtime,axis=0))/np.std(all_worm_devtime,axis=0)
    
    ##########################
    fig_h, ax_h = plt.subplots(1,4)
    for stg in range(4):
        ax_h[stg].set_ylabel('Total eggs laid')
        ax_h[stg].set_title('Eggs laid vs. L{} size\n Spearman r:{:.2f} p={:.2f}'.format(stg+1,*scipy.stats.spearmanr(all_worm_sizes[:,stg+1],totaleggs)))
        ax_h[stg].set_xlabel('Size')
        ax_h[stg].scatter(all_worm_sizes[:,stg+1],totaleggs)
        
    fig_h, ax_h = plt.subplots(1,4)
    for stg in range(4):
        ax_h[stg].set_ylabel('Eggs 1st day')
        ax_h[stg].set_title('Eggs 1st day vs. L{} size\n Spearman r:{:.2f} p={:.2f}'.format(stg+1,*scipy.stats.spearmanr(all_worm_sizes[:,stg+1],eggs_firstday)))
        ax_h[stg].set_xlabel('Size')
        ax_h[stg].scatter(all_worm_sizes[:,stg+1],eggs_firstday)
        
    fig_h, ax_h = plt.subplots(1,4)
    for stg in range(4):
        ax_h[stg].set_ylabel('Total eggs laid')
        ax_h[stg].set_title('Eggs laid vs. L{} growth\n Spearman r:{:.2f} p={:.2f}'.format(stg+1,*scipy.stats.spearmanr(all_worm_rates[:,stg],totaleggs)))
        ax_h[stg].set_xlabel('Average growth rate in stage')
        ax_h[stg].scatter(all_worm_rates[:,stg],totaleggs)
    
    fig_h, ax_h = plt.subplots(1,4)
    for stg in range(4):
        ax_h[stg].set_ylabel('Eggs 1st day')
        ax_h[stg].set_title('Eggs 1st day vs. L{} growth\n Spearman r:{:.2f} p={:.2f}'.format(stg+1,*scipy.stats.spearmanr(all_worm_rates[:,stg],eggs_firstday)))
        ax_h[stg].set_xlabel('Average growth rate in stage')
        ax_h[stg].scatter(all_worm_rates[:,stg],eggs_firstday)


    ####################
    '''
    # For pickling
    with open('/media/Data/Work/ZPLab/Analysis/DevelopmentalVariability/Size/HumanMasks_WithCensoring/WormGrowthData_ForZach.pickle','wb') as my_file:
        pickle.dump([all_worm_sizes, all_worm_sizechgs, all_worm_devtime, all_worm_rates, goodworms], my_file)
    
    # For unpickling
    with open('/media/Data/Work/ZPLab/Analysis/DevelopmentalVariability/Size/HumanMasks_WithCensoring/WormGrowthData_ForZach.pickle','rb') as my_file:
        [all_worm_sizes, all_worm_sizechgs, all_worm_devtime, all_worm_rates, goodworms] = pickle.load(my_file)
    '''




