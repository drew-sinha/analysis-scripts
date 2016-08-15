import numpy as np
import matplotlib.pyplot as plt
import itertools

import process_data
import annotation_file

if __name__=='__main__':
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
        '/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarC/',]
    #ann_fps = ['/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12A/2016.02.29 spe-9 12A Worm Annotations - WormNotes.tsv',
        #'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12B/2016.02.29 spe-9 12B Worm Annotations - WormNotes.tsv']

    #expt_mds = [
        #{'D':'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/experiment_metadata.json',
        #'W':'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12A/experiment_metadata.json'},
        #{'D':'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB/experiment_metadata.json',
        #'W':'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12B/experiment_metadata.json'}
    #]
    #expt_paths = ['/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/']
    #expt_paths = ['/media/Data/Work/ZPLab/WormImages/20160229_spe9Acquisition_DevVarA_temp/',
        #'/media/Data/Work/ZPLab/WormImages/20160229_spe9Acquisition_DevVarB_temp/'
    #]
    #expt_paths = ['/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/',
        #'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB/'
    #]
    bad_worm_kws = ['FERTILITY', 'Nh', 'NO HATCH', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS']
    
    plt.show()
    plt.ion()
    
    
    # Pull annotations and convert them to timestamps
    my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix='D') for ann_fp in ann_fps]
    timestamped_data = [[] for idx in range(len(expt_paths))]
    for expt_num, [expt_md_fps, ann_file] in enumerate(zip(expt_mds, my_ann_files)):
        timestamped_data[expt_num] = ann_file.data_as_timestamps(expt_md_fps)

    # For each experiment, grab sizes and rates of growth for worms, as well as expt_times
    compiled_expts = [process_data.Experiment(expt_path) for expt_path in expt_paths]
    compiled_times = [expt.get_expt_time_offsets() for expt in compiled_expts]
    compiled_sizes = [expt.measure_size_forexpt_batch() for expt in compiled_expts]
    compiled_rates = [expt.measure_growthrate_forexpt_batch() for expt in compiled_expts]
    
    '''
    # Plot worms in aggregate by expt
    fig_h, ax_h = plt.subplots(2,1)
    fig_h.suptitle('Size and growth by expt')
    for expt_sizes, expt_grates, expt_times in zip(compiled_sizes, compiled_rates, compiled_times):
        #ax_h[0].plot(expt_times, expt_sizes.T)
        [ax_h[0].plot(expt_times[0:len(worm_size)], worm_size) for worm_size in expt_sizes]
        ax_h[0].set_xlabel('Time from expt_start (hr)')
        ax_h[0].set_ylabel('Size (px)')
        
        [ax_h[1].plot(expt_times[0:len(worm_grate)], worm_grate) for worm_grate in expt_grates]
        ax_h[1].set_xlabel('Time from expt start (hr)')
        ax_h[1].set_ylabel('Rate of growth (px/hr)')
    '''
    
    #Plot worm growth in aggregate by expt adjusted to birth
    fig_h, ax_h = plt.subplots(2,1, sharex=True)
    fig_h.suptitle('growth adjusted to birth')
    for (expt_sizes, expt_rates, expt_times, expt_events) in zip(compiled_sizes, compiled_rates, compiled_times, timestamped_data):
        #ax_h[0].plot(
            #np.tile(np.array([expt_times]).T,(1,len(expt_sizes)))-expt_events['Hatch'][0]/3600,expt_sizes.T)
        [ax_h[0].plot((expt_times[0:len(worm_size)]-expt_events['Hatch'][worm_num])/3600,worm_size) for worm_num,worm_size in enumerate(expt_sizes)]
        ax_h[0].set_xlabel('Time post hatch (hr)')
        ax_h[0].set_ylabel('Size (px)')
        
        #ax_h[1].plot(
            #(np.tile(np.array([expt_times]).T,(1,len(expt_sizes)))-expt_events['Hatch'][0]/3600)[0:-1],
            #expt_rates.T)
        [ax_h[1].plot((expt_times[0:len(worm_grate)]-expt_events['Hatch'][worm_num])/3600,worm_grate) for worm_num, worm_grate in enumerate(expt_rates)]
        ax_h[1].set_xlabel('Time post hatch (hr)')
        ax_h[1].set_ylabel('Rate of growth (px/hr)')
    
    # Compile data on timepoints for each individual
    dev_labels = np.array([expt.label_time_perstage(ann_file,expt_md_fps) for [expt, expt_md_fps, ann_file] in zip(compiled_expts, expt_mds, my_ann_files)])
    size_stats = [lambda data:data[-1]]
    rate_stats = [np.mean, lambda data: np.percentile(data, 90)]

    size_stat_allexpts = [[] for idx in range(len(expt_paths))]
    rate_stat_allexpts = [[] for idx in range(len(expt_paths))]
    lifespan_allexpts = [[] for idx in range(len(expt_paths))]
    
    for expt_num, [expt, expt_events, expt_sizes, expt_rates, expt_dev_label] in enumerate(zip(compiled_expts, timestamped_data, compiled_sizes, compiled_rates, dev_labels)):
        acq_worms = expt.get_acquired_wells()
        first_worm_num = int(expt_events['Worm'][0][-2:])   # Add in adjustment for the first index of the worms not being at 0
        viable_worm = (expt_events['Hatch']!=-1) \
            & (expt_events['Death']!=-1) \
            & np.array([not any([kw in note for kw in bad_worm_kws]) for note in expt_events['Notes']])
        hatched_on_corral = timestamped_data[expt_num]['Hatch'] != 0
        worm_was_acquired = [str(worm_num+first_worm_num).zfill(len(expt.get_acquired_wells()[0])) in expt.get_acquired_wells() for worm_num in np.arange(len(expt_events['Hatch']))]

        print(np.count_nonzero(viable_worm & hatched_on_corral & worm_was_acquired))
        
        # Calculate lifespan for worms
        lifespan_allexpts[expt_num] = (expt_events['Death']-expt_events['Hatch'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        
        size_stat_allexpts[expt_num] = np.zeros((len(size_stats), np.count_nonzero(viable_worm & hatched_on_corral & worm_was_acquired), 4))
        rate_stat_allexpts[expt_num] = np.zeros((len(rate_stats), np.count_nonzero(viable_worm & hatched_on_corral & worm_was_acquired), 4))
        for worm_num, worm in enumerate(np.where(viable_worm & hatched_on_corral & worm_was_acquired)[0]):
            worm_idx = acq_worms.index(str(worm+first_worm_num).zfill(len(acq_worms[0])))  # idx in sizes and labels matrix
            #print(worm_idx)
            for stg in np.arange(start=1, stop=5): # Only larval
                #print(stg)
                for stat_num, stat in enumerate(size_stats):
                    size_stat_allexpts[expt_num][stat_num,worm_num,stg-1] = stat(expt_sizes[worm_idx][expt_dev_label[worm_idx,0:len(expt_sizes[worm_idx])]==stg])
                
                for stat_num, stat in enumerate(rate_stats):
                    rate_stat_allexpts[expt_num][stat_num,worm_num,stg-1] = stat(expt_rates[worm_idx][expt_dev_label[worm_idx,0:len(expt_rates[worm_idx])]==stg])
    
    # Collapse list across experiments to get a matrix idx'd by (stat, worm, stage)
    size_stat_comb = np.concatenate(size_stat_allexpts,axis=1)
    rate_stat_comb = np.concatenate(rate_stat_allexpts,axis=1)
    lifespan = np.concatenate(lifespan_allexpts,axis=0)
    print(lifespan.shape)
    
    # Plot ending size
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('Size at end of stage')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(size_stat_comb[0,:,stg],lifespan)
        ax_h[stg].set_xlabel('Size (px)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
    
    # Plot average growth rate over stage
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('Mean growth rate in stage')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(rate_stat_comb[0,:,stg],lifespan)
        ax_h[stg].set_xlabel('Growth rate (px/hr)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
    
    # Plot 90 percentile growth rate over time
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('90th percentile growth rate in stage')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(rate_stat_comb[1,:,stg],lifespan)
        ax_h[stg].set_xlabel('Size (px)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
        
    
'''
    Things I want to do do:
        1. Look at variability in size throughout the entire population
            Individuals (per worm and/or in aggregate); per expt DONE
            Align to birth DONE
            Median +/- Iqr
        2. Using annotations, curate final_size, average rate of growth, and percentile rate of growth (90th percentile) for each animal in each larval stage DONE
        3. Scatter final size and rate of growth against lifespan (pulled from annotations) DONE
'''
    
