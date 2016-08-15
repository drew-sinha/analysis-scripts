import numpy as np
import matplotlib.pyplot as plt
import annotation_file
from scipy import stats
import os

if __name__=='__main__':
    do_varplot = False
    do_scatterplots = False
    do_indivplot = True
    
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
    bad_worm_kws = ['FERTILITY', 'Nh', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS']
    
    timestamped_data = {}
    
    my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix='D') for ann_fp in ann_fps]
    [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    for [expt_md_fps, ann_file] in zip(expt_mds, my_ann_files):
        ann_file_data = ann_file.data_as_timestamps(expt_md_fps)
        for expt_key in timestamped_data.keys():
            timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
            
    print(timestamped_data) # Check whether string or ints....
    # Define some filters for data
    viable_worm = (timestamped_data['Hatch']!=-1) \
        & (timestamped_data['Death']!=-1) \
        & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    hatched_on_corral = timestamped_data['Hatch'] != 0
    
    # Define durations used to look at variability in each stage (combining all individuals that were viable, i.e. hatched and developed); durations in hours
    L1_durations = (timestamped_data['L1 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    L2_durations = (timestamped_data['L2 ecdysis']-timestamped_data['L1 ecdysis'])[viable_worm]/3600
    L3_durations = (timestamped_data['L3 ecdysis']-timestamped_data['L2 ecdysis'])[viable_worm]/3600
    L4_durations = (timestamped_data['L4 ecdysis']-timestamped_data['L3 ecdysis'])[viable_worm]/3600
    larval_durations = (timestamped_data['L4 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    
    # Do analysis on only those guys that we could see hatch on the corrals; durations in COLUMNS
    combined_durations = np.array([(timestamped_data['L1 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600,
        (timestamped_data['L2 ecdysis']-timestamped_data['L1 ecdysis'])[viable_worm & hatched_on_corral]/3600,
        (timestamped_data['L3 ecdysis']-timestamped_data['L2 ecdysis'])[viable_worm & hatched_on_corral]/3600,
        (timestamped_data['L4 ecdysis']-timestamped_data['L3 ecdysis'])[viable_worm & hatched_on_corral]/3600]).T

    # Look at lifespan
    lifespan = (timestamped_data['Death']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    egg_maturity = (timestamped_data['First Egg Laid']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600

    plt.show()
    plt.ion()

    if do_varplot:
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

    
    
    #Correlation analysis on durations
    print('spearman')
    (corr_vals,corr_p_vals) = stats.spearmanr(combined_durations,axis=0)
    print(corr_vals)
    print(corr_p_vals)
    
    print('pearson')
    corr_vals = np.corrcoef(combined_durations,rowvar=0)
    print(corr_vals)
    
    print('Kendall-Tau')
    corr_vals = np.zeros((4,4))
    p_vals = np.zeros((4,4))
    for ii in range(4):
        for jj in range(4):
            corr_vals[ii,jj], p_vals[ii,jj] = stats.kendalltau(combined_durations[:,ii],combined_durations[:,jj])
    print(corr_vals)
    print(p_vals)
    
    if do_scatterplots:
        fig_h, ax_h = plt.subplots(4,4)
        for stg_1 in range(3):  # 1st 3 stages (4th doesn't give a new comparison)
            for stg_2 in range(stg_1+1,4):
                ax_h[stg_1,stg_2].scatter(combined_durations[:,stg_1],combined_durations[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} vs. L{}'.format(stg_1+1,stg_2+1))
    
    if do_indivplot:
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(range(len(worm_stgdata)),worm_stgdata,marker='o') for worm_stgdata in combined_durations]
        #ax_h.plot([0,3],[0,0],'k',linestyle='--',linewidth=3.5)
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        ax_h.set_ylabel('Time in Larval Stage')
        
        # Plot larval stage durations by individual (Z-score)
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(range(len(worm_stgdata)),worm_stgdata,marker='o') for worm_stgdata in (combined_durations-np.mean(combined_durations,axis=0))/np.std(combined_durations,axis=0)]
        ax_h.plot([0,3],[0,0],'k',linestyle='--',linewidth=3.5)
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        ax_h.set_ylabel('Time in Larval Stage (z-score w.r.t. population)')
        #fig_h.savefig('./20160227_DevVarIndividualTiming.png')
        
    # First look at total larval growth against lifespan
    fig_h = plt.figure()
    ax_h = fig_h.gca()
    ax_h.scatter(larval_durations, lifespan)
    ax_h.set_xlabel('Total time in larval development (hr)')
    ax_h.set_ylabel('Total lifespan (hr)')
    
    # First look at total larval growth against lifespan
    fig_h = plt.figure()
    ax_h = fig_h.gca()
    ax_h.scatter(egg_maturity, lifespan)
    ax_h.set_xlabel('Time to first egg laid (hr)')
    ax_h.set_ylabel('Total lifespan (hr)')
    
    # Look at growth for each stage against lifespan
    fig_h, ax_h =plt.subplots(1,4)
    for stg in range(4):
        ax_h[stg].scatter(combined_durations[:,stg], lifespan)
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].set_xlabel('Time in stage (hr)')
        ax_h[stg].set_ylabel('Lifespan (hr)')
        
    # Look at total time spent away from mean lifespan
    fig_h = plt.figure()
    ax_h = fig_h.gca()
    combined_durations_mc = (combined_durations-np.mean(combined_durations,axis=0))
    ax_h.scatter(np.sum(np.abs(combined_durations_mc),axis=1), lifespan)
    ax_h.set_xlabel('total deviation from mean development throughout life (hr)')
    ax_h.set_ylabel('Lifespan (hr)')
