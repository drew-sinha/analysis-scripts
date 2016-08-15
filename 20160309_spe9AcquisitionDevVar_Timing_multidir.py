import numpy as np
import matplotlib.pyplot as plt
import annotation_file
from scipy import stats
import os

if __name__=='__main__':
    do_varplot = True
    do_scatterplots = True
    do_indivplot = True
    
    expt_dirs = ['/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarC/',
        '/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarB/', 
        '/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/', 
        '/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB/']
    
    timestamped_data = {}
    
    my_ann_files = [annotation_file.AnnotationFile(expt_dir+os.path.sep+'worm_annotations.tsv') for expt_dir in expt_dirs]
    [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    for [expt_dir, ann_file] in zip(expt_dirs, my_ann_files):
        print(expt_dir)
        ann_file_data = ann_file.data_as_timestamps_simple(expt_dir+os.path.sep+'experiment_metadata.json')
        for expt_key in timestamped_data.keys():
            timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
            
    print(timestamped_data) # Check whether string or ints....
    # Define some filters for data
    viable_worm = timestamped_data['Hatch']!=-1
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
        # Plot larval stage durations by individual
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(range(len(worm_stgdata)),worm_stgdata,marker='o') for worm_stgdata in (combined_durations-np.mean(combined_durations,axis=0))/np.std(combined_durations,axis=0)]
        ax_h.plot([0,3],[0,0],'k',linestyle='--',linewidth=3.5)
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        ax_h.set_ylabel('Time in Larval Stage (z-score w.r.t. population)')
        #fig_h.savefig('./20160227_DevVarIndividualTiming.png')
