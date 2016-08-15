import freeimage
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import csv

import annotation_file
import process_data

import sklearn.linear_model

'''
    Notes:
        - Censoring happens two places(!) Once when loading data from images originally, and the second when loading the lifespan data
'''


####
def partialcoef(var1, var2, control, method = 'spearman'):
    # Handles univariate data; calculates partial r for two variables holding one variable constant
    
    if method == 'spearman':
        return partialcoef(scipy.stats.rankdata(var1,method='ordinal'),
            scipy.stats.rankdata(var2,method='ordinal'),
            scipy.stats.rankdata(control,method='ordinal'),
            method='pearson')
    
    partial_r = (scipy.stats.pearsonr(var1,var2)[0]-scipy.stats.pearsonr(var1,control)[0]*scipy.stats.pearsonr(var2,control)[0])/ \
        np.sqrt((1-scipy.stats.pearsonr(var1,control)[0])**2 * (1-scipy.stats.pearsonr(var2,control)[0])**2)
    
    # Significance (approximately distributed as a t-distribution with n-2-k df
    t_obs = partial_r*np.sqrt((len(var1)-2-1)/(1-partial_r**2))
    p_val = scipy.stats.t.sf(np.abs(t_obs), len(var1)-2-1)   #  Two tailed (uses abs)
    return (partial_r, p_val)
    
def find_char(string,ch):
    return [idx for (idx, letter) in enumerate(string) if letter == ch]


if __name__=='__main__':
    do_sizes = True
    do_sizechgs = True
    do_times = True
    do_rates = True
    
    do_totalcorrelation=True
    do_withincorrelation=True
    do_betweencorrelation=True
    
    pxsize_conversion=1.304/1000    # mm; i.e. 1 px size is ~1.304 um
    
    plt.ion()
    plt.show()
    
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
            
            all_worm_sizes = np.append(all_worm_sizes, [[np.count_nonzero(w_img)*(pxsize_conversion**2) for w_img in worm_images]],axis=0)
            all_worm_sizechgs = np.append(all_worm_sizechgs, [np.diff(all_worm_sizes[-1,:])], axis=0)
            all_worm_devtime = np.append(all_worm_devtime, [np.diff(time_diffs)],axis=0)
            all_worm_rates = np.append(all_worm_rates, [np.diff(all_worm_sizes[-1,:])/all_worm_devtime[-1,:]],axis=0)
            

            
            # Censor worms with too long taken to develop
            #all_worm_sizes = all_worm_sizes[all_worm_devtime[:,0]<=15,:]
            #all_worm_sizechgs = all_worm_sizechgs[all_worm_devtime[:,0]<=15,:]
            #all_worm_rates = all_worm_rates[all_worm_devtime[:,0]<=15,:]
            #all_worm_devtime = all_worm_devtime[all_worm_devtime[:,0]<=15,:]
    #print(all_worm_sizechgs)
    
    sizes_zscored = (all_worm_sizes-np.mean(all_worm_sizes,axis=0))/np.std(all_worm_sizes,axis=0)
    rates_zscored = (all_worm_rates-np.mean(all_worm_rates,axis=0))/np.std(all_worm_rates,axis=0)
    sizechgs_zscored = (all_worm_sizechgs-np.mean(all_worm_sizechgs,axis=0))/np.std(all_worm_sizechgs,axis=0)
    devtime_zscored = (all_worm_devtime-np.mean(all_worm_devtime,axis=0))/np.std(all_worm_devtime,axis=0)
    
    sizechgs_relative = all_worm_sizechgs/all_worm_sizes[:,1:]
    rates_relative = all_worm_rates/all_worm_sizes[:,1:]
    sizechgs_relative_zscored = (sizechgs_relative-np.mean(sizechgs_relative,axis=0))/np.std(sizechgs_relative,axis=0)
    rates_relative_zscored = (rates_relative-np.mean(rates_relative,axis=0))/np.std(rates_relative,axis=0)
    
    if do_withincorrelation:
        ###############
        # Plot overall distribution of worm sizes per stage
        plt.close('all')
        plt.gcf().clf()
        fig_h,ax_h = plt.subplots(4,1)
        fig_h.suptitle('Size distributions')
        for stg in range(4):
            ax_h[stg].hist(all_worm_sizes[:,stg+1])  #Hatch is first column
            ax_h[stg].set_title('L{} - Mean: {:.4f}, Std: {:.4f}'.format(stg+1,np.mean(all_worm_sizes[:,stg+1]),np.std(all_worm_sizes[:,stg+1])))
            #ax_h[stg].set_title('L{}'.format(stg+1))
            if stg == 3:
                ax_h[stg].set_xlabel('Size (mm^2)')
                ax_h[stg].set_ylabel('Frequency')
        
        # Plot sizes over time per individual
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_sizes) for worm_sizes in all_worm_sizes]
        ax_h.set_title('Size over time per individual')
        ax_h.set_ylabel('Mask size (mm^2)')
        ax_h.set_xlim([0,4])
        ax_h.set_xticks(np.linspace(0,4,5))
        ax_h.xaxis.set_ticklabels(['Hatch','End_L1','E_L2','E_L3','E_L4'])
        
        # Plot sizes over time as adjusted z_scores
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_sizes) for worm_sizes in (all_worm_sizes-np.mean(all_worm_sizes, axis=0))/np.std(all_worm_sizes,axis=0)]
        ax_h.set_title('Size over time (z-scored)')
        ax_h.set_ylabel('Z-score')
        ax_h.set_xlim([0,4])
        ax_h.set_xticks(np.linspace(0,4,5))
        ax_h.xaxis.set_ticklabels(['Hatch','End_L1','E_L2','E_L3','E_L4'])
        
        # Spearman correlations on sizes
        print('Spearman (sizes)')
        (corr_vals, corr_p_vals) = scipy.stats.spearmanr(all_worm_sizes[:,1:])
        print(corr_vals)
        print(corr_p_vals)
        
        bob = input('To move to size changes')
        
        bob = input('To move to stage durations')
        
        #######################
        # Plot distribution of stage durations over all individuals
        plt.close('all')
        plt.gcf().clf()
        fig_h,ax_h = plt.subplots(4,1)
        fig_h.suptitle('Larval stage duration')
        for stg in range(4):
            ax_h[stg].hist(all_worm_devtime[:,stg])
            ax_h[stg].set_title('L{} - Mean: {:.2f}, Std: {:.2f}'.format(stg+1,np.mean(all_worm_devtime[:,stg]),np.std(all_worm_devtime[:,stg])))
            #ax_h[stg].set_title('L{}'.format(stg+1))
            if stg == 0:
                ax_h[stg].set_xlabel('Stage duration (hr)')
                ax_h[stg].set_ylabel('Frequency')
        
        # Plot stage durations per individual
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_devtime) for worm_devtime in all_worm_devtime]
        ax_h.set_title('Dev. stage duration per individual')
        ax_h.set_ylabel('Stage duration (hr)')
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        
        # Plot sizes over time as adjusted z_scores
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_devtime) for worm_devtime in (all_worm_devtime-np.mean(all_worm_devtime, axis=0))/np.std(all_worm_devtime,axis=0)]
        ax_h.set_title('Dev. stage duration (z-scored)')
        ax_h.set_ylabel('Duration (Z-score)')
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        
        # Scatter stage durations against each other
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Stage durations')
        for stg_1 in range(3):
            for stg_2 in range(stg_1+1,4):
                ax_h[stg_1,stg_2].scatter(all_worm_devtime[:,stg_1],all_worm_devtime[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} vs. L{}'.format(stg_1+1,stg_2+1))
        
        # Spearman correlations on size changes
        print('Correlation on stage durations (Spearman)')
        (corr_vals, corr_p_vals) = scipy.stats.spearmanr(all_worm_devtime)
        print(corr_vals)
        print(corr_p_vals)
        
        bob = input('To move to growth rates')
        
        ########################
        # Plot overall distribution of worm rates per stage
        plt.close('all')
        plt.gcf().clf()
        fig_h,ax_h = plt.subplots(4,1)
        for stg in range(4):
            ax_h[stg].hist(all_worm_rates[:,stg])  #Hatch is first column
            #ax_h[stg].set_title('L{} - Mean: {}, Std: {}'.format(stg+1,np.mean(all_worm_sizes[:,stg]),np.std(all_worm_sizes[:,stg])))
            ax_h[stg].set_title('L{}'.format(stg+1))
            if stg == 0:
                ax_h[stg].set_xlabel('Growth rate (px/hr)')
                ax_h[stg].set_ylabel('Frequency')
        
        # Plot rates over time per individual
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_rates) for worm_rates in all_worm_rates]
        ax_h.set_title('Growth rate per individual')
        ax_h.set_ylabel('Growth rate (px/hr)')
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        
        # Plot mean centered rates over time per individual
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_rates) for worm_rates in (all_worm_rates-np.mean(all_worm_rates, axis=0))]
        ax_h.set_title('Growth rate per stage (mean-centered)')
        ax_h.set_ylabel('Size relative to mean (px/hr)')
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4']) 
        
        # Plot rates over time as adjusted z_scores
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        [ax_h.plot(worm_rates) for worm_rates in (all_worm_rates-np.mean(all_worm_rates, axis=0))/np.std(all_worm_rates,axis=0)]
        ax_h.set_title('Growth rate over time (z-scored)')
        ax_h.set_ylabel('Growth rate (px/hr)')
        ax_h.set_xlim([0,3])
        ax_h.set_xticks(np.linspace(0,3,4))
        ax_h.xaxis.set_ticklabels(['L1','L2','L3','L4'])
        
        # Scatter rates against each other
        fig_h, ax_h = plt.subplots(4,4)
        for stg_1 in range(3):
            for stg_2 in range(stg_1+1,4):
                ax_h[stg_1,stg_2].scatter(all_worm_rates[:,stg_1],all_worm_rates[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} vs. L{}'.format(stg_1+1,stg_2+1))
        
        # Spearman correlations on rates
        print('Spearman (rates)')
        (corr_vals, corr_p_vals) = scipy.stats.spearmanr(all_worm_rates)
        print(corr_vals)
        print(corr_p_vals)

        bob = input('To correlations b/t variables')

    ###############
    # Correlations between separate variables (analyze temporal causality)
    plt.close('all')
    
    #NOTE: Remember need to add one to everywhere referring to sizes since the first element corresponds to hatch size
    
    if do_betweencorrelation:
        # Size and size changes
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Size changes vs sizes')
        for stg_1 in range(4):
            ax_h[stg_1,stg_1].scatter(all_worm_sizes[:,stg_1+1], all_worm_sizechgs[:,stg_1])
            for stg_2 in range(stg_1+1,4):
                # Above the diagnoal: size[l]->size_chg[l+1]
                ax_h[stg_1,stg_2].scatter(all_worm_sizes[:,stg_1+1], all_worm_sizechgs[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} size chgs vs. L{} sizes'.format(stg_2+1,stg_1+1))
                
                # Below the diagonal: size_chg[l]->size[l+1]
                ax_h[stg_2,stg_1].scatter(all_worm_sizechgs[:,stg_1], all_worm_sizes[:,stg_2+1])
                ax_h[stg_2,stg_1].set_title('L{} sizes vs. L{} size chgs'.format(stg_2+1,stg_1+1))
        [stats, p_vals] = scipy.stats.spearmanr(all_worm_sizes[:,1:],all_worm_sizechgs)
        stats = stats[:4,4:]
        p_vals = p_vals[:4,4:]
        print('Size (row) and size changes (columns)')
        print(stats)
        print(p_vals)
                
        # Size and stage times
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Stg durations vs sizes')
        for stg_1 in range(4):
            ax_h[stg_1,stg_1].scatter(all_worm_sizes[:,stg_1+1], all_worm_devtime[:,stg_1])
            for stg_2 in range(stg_1+1,4):
                # Above the diagnoal: size[l]->duration[l+1]
                ax_h[stg_1,stg_2].scatter(all_worm_sizes[:,stg_1+1], all_worm_devtime[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} stg duration vs. L{} sizes'.format(stg_2+1,stg_1+1))
                
                # Below the diagonal: duration[l]->size[l+1]
                ax_h[stg_2,stg_1].scatter(all_worm_devtime[:,stg_1], all_worm_sizes[:,stg_2+1])
                ax_h[stg_2,stg_1].set_title('L{} sizes vs. L{} stg duration'.format(stg_2+1,stg_1+1))
        [stats, p_vals] = scipy.stats.spearmanr(all_worm_sizes[:,1:],all_worm_devtime)
        stats = stats[:4,4:]
        p_vals = p_vals[:4,4:]
        print('Size (row) and time (columns)')
        print(stats)
        print(p_vals)
        
        # Size and average growth rates
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Avg growth rates vs sizes')
        for stg_1 in range(4):
            ax_h[stg_1,stg_1].scatter(all_worm_sizes[:,stg_1+1], all_worm_rates[:,stg_1])
            for stg_2 in range(stg_1+1,4):
                # Above the diagnoal: size[l]->rate[l+1]
                ax_h[stg_1,stg_2].scatter(all_worm_sizes[:,stg_1+1], all_worm_rates[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} rates vs. L{} sizes'.format(stg_2+1,stg_1+1))
                
                # Below the diagonal: rate[l]->size_chg[l+1]
                ax_h[stg_2,stg_1].scatter(all_worm_rates[:,stg_1], all_worm_sizes[:,stg_2+1])
                ax_h[stg_2,stg_1].set_title('L{} sizes vs. L{} rates'.format(stg_2+1,stg_1+1))
        [stats, p_vals] = scipy.stats.spearmanr(all_worm_sizes[:,1:],all_worm_rates)
        stats = stats[:4,4:]
        p_vals = p_vals[:4,4:]
        print('Size (row) and rates (columns)')
        print(stats)
        print(p_vals)
        
        # Size changes and stage times
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Stg durations vs size chgs')
        for stg_1 in range(4):
            ax_h[stg_1,stg_1].scatter(all_worm_sizechgs[:,stg_1], all_worm_devtime[:,stg_1])
            for stg_2 in range(stg_1+1,4):
                # Above the diagnoal: size_chg[l]->duration[l+1]
                ax_h[stg_1,stg_2].scatter(all_worm_sizechgs[:,stg_1], all_worm_devtime[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} stg duration vs. L{} size chg'.format(stg_2+1,stg_1+1))
                
                # Below the diagonal: size_chg[l]->size_chg[l+1]
                ax_h[stg_2,stg_1].scatter(all_worm_devtime[:,stg_1], all_worm_sizechgs[:,stg_2])
                ax_h[stg_2,stg_1].set_title('L{} size chg vs. L{} stg duration'.format(stg_2+1,stg_1+1))
        [stats, p_vals] = scipy.stats.spearmanr(all_worm_sizechgs,all_worm_devtime)
        stats = stats[:4,4:]
        p_vals = p_vals[:4,4:]
        print('Size change (row) and stage durations (columns)')
        print(stats)
        print(p_vals)
        
        # Size changes and average growth rates
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Avg growth rates vs size chgs')
        for stg_1 in range(4):
            ax_h[stg_1,stg_1].scatter(all_worm_sizechgs[:,stg_1], all_worm_rates[:,stg_1])
            for stg_2 in range(stg_1+1,4):
                # Above the diagnoal: size_chg[l]->rate[l+1]
                ax_h[stg_1,stg_2].scatter(all_worm_sizechgs[:,stg_1], all_worm_rates[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} rates vs. L{} size chg'.format(stg_2+1,stg_1+1))
                
                # Below the diagonal: rate[l]->size_chg[l+1]
                ax_h[stg_2,stg_1].scatter(all_worm_rates[:,stg_1], all_worm_sizechgs[:,stg_2])
                ax_h[stg_2,stg_1].set_title('L{} size chg vs. L{} rates'.format(stg_2+1,stg_1+1))
        [stats, p_vals] = scipy.stats.spearmanr(all_worm_sizechgs,all_worm_rates)
        stats = stats[:4,4:]
        p_vals = p_vals[:4,4:]
        print('Size changes (row) and rates (columns)')
        print(stats)
        print(p_vals)
        
        # Stage times and average growth rates
        fig_h, ax_h = plt.subplots(4,4)
        fig_h.suptitle('Avg growth rates vs stage durations')
        for stg_1 in range(4):
            ax_h[stg_1,stg_1].scatter(all_worm_devtime[:,stg_1], all_worm_rates[:,stg_1])
            for stg_2 in range(stg_1+1,4):
                # Above the diagnoal: duration[l]->rate[l+1]
                ax_h[stg_1,stg_2].scatter(all_worm_devtime[:,stg_1], all_worm_rates[:,stg_2])
                ax_h[stg_1,stg_2].set_title('L{} rates vs. L{} stg duration'.format(stg_2+1,stg_1+1))
                
                # Below the diagonal: rate[l]->duration[l+1]
                ax_h[stg_2,stg_1].scatter(all_worm_rates[:,stg_1], all_worm_devtime[:,stg_2])
                ax_h[stg_2,stg_1].set_title('L{} stg duration vs. L{} rates'.format(stg_2+1,stg_1+1))
        [stats, p_vals] = scipy.stats.spearmanr(all_worm_devtime,all_worm_rates)
        stats = stats[:4,4:]
        p_vals = p_vals[:4,4:]
        print('Stage durations (row) and rates (columns)')
        print(stats)
        print(p_vals)

        bob = input('To total correlation analysis')
    
    #######################
    # Total correlation across phenotypes
    if do_totalcorrelation:
        # Does it make sense to do this analysis...?
        plt.close('all')
        plt.gcf().clf()
        
        # Raw
        # Size and size changes
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(all_worm_sizes[:,1:-1].flatten(), all_worm_sizechgs[:,1:].flatten())
        ax_h[0].set_title('Size changes[stg+1] vs. size[stg]')
        ax_h[1].scatter(all_worm_sizes[:,1:].flatten(), all_worm_sizechgs.flatten())
        ax_h[1].set_title('Size changes[stg] vs. size[stg]')
        ax_h[2].scatter(all_worm_sizechgs[:,:-1].flatten(), all_worm_sizes[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. size_changes[stg]')
        
        # Size and average growth rate
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(all_worm_sizes[:,1:-1].flatten(), all_worm_rates[:,1:].flatten())
        ax_h[0].set_title('Growth rate[stg+1] vs. sizes[stg]')
        ax_h[1].scatter(all_worm_sizes[:,1:].flatten(), all_worm_rates.flatten())
        ax_h[1].set_title('Growth rate[stg] vs. sizes[stg]')
        ax_h[2].scatter(all_worm_rates[:,:-1].flatten(), all_worm_sizes[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. rates[stg]')
        
        # Size changes and average growth rate
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(all_worm_sizechgs[:,:-1].flatten(), all_worm_rates[:,1:].flatten())
        ax_h[0].set_title('Growth rate[stg+1] vs. size changes[stg]')
        ax_h[1].scatter(all_worm_sizechgs.flatten(), all_worm_rates.flatten())
        ax_h[1].set_title('Growth rate[stg] vs. size changes[stg]')
        ax_h[2].scatter(all_worm_rates[:,:-1].flatten(), all_worm_sizechgs[:,1:].flatten())
        ax_h[2].set_title('Size changes[stg+1] vs. rates[stg]')
        
        # Rate and timing
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(all_worm_devtime[:,:-1].flatten(), all_worm_rates[:,1:].flatten())
        ax_h[0].set_title('Rate[stg+1] vs. Timing[stg]')
        ax_h[1].scatter(all_worm_rates.flatten(), all_worm_devtime.flatten())
        ax_h[1].set_title('Timing[stg] vs. rate[stg]')
        ax_h[2].scatter(all_worm_rates[:,:-1].flatten(), all_worm_devtime[:,1:].flatten())
        ax_h[2].set_title('Timing [stg+1] vs. Rate[stg]')
        
        # Size and RELATIVE size changes
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(all_worm_sizes[:,1:-1].flatten(), sizechgs_relative[:,1:].flatten())
        ax_h[0].set_title('Relative size changes[stg+1] vs. size[stg]')
        ax_h[1].scatter(all_worm_sizes[:,1:].flatten(), sizechgs_relative.flatten())
        ax_h[1].set_title('Relative size changes[stg] vs. size[stg]')
        ax_h[2].scatter(sizechgs_relative[:,:-1].flatten(), all_worm_sizes[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. relative size_changes[stg]')
        
        # Size and RELATIVE average growth rate
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(all_worm_sizes[:,1:-1].flatten(), rates_relative[:,1:].flatten())
        ax_h[0].set_title('Relative growth rate[stg+1] vs. sizes[stg]')
        ax_h[1].scatter(all_worm_sizes[:,1:].flatten(), rates_relative.flatten())
        ax_h[1].set_title('Relative rowth rate[stg] vs. sizes[stg]')
        ax_h[2].scatter(rates_relative[:,:-1].flatten(), all_worm_sizes[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. relative growth rates[stg]')
        
        '''
        bob = input('Mean-centered')
        
        # Mean-centered
        plt.close('all')
        plt.gcf().clf()
        
        # Size and size changes
        fig_h, ax_h = plt.subplots(1,2)
        ax_h[0].scatter((all_worm_sizes[:,1:-1]-np.mean(all_worm_sizes[:,1:-1],axis=0)).flatten(), (all_worm_sizechgs[:,1:]-np.mean(all_worm_sizechgs[:,1:],axis=0)).flatten())
        ax_h[0].set_title('Size changes[stg+1] vs. size[stg]')
        ax_h[1].scatter((all_worm_sizechgs[:,:-1]-np.mean(all_worm_sizechgs[:,:-1],axis=0)).flatten(), (all_worm_sizes[:,2:]-np.mean(all_worm_sizes[:,2:],axis=0)).flatten())
        ax_h[1].set_title('Size[stg+1] vs. size_changes[stg]')
        
        # Size and average growth rate
        fig_h, ax_h = plt.subplots(1,2)
        ax_h[0].scatter((all_worm_sizes[:,1:-1]-np.mean(all_worm_sizes[:,1:-1],axis=0)).flatten(), (all_worm_rates[:,1:]-np.mean(all_worm_rates[:,1:],axis=0)).flatten())
        ax_h[0].set_title('Growth rate[stg+1] vs. sizes[stg]')
        ax_h[1].scatter((all_worm_rates[:,:-1]-np.mean(all_worm_rates[:,:-1],axis=0)).flatten(), (all_worm_sizes[:,2:]-np.mean(all_worm_sizes[:,2:],axis=0)).flatten())
        ax_h[1].set_title('Size[stg+1] vs. rates[stg]')
        
        # Size changes and average growth rate
        fig_h, ax_h = plt.subplots(1,2)
        ax_h[0].scatter((all_worm_sizechgs[:,:-1]-np.mean(all_worm_sizechgs[:,:-1],axis=0)).flatten(), (all_worm_rates[:,1:]-np.mean(all_worm_rates[:,1:],axis=0)).flatten())
        ax_h[0].set_title('Growth rate[stg+1] vs. size changes[stg]')
        ax_h[1].scatter((all_worm_rates[:,:-1]-np.mean(all_worm_rates[:,:-1],axis=0)).flatten(), (all_worm_sizechgs[:,1:]-np.mean(all_worm_sizechgs[:,1:],axis=0)).flatten())
        ax_h[1].set_title('Size changes[stg+1] vs. rates[stg]')
        '''
        
        ob = input('Z-scored')
        
        # Z-scored
        plt.close('all')
        plt.gcf().clf()        
        
        # Size and size changes
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(sizes_zscored[:,1:-1].flatten(), sizechgs_zscored[:,1:].flatten())
        (reg_slope,reg_int,reg_r,reg_p,reg_se) = scipy.stats.linregress(sizes_zscored[:,1:-1].flatten(), sizechgs_zscored[:,1:].flatten())
        ax_h[0].set_title('Size changes[stg+1] vs. size[stg]\n sizechg[stg+1] = {:1.3f}*size[stg]+{:1.3f}'.format(reg_slope,reg_int))
        ax_h[1].scatter(sizes_zscored[:,1:].flatten(), sizechgs_zscored.flatten())
        (reg_slope,reg_int,reg_r,reg_p,reg_se) = scipy.stats.linregress(sizes_zscored[:,1:].flatten(), sizechgs_zscored.flatten())
        ax_h[1].set_title('Size changes[stg] vs. size[stg]\n sizechg[stg] = {:1.3f}*size[stg]+{:1.3f}'.format(reg_slope,reg_int))
        ax_h[2].scatter(sizechgs_zscored[:,:-1].flatten(), sizes_zscored[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. size_changes[stg]')
        
        # Size and average growth rate
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(sizes_zscored[:,1:-1].flatten(), rates_zscored[:,1:].flatten())
        ax_h[0].set_title('Growth rate[stg+1] vs. sizes[stg]')
        ax_h[1].scatter(sizes_zscored[:,1:].flatten(), rates_zscored.flatten())
        ax_h[1].set_title('Growth rate[stg] vs. sizes[stg]')
        ax_h[2].scatter(rates_zscored[:,:-1].flatten(), sizes_zscored[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. rates[stg]')
        
        # Size changes and average growth rate
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(sizechgs_zscored[:,:-1].flatten(), rates_zscored[:,1:].flatten())
        ax_h[0].set_title('Growth rate[stg+1] vs. size changes[stg]')
        ax_h[1].scatter(sizechgs_zscored.flatten(), rates_zscored.flatten())
        ax_h[1].set_title('Growth rate[stg] vs. size changes[stg]')
        ax_h[2].scatter(rates_zscored[:,:-1].flatten(), sizechgs_zscored[:,1:].flatten())
        ax_h[2].set_title('Size changes[stg+1] vs. rates[stg]')
        
        # Size and timing
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(sizes_zscored[:,1:-1].flatten(), devtime_zscored[:,1:].flatten())
        ax_h[0].set_title('Timing[stg+1] vs. size[stg]')
        ax_h[1].scatter(sizes_zscored[:,1:].flatten(), devtime_zscored.flatten())
        ax_h[1].set_title('Timing[stg] vs. size[stg]')
        ax_h[2].scatter(devtime_zscored[:,:-1].flatten(), sizes_zscored[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. Timing[stg]')
        
        # Rate and timing
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(devtime_zscored[:,:-1].flatten(), rates_zscored[:,1:].flatten())
        ax_h[0].set_title('Rate[stg+1] vs. Timing[stg]')
        ax_h[1].scatter(rates_zscored.flatten(), devtime_zscored.flatten())
        ax_h[1].set_title('Timing[stg] vs. rate[stg]')
        ax_h[2].scatter(rates_zscored[:,:-1].flatten(), devtime_zscored[:,1:].flatten())
        ax_h[2].set_title('Timing [stg+1] vs. Rate[stg]')
        
        # Size changes and timing
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(devtime_zscored[:,:-1].flatten(), sizechgs_zscored[:,1:].flatten())
        ax_h[0].set_title('Size change [stg+1] vs. Timing[stg]')
        ax_h[1].scatter(sizechgs_zscored.flatten(), devtime_zscored.flatten())
        ax_h[1].set_title('Timing[stg] vs. Size change [stg]')
        ax_h[2].scatter(sizechgs_zscored[:,:-1].flatten(), devtime_zscored[:,1:].flatten())
        ax_h[2].set_title('Timing [stg+1] vs. Size changes [stg]')
        
        #Size vs. size at later time
        #fig_h = plt.figure()
        #ax_h = fig_h.gca()
        #ax_h.scatter(sizes_zscored[:,1:-1].flatten(), sizes_zscored[:,2:].flatten())
        #ax_h.set_title('Sizes [stg+1] vs. Sizes[stg]')
        fig_h, ax_h = plt.subplots(1,3)
        fig_h.suptitle('Autocorrelations in size over time')
        for lag in range(3):
            [ax_h[lag].scatter(sizes_zscored[:,stg+1],sizes_zscored[:,stg+lag+1+1]) for stg in range(4-lag-1)]
            ax_h[lag].set_title('Stages {} away'.format(lag+1))
            ax_h[lag].set_xlabel('Size[stg] (z-score)')
            ax_h[lag].set_ylabel('Size[stg+lag] (z-score)')

        
        #Size change vs. size change at later time
        #fig_h = plt.figure()
        #ax_h = fig_h.gca()
        #ax_h.scatter(sizechgs_zscored[:,:-1].flatten(),sizechgs_zscored[:,1:].flatten())
        #ax_h.set_title('Size changes [stg+1] vs. size changes[stg]')
        fig_h, ax_h = plt.subplots(1,3)
        fig_h.suptitle('Autocorrelations in size change over time')
        for lag in range(3):
            [ax_h[lag].scatter(sizechgs_zscored[:,stg],sizechgs_zscored[:,stg+lag+1]) for stg in range(4-lag-1)]
            print(np.concatenate([sizechgs_zscored[:,stg] for stg in range(4-lag-1)]).shape)
            corr_val, p_val = scipy.stats.spearmanr(
                np.concatenate([sizechgs_zscored[:,stg] for stg in range(4-lag-1)]),
                np.concatenate([sizechgs_zscored[:,stg+lag+1] for stg in range(4-lag-1)]))
            
            ax_h[lag].set_title('Stages {} away \n Spearman r={:.2f}, p={:.3f}'.format(lag+1,corr_val,p_val))
            ax_h[lag].set_xlabel('Size change [stg] (z-score)')
            ax_h[lag].set_ylabel('Size change [stg+lag] (z-score)')
            
        # Growth rate autocorrelations
        #fig_h = plt.figure()
        #ax_h = fig_h.gca()
        #ax_h.scatter(sizechgs_zscored[:,:-1].flatten(),sizechgs_zscored[:,1:].flatten())
        #ax_h.set_title('Size changes [stg+1] vs. size changes[stg]')
        fig_h, ax_h = plt.subplots(1,3)
        fig_h.suptitle('Autocorrelations in growth rate over time')
        for lag in range(3):
            [ax_h[lag].scatter(rates_zscored[:,stg],rates_zscored[:,stg+lag+1]) for stg in range(4-lag-1)]
            #print(np.concatenate([sizechgs_zscored[:,stg] for stg in range(4-lag-1)]).shape)
            corr_val, p_val = scipy.stats.spearmanr(
                np.concatenate([rates_zscored[:,stg] for stg in range(4-lag-1)]),
                np.concatenate([rates_zscored[:,stg+lag+1] for stg in range(4-lag-1)]))
            
            ax_h[lag].set_title('Stages {} away \n Spearman r={:.2f}, p={:.3f}'.format(lag+1,corr_val,p_val))
            ax_h[lag].set_xlabel('Growth rate [stg] (z-score)')
            ax_h[lag].set_ylabel('Growth rate [stg+lag] (z-score)')
        
        #Stage duration vs stage duration at later time
        #fig_h = plt.figure()
        #ax_h = fig_h.gca()
        #ax_h.scatter(sizechgs_zscored[:,:-1].flatten(),sizechgs_zscored[:,1:].flatten())
        #ax_h.set_title('Size changes [stg+1] vs. size changes[stg]')
        fig_h, ax_h = plt.subplots(1,3)
        fig_h.suptitle('Autocorrelations in larval stage duration over time')
        for lag in range(3):
            [ax_h[lag].scatter(devtime_zscored[:,stg],devtime_zscored[:,stg+lag+1]) for stg in range(4-lag-1)]
            #print(np.concatenate([sizechgs_zscored[:,stg] for stg in range(4-lag-1)]).shape)
            corr_val, p_val = scipy.stats.spearmanr(
                np.concatenate([devtime_zscored[:,stg] for stg in range(4-lag-1)]),
                np.concatenate([devtime_zscored[:,stg+lag+1] for stg in range(4-lag-1)]))
            
            ax_h[lag].set_title('Stages {} away \n Spearman r={:.2f}, p={:.3f}'.format(lag+1,corr_val,p_val))
            ax_h[lag].set_xlabel('Stage duration [stg] (z-score)')
            ax_h[lag].set_ylabel('Stage duration [stg+lag] (z-score)')
        
        # Size and relative size changes
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(sizes_zscored[:,1:-1].flatten(), sizechgs_relative_zscored[:,1:].flatten())
        (reg_slope,reg_int,reg_r,reg_p,reg_se) = scipy.stats.linregress(sizes_zscored[:,1:-1].flatten(), sizechgs_relative_zscored[:,1:].flatten())
        ax_h[0].set_title('Relative size changes[stg+1] vs. size[stg]\n sizechg[stg+1] = {:1.3f}*size[stg]+{:1.3f}'.format(reg_slope,reg_int))
        ax_h[1].scatter(sizes_zscored[:,1:].flatten(), sizechgs_relative_zscored.flatten())
        (reg_slope,reg_int,reg_r,reg_p,reg_se) = scipy.stats.linregress(sizes_zscored[:,1:].flatten(), sizechgs_relative_zscored.flatten())
        ax_h[1].set_title('Relative size changes[stg] vs. size[stg]\n sizechg[stg] = {:1.3f}*size[stg]+{:1.3f}'.format(reg_slope,reg_int))
        ax_h[2].scatter(sizechgs_relative_zscored[:,:-1].flatten(), sizes_zscored[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. relative size_changes[stg]')
        
        # Size and RELATIVE average growth rate
        fig_h, ax_h = plt.subplots(1,3)
        ax_h[0].scatter(sizes_zscored[:,1:-1].flatten(), rates_relative_zscored[:,1:].flatten())
        ax_h[0].set_title('Relative rowth rate[stg+1] vs. sizes[stg]')
        ax_h[1].scatter(sizes_zscored[:,1:].flatten(), rates_relative_zscored.flatten())
        ax_h[1].set_title('Relatie growth rate[stg] vs. sizes[stg]')
        ax_h[2].scatter(rates_relative_zscored[:,:-1].flatten(), sizes_zscored[:,2:].flatten())
        ax_h[2].set_title('Size[stg+1] vs. relative growth rates[stg]')
        
        # RELATIVE size changes and timing
        fig_h, ax_h = plt.subplots(1,3) 
        ax_h[0].scatter(devtime_zscored[:,:-1].flatten(), sizechgs_relative_zscored[:,1:].flatten())
        ax_h[0].set_title('Relative size change [stg+1] vs. Timing[stg]')
        ax_h[1].scatter(sizechgs_relative_zscored.flatten(), devtime_zscored.flatten())
        ax_h[1].set_title('Timing[stg] vs. Relative size change [stg]')
        ax_h[2].scatter(sizechgs_relative_zscored[:,:-1].flatten(), devtime_zscored[:,1:].flatten())
        ax_h[2].set_title('Timing [stg+1] vs. Relative size changes [stg]')
        
        #Relative size change vs. relative size change at later time
        fig_h, ax_h = plt.subplots(1,3)
        fig_h.suptitle('Autocorrelations in relative size change over time')
        for lag in range(3):
            [ax_h[lag].scatter(sizechgs_relative_zscored[:,stg],sizechgs_relative_zscored[:,stg+lag+1]) for stg in range(4-lag-1)]
            #print(np.concatenate([sizechgs_zscored[:,stg] for stg in range(4-lag-1)]).shape)
            corr_val, p_val = scipy.stats.spearmanr(
                np.concatenate([sizechgs_relative_zscored[:,stg] for stg in range(4-lag-1)]),
                np.concatenate([sizechgs_relative_zscored[:,stg+lag+1] for stg in range(4-lag-1)]))
            
            ax_h[lag].set_title('Stages {} away \n Spearman r={:.2f}, p={:.3f}'.format(lag+1,corr_val,p_val))
            ax_h[lag].set_xlabel('Relative size change [stg]')
            ax_h[lag].set_ylabel('Relative size change [stg+lag]')

        
        print('Correlations - Z-scored data')
        #(corr_vals, p_vals) = scipy.stats.spearmanr(
            #np.stack((sizes_zscored[:,1:-1].flatten(), rates_zscored[:,:-1].flatten(), sizechgs_zscored[:,:-1].flatten()), axis=1), # Before
            #np.stack((sizes_zscored[:,2:].flatten(), rates_zscored[:,1:].flatten(), sizechgs_zscored[:,1:].flatten()), axis=1))  # After
        print('Within stages')
        print('sizes, size chgs, rates, timing')
        (corr_vals, p_vals) = scipy.stats.spearmanr(
            np.stack((sizes_zscored[:,1:].flatten(),sizechgs_zscored.flatten(),rates_zscored.flatten(),devtime_zscored.flatten()),axis=1))
        print(corr_vals)
        print(p_vals)
        print('Across stages')
        print('sizes[t], size chgs[t], rates[t],timing[t] , sizes [t+1], size chgs [t+1], rates[t+1],timing[t+1]')
        (corr_vals, p_vals) = scipy.stats.spearmanr(
            np.stack((sizes_zscored[:,1:-1].flatten(), sizechgs_zscored[:,:-1].flatten(), rates_zscored[:,:-1].flatten(), devtime_zscored[:,:-1].flatten()), axis=1), # Before
            np.stack((sizes_zscored[:,2:].flatten(), sizechgs_zscored[:,1:].flatten(), rates_zscored[:,1:].flatten(), devtime_zscored[:,1:].flatten()), axis=1))  # After
        print(corr_vals)
        print(p_vals)
        print('r_partial for sizes[t]->sizechgs[t+1]')
        print(partialcoef(sizes_zscored[:,1:-1].flatten(),sizechgs_zscored[:,1:].flatten(),sizechgs_zscored[:,:-1].flatten()))
        print('r_partial for sizechgs[t]->sizechgs[t+1]')
        print(partialcoef(sizechgs_zscored[:,:-1].flatten(),sizechgs_zscored[:,1:].flatten(),sizes_zscored[:,1:-1].flatten()))
        
        bob = input('to trajectories')
        
        ###############################################
        # Trajectories
        plt.close('all')
        plt.gcf().clf()
        
        # Size
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        ax_h.set_title('Size trajectories')
        ax_h.set_xlabel('Size[stg]')
        ax_h.set_ylabel('Size[stg+1]')
        for worm_num in range(len(sizes_zscored)):
            #ax_h.plot(sizes_zscored[worm_num,1:-1],sizes_zscored[worm_num,2:])
            #ax_h.scatter(sizes_zscored[worm_num,1],sizes_zscored[worm_num,2], c='g')
            #ax_h.scatter(sizes_zscored[worm_num,-2],sizes_zscored[worm_num,-1], c='r')
            #bob=input('waiting')
            
            # Raw with hatching
            #ax_h.plot(all_worm_sizes[worm_num,0:-1],
                #all_worm_sizes[worm_num,1:],'o-')
            
            # Raw w/o hatching
            ax_h.plot(all_worm_sizes[worm_num,1:-1],
                all_worm_sizes[worm_num,2:],'o-')
            ax_h.scatter(all_worm_sizes[worm_num,-2],
                all_worm_sizes[worm_num,-1], c='r')
                
            # All sizes aligned at beginning (i.e. subtract [size[stg],size[stg+1]] from all values) w/o hatching
            ax_h.plot(all_worm_sizes[worm_num,1:-1]-(all_worm_sizes[worm_num,1]),
                all_worm_sizes[worm_num,2:]-(all_worm_sizes[worm_num,2]),'o-')
            ax_h.scatter(all_worm_sizes[worm_num,-2]-(all_worm_sizes[worm_num,1]),
                all_worm_sizes[worm_num,-1]-(all_worm_sizes[worm_num,2]), c='r')
            
            # Z-scored unaligned w/o hatching
            #ax_h.plot(all_worm_sizes[worm_num,1:-1],all_worm_sizes[worm_num,2:])
            #ax_h.scatter(all_worm_sizes[worm_num,1],all_worm_sizes[worm_num,2], c='g')
            #ax_h.scatter(all_worm_sizes[worm_num,-2],all_worm_sizes[worm_num,-1], c='r')
            
        # Vector autoregression on sizes 
        #def format_autoregression_sizes(data):
            #predictor_mat = []
            #realized_mat = []
            #for sample in data:
                #print(sample)
                #for item_num in range(len(sample)-2):
                    #predictor_mat.append([sample[item_num],sample[item_num+1]])
                    #realized_mat.append([sample[item_num+1],sample[item_num+2]])
            #return (np.array(predictor_mat), np.array(realized_mat))
        #p_mat, r_mat = format_autoregression(all_worm_sizes[:,1:])
        #glmf = sklearn.linear_model.LinearRegression()
        #glmf.fit(p_mat, r_mat)
        #prediction = glmf.predict(p_mat)
        #ax_h.scatter(prediction[:,0],prediction[:,1],marker="1")    # Upside-down triangle
        
        # Bootstrap for dispersion in distribution
        num_realizations = 10000
        num_samples = all_worm_sizechgs.shape[0]
        iqr = []
        idr = []
        for idx in range(num_realizations):
            bs_sizes = np.hstack(  # Starting from hatch
                (sizechgs[np.random.randint(all_worm_sizechgs.shape[0],size=(num_samples,1))] for sizechgs in all_worm_sizechgs.T)) \
                + all_worm_sizes[np.random.randint(all_worm_sizechgs.shape[0],size=(num_samples,1)),0]
            iqr.append(np.subtract(*np.percentile(bs_sizes, [75, 25],axis=0,interpolation='higher')))
            idr.append(np.subtract(*np.percentile(bs_sizes, [90, 10],axis=0,interpolation='higher')))
        iqr = np.array(iqr)
        idr = np.array(idr)
        
        data_iqr = np.subtract(*np.percentile(all_worm_sizes[:,1:], [75, 25],axis=0,interpolation='higher'))
        data_idr = np.subtract(*np.percentile(all_worm_sizes[:,1:], [90, 10],axis=0,interpolation='higher'))
        
        p_iqr = np.array([np.count_nonzero(stg_iqr<stg_data_iqr)/stg_iqr.size for stg_iqr,stg_data_iqr in zip(iqr.T,data_iqr)])
        print(p_iqr)
        p_idr = np.array([np.count_nonzero(stg_idr<stg_data_idr)/stg_idr.size for stg_idr,stg_data_idr in zip(idr.T,data_idr)])
        print(p_idr)
            
        # Size changes
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        ax_h.set_title('Size change trajectories')
        ax_h.set_xlabel('Size change[stg]')
        ax_h.set_ylabel('Size change[stg+1]')
        for worm_num in range(len(sizes_zscored)):
            #ax_h.plot(sizes_zscored[worm_num,1:-1],sizes_zscored[worm_num,2:])
            #ax_h.scatter(sizes_zscored[worm_num,1],sizes_zscored[worm_num,2], c='g')
            #ax_h.scatter(sizes_zscored[worm_num,-2],sizes_zscored[worm_num,-1], c='r')
            #bob=input('waiting')
            
            # Raw with hatching
            #ax_h.plot(all_worm_sizes[worm_num,0:-1],
                #all_worm_sizes[worm_num,1:],'o-')
            
            # Raw w/o hatching
            #ax_h.plot(all_worm_sizechgs[worm_num,1:-1],
                #all_worm_sizechgs[worm_num,2:],'o-')
            #ax_h.scatter(all_worm_sizechgs[worm_num,-2],
                #all_worm_sizechgs[worm_num,-1], c='r')
                
            # All sizes aligned at beginning (i.e. subtract [size[stg],size[stg+1]] from all values) w/o hatching
            #ax_h.plot(all_worm_sizechgs[worm_num,:-1]-(all_worm_sizechgs[worm_num,0]),
                #all_worm_sizechgs[worm_num,1:]-(all_worm_sizechgs[worm_num,1]),'o-')
            #ax_h.scatter(all_worm_sizechgs[worm_num,-2]-(all_worm_sizechgs[worm_num,0]),
                #all_worm_sizechgs[worm_num,-1]-(all_worm_sizechgs[worm_num,1]), c='r')
            
            # Z-scored unaligned w/o hatching
            ax_h.plot(sizechgs_zscored[worm_num,:-1],sizechgs_zscored[worm_num,1:])
            ax_h.scatter(sizechgs_zscored[worm_num,0],sizechgs_zscored[worm_num,1], c='g')
            ax_h.scatter(sizechgs_zscored[worm_num,-2],sizechgs_zscored[worm_num,-1], c='r')
        
        
        # Size and size changes
        fig_h = plt.figure()
        ax_h = fig_h.gca()
        ax_h.set_title('Size+size change trajectories')
        ax_h.set_xlabel('Size[stg]')
        ax_h.set_ylabel('Size change[stg]]')
        for worm_num in range(len(sizes_zscored)):       
        #for worm_num in range(8):       
            # Raw w/o hatching
            #ax_h.plot(all_worm_sizes[worm_num,1:],
                #all_worm_sizechgs[worm_num,:],'o-')
            #ax_h.scatter(all_worm_sizes[worm_num,-1],
                #all_worm_sizechgs[worm_num,-1], c='r')
                
            # All sizes aligned at beginning (i.e. subtract [size[stg],size[stg+1]] from all values) w/o hatching
            #ax_h.plot(all_worm_sizes[worm_num,1:]-(all_worm_sizes[worm_num,1]),
                #all_worm_sizechgs[worm_num,:]-(all_worm_sizechgs[worm_num,0]),'o-')
            #ax_h.scatter(all_worm_sizes[worm_num,-1]-(all_worm_sizes[worm_num,1]),
                #all_worm_sizechgs[worm_num,-1]-(all_worm_sizechgs[worm_num,0]), c='r')
                
            # All sizes aligned with padding
            ax_h.plot(np.pad(all_worm_sizes[worm_num,1:],(1,0),mode='constant',constant_values=0),
                np.pad(all_worm_sizechgs[worm_num,:],(1,0),mode='constant',constant_values=0),'o-')
            ax_h.scatter(0,0, c='g')
            ax_h.scatter(all_worm_sizes[worm_num,-1],
                all_worm_sizechgs[worm_num,-1], c='r')
            
            # Z-scored unaligned w/o hatching
            #ax_h.plot(sizes_zscored[worm_num,1:],sizechgs_zscored[worm_num,:])
            #ax_h.scatter(sizes_zscored[worm_num,1],sizechgs_zscored[worm_num,0], c='g')
            #ax_h.scatter(sizes_zscored[worm_num,-1],sizechgs_zscored[worm_num,-1], c='r')
            
            # Z-scored aligned
            #ax_h.plot(sizes_zscored[worm_num,1:]-sizes_zscored[worm_num,1],sizechgs_zscored[worm_num,:]-sizechgs_zscored[worm_num,0])
            #ax_h.scatter(sizes_zscored[worm_num,1]-sizes_zscored[worm_num,1],sizechgs_zscored[worm_num,0]-sizechgs_zscored[worm_num,0], c='g')
            #ax_h.scatter(sizes_zscored[worm_num,-1]-sizes_zscored[worm_num,1],sizechgs_zscored[worm_num,-1]-sizechgs_zscored[worm_num,0], c='r')
            
            # Z-scored aligned with padding
            ax_h.plot(np.pad(sizes_zscored[worm_num,1:],(1,0),mode='constant',constant_values=0),    \
                np.pad(sizechgs_zscored[worm_num,:],(1,0),mode='constant',constant_values=0))
            ax_h.scatter(0,0, c='g')
            ax_h.scatter(sizes_zscored[worm_num,-1],sizechgs_zscored[worm_num,-1], c='r')
        
        '''
        def format_autoregression_sizesandsizechgs(sizes,sizechgs):
            predictor_mat = []
            realized_mat = []
            for stg_size, stg_sizechg in zip(sizes,sizechgs):
                for size_item, sizechg_item in zip(stg_size,stg_sizechg):
                    predictor_mat.append([size_item,sizechg_item])
                    realized_mat.append([size_item,sizechg_item])
            return (np.array(predictor_mat), np.array(realized_mat))
        p_mat, r_mat = format_autoregression(all_worm_sizes[:,1:])
        glmf = sklearn.linear_model.LinearRegression()
        glmf.fit(p_mat, r_mat)
        prediction = glmf.predict(p_mat)
        ax_h.scatter(prediction[:,0],prediction[:,1],marker="1")    # Upside-down triangle
        '''
        
        
        bob = input('To modeling')
        
    bob = input('To lifespan analysis')

    ###############
    # Load lifespans
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

    # For each experiment, grab sizes and rates of growth for worms, as well as expt_times
    compiled_expts = [process_data.Experiment(expt_path) for expt_path in expt_paths]
    lifespan_allexpts = [[] for idx in range(len(expt_paths))]
    goodworms_allexpts = []
    
    for expt_num, [expt, expt_events] in enumerate(zip(compiled_expts, timestamped_data)):
        acq_worms = expt.get_acquired_wells()
        first_worm_num = int(expt_events['Worm'][0][-2:])   # Add in adjustment for the first index of the worms not being at 0
        viable_worm = (expt_events['Hatch']!=-1) \
            & (expt_events['Death']!=-1) \
            & np.array([not any([kw in note for kw in bad_worm_kws]) for note in expt_events['Notes']])
        #& ((expt_events['L1 ecdysis']-expt_events['Hatch'])/3600 <= 15) \
        hatched_on_corral = timestamped_data[expt_num]['Hatch'] != 0
        worm_was_acquired = [str(worm_num+first_worm_num).zfill(len(expt.get_acquired_wells()[0])) in expt.get_acquired_wells() for worm_num in np.arange(len(expt_events['Hatch']))]

        print(np.count_nonzero(viable_worm & hatched_on_corral & worm_was_acquired))
        
        # Calculate lifespan for worms
        lifespan_allexpts[expt_num] = (expt_events['Death']-expt_events['Hatch'])[viable_worm & hatched_on_corral & worm_was_acquired]/3600
        goodworms_allexpts.append([worm_name.replace('/','') for worm_name in expt_events['Worm_FullName'][viable_worm & hatched_on_corral & worm_was_acquired]])
    
    lifespans = np.concatenate(lifespan_allexpts,axis=0)   
    goodworms = [worm for good_group in goodworms_allexpts for worm in good_group]
    
    # Now scatter data against lifespan
    plt.close('all')
    
    # Plot ending size against lifespan
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('Size at end of stage')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(all_worm_sizes[:,stg+1],lifespans)
        ax_h[stg].set_xlabel('Size (mm^2)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
    
    # Plot size change per stage against lifespan
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('Size change per stage')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(all_worm_sizechgs[:,stg],lifespans)
        ax_h[stg].set_xlabel('Size change (mm^2)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
    
    # Plot duration of larval stages against lifespan
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('Duration of larval stage')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(all_worm_devtime[:,stg],lifespans)
        ax_h[stg].set_xlabel('Larval stage duration (hr)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
    
    # Plot average growth rate over stage
    fig_h, ax_h = plt.subplots(4,1)
    fig_h.suptitle('Average growth rate across stages')
    for stg in range(4):
        ax_h[stg].set_title('L{}'.format(stg+1))
        ax_h[stg].scatter(all_worm_rates[:,stg],lifespans)
        ax_h[stg].set_xlabel('Growth rate (mm^2/hr)')
        ax_h[stg].set_ylabel('Lifespan (hr.)')
        ax_h[stg].set_xlim([-0.001,0.003])
    
    # Lifespan vs. Variation in Growth
    
    # Size
    fig_h= plt.figure()
    ax_h = fig_h.gca()
    ax_h.set_title('Lifespan vs. Variation in Size')
    ax_h.scatter(np.sum(np.abs(sizes_zscored[:,1:]),axis=1),lifespans)
    (corr_val, p_val) = scipy.stats.spearmanr(np.sum(np.abs(sizes_zscored[:,1:]),axis=1),lifespans)
    print(corr_val)
    print(p_val)
    
    # Size changes
    fig_h= plt.figure()
    ax_h = fig_h.gca()
    ax_h.set_title('Lifespan vs. Variation in Size Change')
    ax_h.scatter(np.sum(np.abs(sizechgs_zscored),axis=1),lifespans)
    (corr_val, p_val) = scipy.stats.spearmanr(np.sum(np.abs(sizechgs_zscored),axis=1),lifespans)
    print(corr_val)
    print(p_val)
    
    # Timing
    fig_h= plt.figure()
    ax_h = fig_h.gca()
    ax_h.set_title('Lifespan vs. Variation in Timing')
    ax_h.scatter(np.sum(np.abs(devtime_zscored),axis=1),lifespans)
    (corr_val, p_val) = scipy.stats.spearmanr(np.sum(np.abs(devtime_zscored),axis=1),lifespans)
    print(corr_val)
    print(p_val)
    
    # Rates
    fig_h= plt.figure()
    ax_h = fig_h.gca()
    ax_h.set_title('Lifespan vs. Variation in Growth Rate')
    ax_h.scatter(np.sum(np.abs(rates_zscored),axis=1),lifespans)
    (corr_val, p_val) = scipy.stats.spearmanr(np.sum(np.abs(rates_zscored),axis=1),lifespans)
    print(corr_val)
    print(p_val)
        
    bob = input('To healthspan analysis')
    
    #############################
    
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
    
    # Use autofluorescence(intensity_90), bulk_movement, stimulated_rate as regressors
    health_params = ['intensity_90', 'bulk_movement', 'stimulated_rate_b', 'total_size']
    
    # Intraindividual specific cutoffs
    cutoffs_frac = [0.5, 0.8, 0.95,0.9]    # Fixed cutoff as percentile of range of parameter
    worm_intra_events = np.zeros((len(worm_health),len(health_params)))
    #for worm_num,worm in enumerate(worm_health):
        #last_timepoint = np.where(worm['intensity_90']==-1.0)[0][0] if any(worm['intensity_90']==-1.0) else len(worm['intensity_90'])
        #worm_intra_events[worm_num] = np.array([
            #worm['time'][np.where(worm['intensity_90']<cutoffs_frac[0]*np.ptp(worm['intensity_90'][:last_timepoint]))[0][0]],
            #worm['time'][np.where(
                #(worm['bulk_movement']<cutoffs_frac[1]*np.ptp(worm['bulk_movement'][:last_timepoint])) \
                #& (worm['time']>worm['time'][np.where(worm['bulk_movement']==worm['bulk_movement'].max())[0][0]]))[0][0]],
            #worm['time'][np.where(
                #(worm['stimulated_rate_b']<cutoffs_frac[1]*np.ptp(worm['stimulated_rate_b'][:last_timepoint])) \
                #& (worm['time']>worm['time'][np.where(worm['stimulated_rate_b']==worm['stimulated_rate_b'].max())[0][0]]))[0][0]] \
        #])
        
    for worm_num, worm in enumerate(worm_health):
        last_timepoint = np.where(worm['intensity_90']==-1.0)[0][0] if any(worm['intensity_90']==-1.0) else len(worm['intensity_90'])
        for param_num, phys_param in enumerate(health_params):
            worm_max = worm[phys_param].max()
            worm_range = np.ptp(worm[phys_param][:last_timepoint])
            if phys_param != 'intensity_90':    # Decreasing over time
                worm_intra_events[worm_num,param_num] = worm['time'][np.where(
                    (worm[phys_param]<(worm_max-cutoffs_frac[param_num]*worm_range))\
                    & (worm['time']>worm['time'][np.where(worm[phys_param]==worm[phys_param].max())[0][0]]))[0][0]]
            else:
                worm_intra_events[worm_num,param_num] = worm['time'][np.where(
                    worm[phys_param]>(worm_max-cutoffs_frac[param_num]*worm_range))[0][0]]
        
    worm_intra_events = worm_intra_events * 24 # Convert to hours
        
        
    # Scatter and calculate correlations
    size_xmin = [0,0.005,0.01,0.015]
    size_xmax = [0.01,0.015,0.025,0.05]
    sizechg_xmax = [0.01,0.01,0.015,0.0305]
    rate_xmax = [0.0005, 0.001, 0.002, 0.003]
    
    plt.close('all')
    
    # Plot ending size against healthspans
    fig_h, ax_h = plt.subplots(4,len(health_params))
    fig_h.suptitle('Size at end of stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_intra_events.T,health_params)):
        for stg in range(4):
            ax_h[stg,event_num].scatter(all_worm_sizes[:,stg+1],health_event)
            if stg==3: ax_h[stg,event_num].set_xlabel('Size (mm^2)')
            #ax_h[stg,event_num].set_ylabel(health_param)
    print('Correlations - size vs. intraindividual health')
    
    #(corr_vals, p_vals) = scipy.stats.spearmanr(
        #np.concatenate((all_worm_sizes, worm_inter_events), axis = 1))
    #corr_vals = corr_vals[0:5,5:8]
    #p_vals = p_vals[0:5,5:8]
    #print(corr_vals)
    #print(p_vals)
    '''
    # Plot size change per stage against healthspans
    fig_h, ax_h = plt.subplots(4,len(health_params))
    fig_h.suptitle('Rate')
    for event_num, (health_event,health_param) in enumerate(zip(worm_intra_events.T,health_params)):
        for stg in range(4):
            ax_h[stg,event_num].scatter(all_worm_rates[:,stg],health_event)
            if stg == 3: ax_h[stg,event_num].set_xlabel('Growth Rate (mm^2/hr)')
            #ax_h[stg,event_num].set_ylabel(health_param)
    
    # Plot size change per stage against healthspans
    fig_h, ax_h = plt.subplots(4,len(health_params))
    fig_h.suptitle('Size change per stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_intra_events.T,health_params)):
        for stg in range(4):
            ax_h[stg,event_num].scatter(all_worm_sizechgs[:,stg],health_event)
            if stg == 3: ax_h[stg,event_num].set_xlabel('Size change (mm^2)')
            #ax_h[stg,event_num].set_ylabel(health_param)
    '''
    # Plot duration of larval stages against healthspans
    fig_h, ax_h = plt.subplots(4,len(health_params))
    fig_h.suptitle('Duration of larval stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_intra_events.T,health_params)):
        for stg in range(4):
            ax_h[stg,event_num].scatter(all_worm_devtime[:,stg],health_event)
            if stg == 3: ax_h[stg,event_num].set_xlabel('Larval stage duration (hr)')
            #ax_h[stg,event_num].set_ylabel(health_param)
    
    bob = input('To interindividual cutoffs')
    
    # Interindividual cutoffs (use some std deviation of parameter across population and times
    # TODO: THINK OF A RIGOROUS WAY TO TREAT TIME TO EVENTS (CUTOFFS NOT JUST MEAN AND STD)
    cutoffs_percentile = 70
    cutoffs_pop = []    # Holds cutoffs over population
    worm_inter_events = np.zeros((len(worm_health),len(health_params)))
    
    for phys_param in health_params:
        #total_mean = np.median(np.array([worm[phys_param][np.where(worm['intensity_90']==-1.0)[0][0] if any(worm['intensity_90']==-1.0) else len(worm['intensity_90'])] \
            #for worm in worm_health]))
        #total_mean = np.median(np.array([worm[phys_param][np.where(worm['intensity_90']==-1.0)[0][0] if any(worm['intensity_90']==-1.0) else len(worm['intensity_90'])] \
            #for worm in worm_health]))
        #total_std = np.std(np.array([worm[phys_param] for worm in worm_health]))
        compiled_pop_data = np.array([worm[phys_param] for worm in worm_health])
        if phys_param != 'intensity_90':
            total_percentile = np.percentile(compiled_pop_data[compiled_pop_data>0],100-cutoffs_percentile)
        else:
            total_percentile = np.percentile(compiled_pop_data[compiled_pop_data>0], cutoffs_percentile)
        cutoffs_pop.append(total_percentile)
        
    #for worm_num,worm in enumerate(worm_health):
        #worm_inter_events[worm_num] = np.array([
            #worm['time'][np.where(worm['intensity_90']<cutoffs_pop[0]*np.ptp(worm['intensity_90']))[0][0]],
            #worm['time'][np.where(
                #(worm['bulk_movement']<cutoffs_pop[1]*np.ptp(worm['bulk_movement'])) \
                #& worm['time']>worm['time'][np.where(worm['bulk_movement']==worm['bulk_movement'].max())[0][0]])],
            #worm['time'][np.where(
                #(worm['stimulated_rate_b']<cutoffs_pop[1]*np.ptp(worm['stimulated_rate_b'])) \
                #& worm['time']>worm['time'][np.where(worm['stimulated_rate_b']==worm['stimulated_rate_b'].max())[0][0]])]
        #])
    for worm_num, worm in enumerate(worm_health):
        last_timepoint = np.where(worm['intensity_90']==-1.0)[0][0] if any(worm['intensity_90']==-1.0) else len(worm['intensity_90'])
        for param_num, (phys_param,cutoff_param) in enumerate(zip(health_params,cutoffs_pop)):
            if phys_param != 'intensity_90':    # Decreasing over time
                try:
                    worm_inter_events[worm_num,param_num] = worm['time'][np.where(
                        (worm[phys_param]<cutoff_param)\
                        & (worm[phys_param]!=-1)
                        & (worm['time']>worm['time'][np.where(worm[phys_param]==worm[phys_param].max())[0][0]]))[0][0]]
                except IndexError as ie:
                    print(worm_num)
                    worm_inter_events[worm_num,param_num]=-1
            else:
                try:
                    worm_inter_events[worm_num,param_num] = worm['time'][np.where(
                        worm[phys_param]>cutoff_param)[0][0]]
                except IndexError as ie:
                    print(worm_num)
                    worm_inter_events[worm_num,param_num] = -1
    worm_inter_events = worm_inter_events * 24 # Convert to hours
    
    # Scatter and calculate correlations
    plt.close('all')
    
    # Plot ending size against healthspans
    fig_h, ax_h = plt.subplots(len(health_params),4)
    fig_h.suptitle('Size at end of stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_inter_events.T,health_params)):
        for stg in range(4):
            ax_h[event_num,stg].scatter(all_worm_sizes[:,stg+1],health_event)
            if event_num==0: ax_h[event_num,stg].set_title('L{}'.format(stg+1))
            if event_num == 3: ax_h[event_num,stg].set_xlabel('Size (mm^2)')
            if stg == 0: ax_h[event_num,stg].set_ylabel(health_param)
            ax_h[event_num,stg].set_xlim([size_xmin[stg],size_xmax[stg]])
    (corr_vals, p_vals) = scipy.stats.spearmanr(
        np.concatenate((all_worm_sizes[:,1:], worm_inter_events), axis = 1))
    print(corr_vals[0:4,4:8])
    print(p_vals[0:4,4:8])
    
    
    # Plot size change per stage against healthspans
    fig_h, ax_h = plt.subplots(len(health_params),4)
    fig_h.suptitle('Size change per stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_inter_events.T,health_params)):
        for stg in range(4):
            ax_h[event_num,stg].scatter(all_worm_sizechgs[:,stg],health_event)
            if event_num==0: ax_h[event_num,stg].set_title('L{}'.format(stg+1))
            if event_num == 3: ax_h[event_num,stg].set_xlabel('Size change (mm^2)')
            if stg == 0: ax_h[event_num,stg].set_ylabel(health_param)
            ax_h[event_num,stg].set_xlim([0,sizechg_xmax[stg]])
    (corr_vals, p_vals) = scipy.stats.spearmanr(
        np.concatenate((all_worm_sizechgs, worm_inter_events), axis = 1))
    print(corr_vals[0:4,4:8])
    print(p_vals[0:4,4:8])
    
    # Plot growth rate per stage against healthspans
    fig_h, ax_h = plt.subplots(len(health_params),4)
    fig_h.suptitle('Growth rate per stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_inter_events.T,health_params)):
        for stg in range(4):
            ax_h[event_num,stg].scatter(all_worm_rates[:,stg],health_event)
            if event_num==0: ax_h[event_num,stg].set_title('L{}'.format(stg+1))
            if event_num == 3: ax_h[event_num,stg].set_xlabel('Growth rate (mm^2/hr)')
            if stg == 0: ax_h[event_num,stg].set_ylabel(health_param)
            ax_h[event_num,stg].set_xlim([0,rate_xmax[stg]])
    (corr_vals, p_vals) = scipy.stats.spearmanr(
        np.concatenate((all_worm_rates, worm_inter_events), axis = 1))
    print(corr_vals[0:4,4:8])
    print(p_vals[0:4,4:8])
    
    # Plot duration of larval stages against healthspans
    fig_h, ax_h = plt.subplots(len(health_params),4)
    fig_h.suptitle('Duration of larval stage')
    for event_num, (health_event,health_param) in enumerate(zip(worm_inter_events.T,health_params)):
        for stg in range(4):
            ax_h[event_num,stg].scatter(all_worm_devtime[:,stg],health_event)
            if event_num==0: ax_h[event_num,stg].set_title('L{}'.format(stg+1))
            if event_num == 3: ax_h[event_num,stg].set_xlabel('Larval stage duration (hr)')
            if stg == 0: ax_h[event_num,stg].set_ylabel(health_param)
    (corr_vals, p_vals) = scipy.stats.spearmanr(
        np.concatenate((all_worm_devtime, worm_inter_events), axis = 1))
    print(corr_vals[0:4,4:8])
    print(p_vals[0:4,4:8])

    bob = input('To trajectory analysis')
    
    
