import pandas as pd
import plotting_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

def plot_spanseries(spans,ax_h = None,**kws):
    '''
        Helper for plotting a single set of spans (e.g. lifespans) on an axis)
    '''
    if ax_h is None:
        fig_h, ax_h = plt.subplots(1,1)
        ax_provided = False
    else:
        ax_provided = True
    sorted_s = np.sort(spans[~np.isnan(spans)]) #handling the nan eliminates potential leak-through of alive animals when used with annotation code
    prop_alive = (1 - (np.arange(start=0,stop=np.size(spans[~np.isnan(spans)])))/np.size(spans))
    
    # Generate all the points needed for plotting
    # Plot the first point (0,1), then a straight line across and down for each event/death
    x_vals = np.sort(np.concatenate(([0],sorted_s,sorted_s)))
    y_vals = np.sort(np.concatenate(([1,1],prop_alive[:-1],prop_alive[:-1],[prop_alive[-1]])))[::-1] # Reverse this.
    
    data = ax_h.plot(x_vals,y_vals,**kws)
    
    if not ax_provided:
        return [fig_h,ax_h,data]
    else:
        return data

def plot_manual_ls(annotation_fns):
    '''
        Plotting tool for plotting survival curves from manually curated tsv files (e.g. identified individuals in worm corrals)
    '''
    compiled_ann_data = [pd.read_csv(ann_file,sep='\t') for ann_file in annotation_fns]
    fig_h, ax_h = plt.subplots(1,1)
    data_series = []
    metadata = []
    compiled_ls = []
    for ann_data, ann_color in zip(compiled_ann_data,plotting_tools.qual_colors):
        good_worms = pd.Series(['DEAD' in note for note in ann_data['Notes']])
        good_ls = ann_data['Death'][good_worms]-ann_data['Hatch'][good_worms]
        data_series.append(plot_survival(good_ls, ax_h=ax_h,color=ann_color))
        metadata.append({'n':len(good_ls)})
        compiled_ls.append(good_ls)
    ax_h.set_ylim([0,1])
    ax_h.set_xlabel('Time Post-Transfer (d)')
    ax_h.set_ylabel('Proportion Surviving')
    
    print('Means:')
    [print(np.nanmean(ls)) for ls in compiled_ls] # Excludes alive/non-annotated worms!
    
    print('Medians:')
    [print(np.nanmedian(ls)) for ls in compiled_ls] # Excludes alive/non-annotated worms!
    
    print('T-test for lifespan of *dead* animals')
    clean_compiled_ls = [ls[~np.isnan(ls)] for ls in compiled_ls] #Take care of worms that are still alive
    print(scipy.stats.ttest_ind(*clean_compiled_ls,equal_var=False))

    return [fig_h, ax_h, data_series, metadata]

#~ def plot_lifespan(ann_fps, expt_mds, annotation_prefix_list = [], bad_worm_kws=[],debug_mode=False,plot_mode='both',hist_mode='kde'):
    #~ def draw_hist(data, my_ax, hist_mode = 'kde'):
        #~ if hist_mode is 'kde':
            #~ interval_supp,interval_data,interval_obj = zplib.scalar_stats.kde.kd_distribution(data)
            #~ my_ax.plot(interval_supp,interval_data)
        #~ elif hist_mode is 'bar':
            #~ my_ax.hist(data)
        #~ return my_ax
        
    #~ if annotation_prefix_list == []:
        #~ my_ann_files = [annotation_file.AnnotationFile(ann_fp) for ann_fp in ann_fps]
    #~ else:
        #~ my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix=prefix) for (ann_fp,prefix) in zip(ann_fps,annotation_prefix_list)]
    
    #~ timestamped_data = {}
    #~ [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in list(my_ann_files[0].data.keys())]
    #~ for [expt_md_fp, ann_file] in zip(expt_mds, my_ann_files):
        #~ if debug_mode: print(expt_md_fp)
        #~ if type(expt_mds[0]) == str:    # Simple - one md file
            #~ ann_file_data = ann_file.data_as_timestamps_simple(expt_md_fp)
        #~ if type(expt_mds[0]) == dict:   # Need to link multiple md files
            #~ ann_file_data = ann_file.data_as_timestamps(expt_md_fp)
            
        #~ for expt_key in timestamped_data.keys():
            #~ timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
    
    #~ viable_worm = (timestamped_data['Hatch']!=-1) \
        #~ & (timestamped_data['Death']!=-1) \
        #~ & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    #~ print(timestamped_data['Worm'][viable_worm])
    
    #~ lifespan = (timestamped_data['Death']-timestamped_data['Hatch'])[viable_worm]/(3600*24) # Days
    #~ sorted_ls = sorted(lifespan)
    #~ prop_alive = 1 - (np.arange(start=0,stop=np.size(lifespan)))/np.size(lifespan)
    
    #~ plt.show()
    #~ plt.ion()
    
    #~ plt.gcf().clf()
    #~ if plot_mode is 'both':
        #~ fig_h, ax_h = plt.subplots(2,1,sharex=True)
        #~ ax_h[0].plot(np.append([0],sorted_ls),np.append([1],prop_alive))
        #~ ax_h[0].set_xlabel('Time since expt. start (d)')
        #~ ax_h[0].set_ylabel('Proportion alive')
        #~ ax_h[0].set_title('Survival curve - n = {}'.format(np.size(lifespan)))
        
        #~ ax_h[1].hist(lifespan)
        #~ ax_h[1].set_xlabel('Time to death (d)')
        #~ ax_h[1].set_ylabel('Frequency')
        #~ ax_h[1].set_title('Mean+/-STD: {:.2f}+/-{:.2f}d\nMedian:{}d'.format(np.mean(lifespan),np.std(lifespan),np.median(lifespan)))
    #~ elif plot_mode is 'lifespan':
        #~ fig_h,ax_h = plt.subplots(1,1)
        #~ draw_hist(lifespan, ax_h)
        #~ ax_h.set_xlabel('Time to death (d)')
        #~ ax_h.set_ylabel('Frequency')
        #~ ax_h.set_title('Mean+/-STD: {:.2f}+/-{:.2f}d\nMedian:{}d'.format(np.mean(lifespan),np.std(lifespan),np.median(lifespan)))
    #~ else:
        #~ fig_h,ax_h = plt.subplots(1,1)
        #~ ax_h.plot(np.append([0],sorted_ls),np.append([1],prop_alive))
        #~ ax_h.set_xlabel('Time since expt. start (d)')
        #~ ax_h.set_ylabel('Proportion alive')
        #~ ax_h.set_title('Survival curve - n = {}'.format(np.size(lifespan)))
    
    #~ return fig_h, ax_h
