import pandas as pd
import plotting_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import pathlib

import corral_annotations.annotation_file as annotation_file

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

def plot_manual_ls(annotation_fns,ax_h=None):
    '''
        Plotting tool for plotting survival curves from manually curated tsv files (i.e. not automated scope experiments)
        
        interval - Conversion from timepoint number (i.e. in 
    '''
    compiled_ann_data = [pd.read_csv(ann_file,sep='\t') for ann_file in annotation_fns]
    data_series = []
    metadata = []
    compiled_ls = []
    
    if ax_h is None:
        fig_h, ax_h = plt.subplots(1,1)
        ax_provided = False
    else: ax_provided = True
    
    for ann_data, ann_color in zip(compiled_ann_data,plotting_tools.qual_colors):
        good_worms = pd.Series(['DEAD' in note for note in ann_data['Notes']])
        good_ls = ann_data['Death'][good_worms]-ann_data['Hatch'][good_worms]
        data_series.append(plot_spanseries(good_ls, ax_h=ax_h,color=ann_color))
        metadata.append({'n':len(good_ls)})
        compiled_ls.append(good_ls)
    ax_h.set_ylim([0,1])
    ax_h.set_xlabel('Days Post-Transfer')
    ax_h.set_ylabel('Proportion Surviving')
    
    print('Means:')
    [print(np.nanmean(ls)) for ls in compiled_ls] # Excludes alive/non-annotated worms!
    
    print('Medians:')
    [print(np.nanmedian(ls)) for ls in compiled_ls] # Excludes alive/non-annotated worms!
    
    print('T-test for lifespan of *dead* animals')
    clean_compiled_ls = [ls[~np.isnan(ls)] for ls in compiled_ls] #Take care of worms that are still alive
    print(scipy.stats.ttest_ind(*clean_compiled_ls,equal_var=False))
    
    if not ax_provided:
        return (fig_h, ax_h, data_series, metadata)
    else:
        return (ax_h, data_series, metadata)
        
def plot_expt_ls(expt_dirs, ax_h=None, calc_adultspan=False,bad_worm_kws=[], **plot_kws):
    '''
        Plot lifespans corresponding to one or more single longitudinal scope experiments
    '''
    
    if type(expt_dirs) is not list: expt_dirs = [expt_dirs] # Single path
    if type(expt_dirs[0]) is str: expt_dirs = list(map(pathlib.Path,expt_dirs))
    
    ax_provided = ax_h is not None
    if not ax_provided: fig_h, ax_h = plt.subplots()

    expt_afs = [annotation_file.AnnotationFile(
        [my_file for my_file in my_dir.iterdir() if my_file.is_file() 
        and '.tsv' == my_file.suffix][0]) 
    for my_dir in expt_dirs]
    
    ts_data = [(my_af.data_as_timestamps_simple(
        my_dir / 'experiment_metadata.json',restricted_list=my_af.get_goodworms(bad_worm_kws=bad_worm_kws))) 
        for my_af,my_dir in zip(expt_afs, expt_dirs)]

    if calc_adultspan: initial_time_label = 'First Egg Laid'
    else: initial_time_label = 'Hatch'
    
    combined_data_series = [plot_spanseries(
            (my_data['Death']-my_data[initial_time_label])/(3600*24),
            ax_h=ax_h)[0] 
        for my_data in ts_data]

    if initial_time_label == 'Hatch':
        ax_h.set_xlabel('Days Post-Hatch')
    elif initial_time_label == 'First Egg Laid':
        ax_h.set_xlabel('Days Post-Adulthood')
    ax_h.set_ylabel('Proportion Surviving')

    # Print stats
    print('Avg Lifespan')
    [print('{}: {:.3f}'.format(my_path.parts[-1],((my_data['Death']-my_data[initial_time_label])/(3600*24)).mean())) for my_path,my_data in zip(expt_dirs,ts_data)]
    print('All: {:.3f}'.format(np.array([
        ((my_data['Death']-my_data[initial_time_label])/(3600*24)).mean() for my_data in ts_data
    ]).mean()))
    
    print ('Median Lifespan')
    [print('{}: {:.3f}'.format(my_path.parts[-1],((my_data['Death']-my_data[initial_time_label])/(3600*24)).median())) for my_path,my_data in zip(expt_dirs,ts_data)]
    
    if ax_provided:
        return ax_h, combined_data_series
    else:
        return fig_h, ax_h, combined_data_series
