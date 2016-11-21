import analyzeHealth
import matplotlib.pyplot as plt
import pickle
import numpy as np
import collections
import os
import numpy as np

plt.ion()
plt.show()
plt.close('all')

do_reg = True
do_cohorts = True
do_scaling = True
make_labels= True
out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/processed_data/age-1_cohorts+regression_20160818/'
#out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/age-1specifichealth_cohorts+regression_20160823/'
#out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/spe-9+age-1perage-1SVM_cohorts+reg_20160824/'
#out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/processed_data/spe-9+age-1percombinedweightedSVM_20161006/'
#out_dir=''
if make_labels and out_dir is not '': 
    out_dir = out_dir+os.path.sep+'labeled'+os.path.sep
    if not os.path.isdir(out_dir): os.mkdir(out_dir)

#health_measures = {
    #'autofluorescence': ['intensity_80'],       
    #'size': ['adjusted_size', 'adjusted_size_rate'],        
    #'eggs': ['cumulative_eggs', 'cumulative_eggs_rate'],        
    #'texture': ['life_texture'],        
    #'movement': ['bulk_movement', 'stimulated_rate_a', 'stimulated_rate_b', 'unstimulated_rate'],
#}

# Unfold health measures all the way
def flatten_list(my_list):
    import collections
    flat_list = []
    
    for my_item in my_list:
        if (not isinstance(my_item, collections.Iterable)) or (type(my_item) is str):
            flat_list.append(my_item)
        else:
            flat_list.extend(flatten_list(my_item))
    return flat_list
health_measure_keys = ['autofluorescence', 'size', 'eggs','texture','movement']
health_measure_values = [
    ['intensity_80'], 
    ['adjusted_size','adjusted_size_rate'], 
    ['cumulative_eggs', 'cumulative_eggs_rate'], 
    ['life_texture'], 
    ['bulk_movement', 'stimulated_rate_a','stimulated_rate_b','unstimulated_rate']]
health_measures = collections.OrderedDict([(k,v) for k,v in zip(health_measure_keys, health_measure_values)])
all_physiology = flatten_list([health_vals for k,health_vals in health_measures.items()])

# 20160926 - Moved to plotting_tools
def clean_plot(my_plot,make_labels=False,suppress_ticklabels=False):
    my_plot.tick_params(axis='both',which='both', top='off', bottom='off', left='off', right='off')
    
    if not suppress_ticklabels:
        full_xticks = my_plot.get_xticks().copy()
        full_xticklabels = my_plot.get_xticklabels().copy()
        my_plot.set_xticks([full_xticks[0], full_xticks[-1]])
        
        full_yticks = my_plot.get_yticks().copy()
        full_yticklabels = my_plot.get_yticklabels().copy()
        my_plot.set_yticks([full_yticks[0], full_yticks[-1]])
    else:
        my_plot.set_xticks([])
        my_plot.set_yticks([])
    
    if not make_labels:
        my_plot.set_xlabel('')
        my_plot.set_ylabel('')
        my_plot.set_title('')

strains = ['spe-9','age-1']
#strains = ['age-1','age-1specific']
#strains = ['spe-9perage-1SVM','age-1specific']
#strains = ['spe-9', 'age-1specific']
#strains = ['spe-9percombinedweightedSVM', 'age-1percombinedweightedSVM']
strain_dfs = []
for strain in strains:
    #with open('/mnt/bulkdata/wzhang/human_dir/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
    with open('/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
        strain_dfs.append(pickle.load(my_file)['adult_df'])
    print(strain+": n = {}".format(len(strain_dfs[-1].worms)))
    print(strain+": CV={:.2f}".format(
        np.std(analyzeHealth.selectData.get_adultspans(strain_dfs[-1]))/np.mean(analyzeHealth.selectData.get_adultspans(strain_dfs[-1]))))

if do_reg:
    biomarker_data = [] #Raw biomarkers
    biomarker_predict_data=[]   # Predicted life left based on biomarkers
    for strain,strain_df in zip(strains,strain_dfs):
        [predicted_health, reg_weights, reg_intercept, actual_health] = analyzeHealth.computeStatistics.multiple_regression_combine(
            strain_df, all_physiology, dependent_variable='ghost_age', return_reg_data=True)    # Mult. reg of health on raw_biomarkers

        r_squared = analyzeHealth.computeStatistics.quick_pearson(
            np.ndarray.flatten(predicted_health), 
            np.ndarray.flatten(actual_health))
        
        # Calculate some r's TODO think about some smart block diagnoal arrangement over all n^2 x-correlations.
        # Build health
        flat_dependent = np.ndarray.flatten(strain_df.mloc(measures = ['ghost_age']))
        all_flats = []
        for a_var in all_physiology:
            flat_var = np.ndarray.flatten(strain_df.mloc(measures = [a_var]))
            all_flats.append(flat_var)
        all_flats = np.array(all_flats).transpose()
        biomarker_data.append([all_flats, flat_dependent])
        
        all_flats_biomarkerpredict = []
        for a_var in health_measure_keys:
            flat_var = np.ndarray.flatten(strain_df.extra_data[a_var])
            all_flats_biomarkerpredict.append(flat_var)
        all_flats_biomarkerpredict = np.array(all_flats_biomarkerpredict).transpose()
        biomarker_predict_data.append([all_flats_biomarkerpredict, flat_dependent])
        
        pearson_rawbiomarker_data = [analyzeHealth.computeStatistics.quick_pearson(my_predictor, flat_dependent) for my_predictor in all_flats.T]
        pearson_predbiomarker_data = [analyzeHealth.computeStatistics.quick_pearson(my_predictor, flat_dependent) for my_predictor in all_flats_biomarkerpredict.T]
        
        # Plot the prognostic vs. actual time to death for animals
        fig_h, ax_h = plt.subplots(1,1)
        ax_h.scatter(predicted_health, actual_health)
        ax_h.set_title(strain+' Pearson mult. reg. r^2={:.3}'.format(r_squared))
        
        fig_h, ax_h = plt.subplots(1,1)
        ax_h.scatter(strain_df.extra_data['health'], actual_health)
        ax_h.set_title(strain+' SVR r^2={:.3}'.format(analyzeHealth.computeStatistics.quick_pearson(    # SVR of actual health on predicted health
            np.ndarray.flatten(strain_df.extra_data['health']), 
            np.ndarray.flatten(actual_health))))
        
        print(strain)
        print(all_physiology)
        print(reg_weights*-1)
        print('Pearson r^2 for raw biomarkers against actual health')
        print(pearson_rawbiomarker_data)
        print('Pearson r^2 for biomarker-based prognosis against actual health')
        print(health_measure_keys)
        print(pearson_predbiomarker_data)
        print('Total r^2s')
        print('Pearson multiple regression raw biomarkers and actual health r^2:{:.3}'.format(r_squared))
        print('SVR mult. reg. b/t predicted health and actual health r^2:{:.3}'.format(analyzeHealth.computeStatistics.quick_pearson(
            np.ndarray.flatten(strain_df.extra_data['health']), 
            np.ndarray.flatten(actual_health))))
    
    # Plot the raw biomarker values against actual lifespan
    rawbiomarker_fig_h, rawbiomarker_ax_h = plt.subplots(len(all_physiology),len(strains))
    for strain_biomarker_data,strain_axs in zip(biomarker_data, rawbiomarker_ax_h.T):
        for (predictor_data, b_ax, b_name) in zip(strain_biomarker_data[0].T,strain_axs,all_physiology):
            b_ax.scatter(predictor_data, strain_biomarker_data[1])
            b_ax.set_ylabel(b_name+'raw',rotation=0)
    if len(out_dir)>0:
        [clean_plot(my_plot, make_labels=make_labels,suppress_ticklabels=True) for my_plot in rawbiomarker_ax_h.flatten()]
        #biomarker_fig_h.savefig(out_dir+'reg_biomarkers.svg')
        rawbiomarker_fig_h.savefig(out_dir+os.path.sep+'reg_biomarkers_raw.png')
    
    predictbiomarker_fig_h, predictbiomarker_ax_h = plt.subplots(len(health_measure_keys),len(strains))
    for strain_biomarker_data,strain_axs in zip(biomarker_predict_data, predictbiomarker_ax_h.T):
        for (predictor_data, b_ax, b_name) in zip(strain_biomarker_data[0].T,strain_axs,health_measure_keys):
            b_ax.scatter(predictor_data, strain_biomarker_data[1])
            b_ax.set_ylabel(b_name+'prognostic',rotation=0)
    if len(out_dir)>0:
        [clean_plot(my_plot, make_labels=make_labels,suppress_ticklabels=True) for my_plot in predictbiomarker_ax_h.flatten()]
        #biomarker_fig_h.savefig(out_dir+'reg_biomarkers.svg')
        predictbiomarker_fig_h.savefig(out_dir+os.path.sep+'reg_biomarkers_prognostic.png')
    
    
    print('done with regression')

##############

def get_prototype_lines(my_plot, feature='linestyle'):
    '''
        Find examples of prototypical lines based on a specified feature
        feature - dict w
    '''
    prototypes = []
    props_observed=[]   
    for my_line in my_plot.get_lines():
        if feature is 'linestyle' and my_line.get_linestyle() not in props_observed:
            prototypes.append(my_line)
            props_observed.append(my_line.get_linestyle())
    return (prototypes, props_observed)
    

if do_cohorts:
    
    # Do global summary
    spe9_bins = np.array([100])
    
    # Inidividual health vars - 1 column
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    health_fig, ax_h = plt.subplots(len(health_vars),1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        for var_plot, var in zip(ax_h.T, health_vars):
            cohort_traces(var_plot,var, strain_health, make_labels=make_labels, bin_width_days=spe9_bins,bin_mode = 'percentile', line_style=line_style, stop_with_death=False)
    if len(out_dir)>0: 
        [clean_plot(var_plot,make_labels) for var_plot in ax_h]
        health_fig.savefig(out_dir+os.path.sep+'health_vars_pop.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_vars_pop.png')
    
    # Overall health - 1 figure
    health_fig, ax_h = plt.subplots(1,1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=spe9_bins, bin_mode='percentile', line_style = line_style, stop_with_death=False)
    if len(out_dir)>0: 
        clean_plot(ax_h,make_labels)
        health_fig.savefig(out_dir+os.path.sep+'health_overall_pop.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_overall_pop.png')
    
    # Do summary for all cohorts (same relative size)
    spe9_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    spe9_bins = 100*np.cumsum(spe9_bins)/sum(spe9_bins) # Make percentile bins
            
    # Individual health vars
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    health_fig, ax_h = plt.subplots(len(health_vars),1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        for var_plot, var in zip(ax_h.T, health_vars):
            cohort_traces(var_plot,var, strain_health, make_labels=make_labels, bin_width_days=spe9_bins,bin_mode = 'percentile', line_style=line_style)
    if len(out_dir)>0: 
        [clean_plot(var_plot,make_labels=make_labels,suppress_ticklabels=not make_labels) for var_plot in ax_h]
        health_fig.savefig(out_dir+os.path.sep+'health_vars_allcohorts.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_vars_allcohorts.png')
    
    # Overall health
    health_fig, ax_h = plt.subplots(1,1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=spe9_bins, bin_mode='percentile', line_style = line_style)
    if len(out_dir)>0: 
        clean_plot(var_plot,make_labels=make_labels,suppress_ticklabels=not make_labels)
        health_fig.savefig(out_dir+os.path.sep+'health_overall_allcohorts.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_overall_allcohorts.png')
        
    # Do summary for 3 cohorts
    spe9_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    spe9_bins = 100*np.cumsum(spe9_bins)/sum(spe9_bins) # Make percentile bins
    cohorts_to_use = [0,int(len(spe9_bins)/2)+1, -1]
    
    # Make some figures with the old analysis
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    health_fig, ax_h = plt.subplots(len(health_vars),1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        for var_plot, var in zip(ax_h.T, health_vars):
            cohort_traces(var_plot,var, strain_health, make_labels=make_labels, bin_width_days=spe9_bins,bin_mode = 'percentile', line_style=line_style, cohorts_to_use=cohorts_to_use)
    if len(out_dir)>0:
        [clean_plot(var_plot) for var_plot in ax_h] 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_3cohorts.svg') 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_3cohorts.png')
    
    health_fig, ax_h = plt.subplots(1,1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=spe9_bins, bin_mode='percentile',line_style=line_style, cohorts_to_use = cohorts_to_use)
    if len(out_dir)>0: 
        clean_plot(ax_h, make_labels)
        health_fig.savefig(out_dir+os.path.sep+'health_overall_3cohorts.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_overall_3cohorts.png')
        
    # Plot global health of mutant vs. spe-9 cohorts
    spe9_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    spe9_bins = 100*np.cumsum(spe9_bins)/sum(spe9_bins) # Make percentile bins
    mt_bins = np.array([100])
    health_fig, ax_h = plt.subplots(1,1)
    cohort_traces(ax_h,'health', strain_dfs[0], make_labels=make_labels, bin_width_days=spe9_bins, bin_mode='percentile', line_style = '-')
    for strain_health, line_style in zip(strain_dfs[1:], ['--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=mt_bins, bin_mode='percentile', line_style = line_style,stop_with_death=False)
    if len(out_dir)>0:
        clean_plot(ax_h)
        health_fig.savefig(out_dir+os.path.sep+'health_vars_mutantvsWT.svg') 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_mutantvsWT.png') 

    
    # Make some figures while not suppressing data from after the first death in cohorts
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    health_fig, ax_h = plt.subplots(len(health_vars),1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        for var_plot, var in zip(ax_h.T, health_vars):
            cohort_traces(var_plot,var, strain_health, make_labels=make_labels, bin_width_days=spe9_bins,bin_mode = 'percentile', line_style=line_style, cohorts_to_use=cohorts_to_use, stop_with_death=False)
    if len(out_dir)>0:
        [clean_plot(var_plot, make_labels) for var_plot in ax_h] 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_3cohorts_unsuppressed.svg') 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_3cohorts_unsuppressed.png')
    
    health_fig, ax_h = plt.subplots(1,1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=spe9_bins, bin_mode='percentile',line_style=line_style, cohorts_to_use = cohorts_to_use, stop_with_death=False)
    if len(out_dir)>0: 
        clean_plot(ax_h,make_labels)
        health_fig.savefig(out_dir+os.path.sep+'health_overall_3cohorts_unsuppressed.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_overall_3cohorts_unsuppressed.png')
        
    # Bin animals in each population into equally sized bins
    #animal_bins = [20,40,60,80,100]
    animal_bins = np.linspace(100/7,100,7)
    
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    health_fig, ax_h = plt.subplots(len(health_vars),1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        for var_plot, var in zip(ax_h.T, health_vars):
            cohort_traces(var_plot,var, strain_health, make_labels=make_labels, bin_width_days=animal_bins,bin_mode = 'percentile', line_style=line_style)
    if len(out_dir)>0:
        [clean_plot(var_plot,make_labels) for var_plot in ax_h] 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_7equalbins.svg') 
        health_fig.savefig(out_dir+os.path.sep+'health_vars_7equalbins.png')
    
    health_fig, ax_h = plt.subplots(1,1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=animal_bins, bin_mode='percentile',line_style=line_style)
    if len(out_dir)>0: 
        #clean_plot(ax_h,make_labels)
        health_fig.savefig(out_dir+os.path.sep+'health_overall_7equalbins.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_overall_7equalbins.png')
    
if do_scaling:
    import scipy.interpolate
    import scipy.optimize
    
    # 20160928 Moved to analyzeHealth.selectData.py
    def get_cohort_data(adult_df, cohort_assignments,variable_to_get='health', stop_with_death=True, skip_conversion=False):
        '''
            cohort_assignments - Output from adult_cohort_bins (lists of indices for worms in adult_df corresponding to each cohort)
            variable_to_get - String corresponding to variable of interest
            normalize_var - Normalize variables into the range [0,1]
        '''
        
        cohort_data = []
        cohort_ages = []
        for a_cohort in cohort_assignments:
            worm_cohort_data = adult_df.mloc(adult_df.worms, [variable_to_get])[a_cohort, 0,:] # Data for individual worms in cohort
            worm_cohort_data = worm_cohort_data[~np.isnan(worm_cohort_data).all(axis=1)]
            if stop_with_death:
                cohort_data.append(np.mean(worm_cohort_data,axis=0))
            else:
                cohort_data.append(np.nanmean(worm_cohort_data,axis=0))
            if not skip_conversion:
                (cohort_data[-1], my_unit, fancy_name) = adult_df.display_variables(cohort_data[-1], variable_to_get)
            cohort_data[-1] = cohort_data[-1][~np.isnan(cohort_data[-1])]
            #if normalize_var: cohort_data[-1] = (cohort_data[-1]-cohort_data[-1][0])/(cohort_data[-1][-1]-cohort_data[-1][0])
            cohort_ages.append(adult_df.ages[:cohort_data[-1].shape[0]])
        return (cohort_data, cohort_ages)
    
    def find_scalefactor(t1,f1,t2,f2):
        def err_fun(t1,f1,t2,f2,time_scale,range_scale):
            int_f1 = scipy.interpolate.interp1d(t1,f1)
            int_f2 = scipy.interpolate.interp1d(time_scale*t2,range_scale*f2)
            
            num_timepts = min(len(t1),len(t2))
            rescaled_t1 = np.linspace(t1.min(),t1.max(),num_timepts)
            rescaled_t2 = time_scale*np.linspace(t2.min(),t2.max(),num_timepts)
            
            return np.sum((int_f2(rescaled_t2)-int_f1(rescaled_t1))**2)
        
        # Make error function and optimize on it
        result = scipy.optimize.minimize(lambda scale_factors: err_fun(t1,f1,t2,f2, scale_factors[0],scale_factors[1]),
            [0.9*t1.max()/t2.max(),f1.max()/f2.max()],
            method='L-BFGS-B',
            bounds=[(0,None), (0,None)])
        return result
        
    def rescale_data(t_ref,f_ref, t_tofit, f_tofit, normalize_fit=True):
        # Pass back equally-spaced data for a reference and new functon associated with its corresponding range
        scale_results = find_scalefactor(t_ref,f_ref,t_tofit,f_tofit)
        print(scale_results.x)
        #domain = np.linspace(t_ref.min(), t_ref.max())
        #if not normalize_fit:
            #return (domain, int_fref(domain), int_ffit(domain))
        #else:
            #return (domain/domain.max(), int_fref(domain)/int_fref(domain).max(), int_ffit(domain)/int_ffit(domain).max())
            
        if normalize_fit:
            num_timepts = min(len(t_ref),len(t_tofit))
            t_ref_resampled = np.linspace(t_ref.min(),t_ref.max(),num_timepts)
            t_tofit_resampled = np.linspace(t_tofit.min(),t_tofit.max(),num_timepts)*scale_results.x[0]
            int_fref = scipy.interpolate.interp1d(t_ref,f_ref/f_ref.max())
            int_ffit = scipy.interpolate.interp1d(t_tofit*scale_results.x[0],f_tofit*scale_results.x[1]/f_ref.max())   # Scale factor stored in 'x'
            
            return(np.linspace(0,1,num_timepts),
                int_fref(t_ref_resampled),
                int_ffit(t_tofit_resampled))
                
        
    # Get cohort assignments from adult_cohort_bins
    animal_bins = np.array([100])
    adult_cohort_bin_data = [analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days=animal_bins,bin_mode='percentile') for strain_df in strain_dfs]
    my_cohort_data_abs = [get_cohort_data(strain_df, cohort_assignments, stop_with_death=False) for strain_df, cohort_assignments in zip(strain_dfs, [bin_data[0] for bin_data in adult_cohort_bin_data])]
    res_scaling_abs = rescale_data(np.array(my_cohort_data_abs[0][1][0]), 
        np.array(my_cohort_data_abs[0][0][0]),
        np.array(my_cohort_data_abs[1][1][0]), 
        np.array(my_cohort_data_abs[1][0][0]))
    print(res_scaling_abs)
    #my_cohort_data_zscored = [get_cohort_data(strain_df, cohort_assignments, stop_with_death=False,skip_conversion=True) for strain_df, cohort_assignments in zip(strain_dfs, [bin_data[0] for bin_data in adult_cohort_bin_data])]
    #res_scaling_zscored = rescale_data(np.array(my_cohort_data_zscored[0][1][0]), 
        #np.array(my_cohort_data_zscored[0][0][0]),
        #np.array(my_cohort_data_zscored[1][1][0]), 
        #np.array(my_cohort_data_zscored[1][0][0]))
    #print(res_scaling_zscored)
    
    # Plot original data
    #fig_h,ax_h=plt.subplots(1,1)
    #ax_h.plot(np.array(my_cohort_data[0][1][0]), 
        #np.array(my_cohort_data[0][0][0]),
        #np.array(my_cohort_data[1][1][0]), 
        #np.array(my_cohort_data[1][0][0]))
    
    # Plot data
    fig_h, ax_h = plt.subplots(1,1)
    plot_data = [(ax_h[0].plot(res_scaling_abs[0], scaled_data))[0] for scaled_data in res_scaling_abs[1:]]
    if make_labels: ax_h[0].legend(plot_data, [strain for strain in strains])
    ax_h[0].set_xlabel('Normalized adult life (rel. to max. lifespan)')
    ax_h[0].set_ylabel('Normalized Prognosis')
        
    # TODO 
    # Significance testing....
    if len(out_dir)>0: 
        [clean_plot(my_ax,make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h]
        fig_h.savefig(out_dir+os.path.sep+'health_rescaledtime.svg')
        fig_h.savefig(out_dir+os.path.sep+'health_rescaledtime.png')
    
