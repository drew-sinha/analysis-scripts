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

do_reg = False
do_cohorts = False
do_scaling = True
make_labels= False
#out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/age-1_cohorts+regression_20160818/'
#out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/age-1specifichealth_cohorts+regression_20160823/'
out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/spe-9+age-1perage-1SVM_cohorts+reg_20160824/'
out_dir=''
if make_labels and out_dir is not '': out_dir = out_dir+os.path.sep+'labeled'+os.path.sep

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

#def clean_subplot_figure(my_fig, make_labels=False):
    #[clean_plot(my_plot,make_labels) for my_plot in my_fig.get_axes()]
    #if make_labels:
        #for row_num, ax_set in my_fig.get_axes():
            #for col_num, my_plot in ax_set:
                #if row_num is not len(ax_set)-1:
                    #my_plot.set_xlabel('')
                #if row_num is not 0:
                    #my_plot.set_title('')
                #if col_num is not 0:
                    #my_plot.set_ylabel('')


#strains = ['spe-9','age-1']
#strains = ['age-1','age-1specific']
#strains = ['spe-9perage-1SVM','age-1specific']
strains = ['spe-9', 'age-1specific']
strain_dfs = []
for strain in strains:
    with open('/mnt/bulkdata/wzhang/human_dir/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
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
def cohort_traces(my_subfigure, a_variable, adult_df, the_title = None, the_xlabel = None, the_ylabel = None, x_normed = False, y_normed = False, skip_conversion = False, zero_to_one = False, only_worms = None, make_labels=True, bin_width_days=2,bin_mode='day', line_style='-', cohorts_to_use=[], stop_with_death=True):
    '''
    Make cohort traces for a_variable.
    '''
    # Make bins of lifespans.
    (life_cohorts, bin_lifes, my_bins, my_colors) = analyzeHealth.selectData.adult_cohort_bins(adult_df, my_worms = adult_df.worms, bin_width_days = bin_width_days,bin_mode=bin_mode)
    if len(cohorts_to_use) >0:
        life_cohorts = [life_cohorts[c_idx] for c_idx in cohorts_to_use]
        bin_lifes = [bin_lifes[c_idx] for c_idx in cohorts_to_use]
        my_bins = [my_bins[c_idx] for c_idx in cohorts_to_use]
        my_colors = [my_colors[c_idx] for c_idx in cohorts_to_use]
    
    #print(my_bins)
    #print([len(cohort) for cohort in life_cohorts])

    # Exclude worms.    
    if only_worms != None:
        life_cohorts = [[a_worm for a_worm in a_cohort if a_worm in only_worms] for a_cohort in life_cohorts]
    else:
        pass
    my_cohorts = life_cohorts

    # Plot the actual stuff.
    if y_normed:
        if type(a_variable) == type(''):
            mean_start = np.mean(adult_df.mloc(adult_df.worms, [a_variable], ['0.0'])[:, 0, 0])
            (mean_start, my_unit, fancy_name) = adult_df.display_variables(mean_start, a_variable)
        else:
            mean_start = np.mean(a_variable[:, 0, 0])
    for i in range(0, len(my_cohorts)):
        if len(my_cohorts[i]) > 0:
            # Figure out the cohort data.
            a_cohort = my_cohorts[i]
            if type(a_variable) == type(''):
                cohort_data = adult_df.mloc(adult_df.worms, [a_variable])[a_cohort, 0, :]
                variable_name = a_variable
            else:
                cohort_data = a_variable[a_cohort, 0, :]
                variable_name = 'health'
            
            cohort_data = cohort_data[~np.isnan(cohort_data).all(axis = 1)]
            if stop_with_death: # Suppress data after the first individual in a cohort dies
                cohort_data = np.mean(cohort_data, axis = 0)
            else:
                cohort_data = np.nanmean(cohort_data, axis=0)
            if not skip_conversion:
                (cohort_data, my_unit, fancy_name) = adult_df.display_variables(cohort_data, variable_name)
            else:
                (whatever, my_unit, fancy_name) = adult_df.display_variables(cohort_data, variable_name)
            cohort_data = cohort_data[~np.isnan(cohort_data)]
            if y_normed:
                if zero_to_one:
                    cohort_data = cohort_data - cohort_data[0]
                    cohort_data = cohort_data/cohort_data[-1]
                else:
                    cohort_data = cohort_data - cohort_data[-1]
                    cohort_data = cohort_data/cohort_data[0]

            # Figure out the cohort ages.           
            cohort_ages = adult_df.ages[:cohort_data.shape[0]]          
            if x_normed:
                cohort_ages = cohort_ages/np.max(cohort_ages)
            my_subfigure.plot(cohort_ages, cohort_data, color = my_colors[i], linewidth = 2, linestyle=line_style)

    # Label the subplot.
    if the_title == None:
        the_title = fancy_name + ' Over Time'
    if make_labels: my_subfigure.set_title(the_title)
    if make_labels: my_subfigure.set_xlabel(the_xlabel)
    if the_ylabel == None:
        the_ylabel = my_unit
    if make_labels: my_subfigure.set_ylabel(the_ylabel) 
    return my_subfigure

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
    
     #Individual Health bars - 2 columns
    #health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    #health_fig, ax_h = plt.subplots(len(health_vars),len(strains))
    #for strain_health, strain_plots in zip(strain_dfs, ax_h.T):
        #for var_plot, var in zip(strain_plots, health_vars):
            #cohort_traces(var_plot,var, strain_health, make_labels=make_labels, bin_width_days=spe9_bins,bin_mode = 'percentile')
            #clean_plot(var_plot)
    
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
        [clean_plot(var_plot,make_labels) for var_plot in ax_h]
        health_fig.savefig(out_dir+os.path.sep+'health_vars_allcohorts.svg')
        health_fig.savefig(out_dir+os.path.sep+'health_vars_allcohorts.png')
    
    # Overall health
    health_fig, ax_h = plt.subplots(1,1)
    for strain_health, line_style in zip(strain_dfs, ['-','--']):
        cohort_traces(ax_h,'health', strain_health, make_labels=make_labels, bin_width_days=spe9_bins, bin_mode='percentile', line_style = line_style)
    if len(out_dir)>0: 
        clean_plot(var_plot,make_labels)
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
        '''
            max(t1)<max(t2)
        '''
        int_f2 = scipy.interpolate.interp1d(t2,f2)
        
        # Make error function and optimize on it
        my_err_fun = lambda interp_fun_tocompare, time_ref, fun_ref, scale_factor: np.sum((interp_fun_tocompare(scale_factor*time_ref)-fun_ref)**2)
        result = scipy.optimize.minimize_scalar(lambda scale: my_err_fun(int_f2,t1,f1, scale),method='bounded',
            bounds=(0,t2.max()/t1.max()))
        return result
        
    def rescale_data(t_ref,f_ref, t_tofit, f_tofit, normalize_time=True):
        # Pass back equally-spaced data for a reference and new functon associated with its corresponding range
        scale_results = find_scalefactor(t_ref,f_ref,t_tofit,f_tofit)
        if scale_results.status is not 0:
            print('(rescale_data) Warning solver didn\'t figure out what it needed to')
        int_fref = scipy.interpolate.interp1d(t_ref,f_ref)
        int_ffit = scipy.interpolate.interp1d(t_tofit/scale_results.x,f_tofit)   # Scale factor stored in 'x'
        
        domain = np.linspace(t_ref.min(), t_ref.max())
        if not normalize_time:
            return (domain, int_fref(domain), int_ffit(domain))
        else:
            return (domain/domain.max(), int_fref(domain), int_ffit(domain))
        
    # Get cohort assignments from adult_cohort_bins
    animal_bins = np.array([100])
    adult_cohort_bin_data = [analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days=animal_bins,bin_mode='percentile') for strain_df in strain_dfs]
    my_cohort_data_abs = [get_cohort_data(strain_df, cohort_assignments, stop_with_death=False) for strain_df, cohort_assignments in zip(strain_dfs, [bin_data[0] for bin_data in adult_cohort_bin_data])]
    res_scaling_abs = rescale_data(np.array(my_cohort_data_abs[0][1][0]), 
        np.array(my_cohort_data_abs[0][0][0]),
        np.array(my_cohort_data_abs[1][1][0]), 
        np.array(my_cohort_data_abs[1][0][0]))
    print(res_scaling_abs)
    my_cohort_data_zscored = [get_cohort_data(strain_df, cohort_assignments, stop_with_death=False,skip_conversion=True) for strain_df, cohort_assignments in zip(strain_dfs, [bin_data[0] for bin_data in adult_cohort_bin_data])]
    res_scaling_zscored = rescale_data(np.array(my_cohort_data_zscored[0][1][0]), 
        np.array(my_cohort_data_zscored[0][0][0]),
        np.array(my_cohort_data_zscored[1][1][0]), 
        np.array(my_cohort_data_zscored[1][0][0]))
    print(res_scaling_zscored)
    
    # Plot original data
    #fig_h,ax_h=plt.subplots(2,1)
    #ax_h.plot(np.array(my_cohort_data[0][1][0]), 
        #np.array(my_cohort_data[0][0][0]),
        #np.array(my_cohort_data[1][1][0]), 
        #np.array(my_cohort_data[1][0][0]))
    
    # Plot data
    fig_h, ax_h = plt.subplots(2,1)
    plot_data = [(ax_h[0].plot(res_scaling_abs[0], scaled_data))[0] for scaled_data in res_scaling_abs[1:]]
    if make_labels: ax_h[0].legend(plot_data, [strain for strain in strains])
    ax_h[0].set_xlabel('Normalized adult life (rel. to max. lifespan)')
    ax_h[0].set_ylabel('Prognosis (d)')
    plot_data = [(ax_h[1].plot(res_scaling_zscored[0], scaled_data))[0] for scaled_data in res_scaling_zscored[1:]]
    if make_labels: ax_h[1].legend(plot_data, [strain for strain in strains])
    ax_h[1].set_xlabel('Normalized adult life (rel. to max. lifespan)')
    ax_h[1].set_ylabel('Prognosis (z-scored)')
        
    # TODO 
    # Significance testing....
    # Redo with the variants of SVMs. Do combined SVM?
    # Work on rescaling range in addition to time as well?
    if len(out_dir)>0: 
        [clean_plot(my_ax,make_labels) for my_ax in ax_h]
        fig_h.savefig(out_dir+os.path.sep+'health_rescaledtime.svg')
        fig_h.savefig(out_dir+os.path.sep+'health_rescaledtime.png')
    
