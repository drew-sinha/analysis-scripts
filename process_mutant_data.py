import matplotlib.pyplot as plt
import matplotlib.lines
import pickle
import numpy as np
import collections
import os
import numpy as np
import scipy.interpolate
import scipy.optimize
import scipy.stats as stats

import analyzeHealth
import graphingFigures
import plotting_tools
import survival_plotting

import zplib.scalar_stats.smoothing as smoothing
import zplib.scalar_stats.kde

plt.ion()
plt.show()
plt.close('all')

def load_strain_data(strains,health_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/'):
    strain_data = []
    for strain in strains:
        #with open('/mnt/bulkdata/wzhang/human_dir/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
        with open(health_dir+os.path.sep+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
            strain_data.append(pickle.load(my_file))
        print(strain+": n = {}".format(len(strain_data[-1]['adult_df'].worms)))
        print(strain+": CV={:.2f}".format(
            np.std(analyzeHealth.selectData.get_adultspans(strain_data[-1]['adult_df']))/np.mean(analyzeHealth.selectData.get_adultspans(strain_data[-1]['adult_df']))))
    return strain_data

def plot_survival(strain_dfs,strains, out_dir='', make_labels=True, plot_mode = 'pop',ax_h=None):
    # Plot survival curves and lifespan distributions for each strain
    
    if plot_mode is 'cohorts':
        # Use percentile bins
        animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
        animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
        
        survival_fig, ax_h = plt.subplots(len(strains),2)
        for strain_health, strain_ax in zip(strain_dfs,ax_h.T):
            graphingFigures.cannedFigures.survival_lifespan(strain_ax[0],strain_ax[1],strain_health, make_labels=make_labels,
                cohort_info=analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = animal_bins,bin_mode='percentile'))
        [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]

    elif plot_mode in ['pop','population']:
        
        if ax_h is None:
            survival_fig, ax_h = plt.subplots(1,1)
            ax_provided = False
        else: ax_provided = True
        max_life = max([max(analyzeHealth.selectData.get_adultspans(strain_df))/24 for strain_df in strain_dfs])//1 + 1
        data_series = []
        
        for strain,strain_df,strain_color in zip(strains,strain_dfs,plotting_tools.qual_colors[:len(strains)]):
            lifespans = analyzeHealth.selectData.get_lifespans(strain_df)/24
            data_series.append(survival_plotting.plot_spanseries(lifespans,ax_h=ax_h,color=strain_color,linewidth=2.0))
            print('Stats for strain '+strain+' +: ({:.2f}+/-{:.2f})'.format(np.mean(lifespans),np.std(lifespans)))
            print('Comparison for average lifespan (t-test): t = {:.4f}, p = {}'.format(*scipy.stats.ttest_ind(analyzeHealth.selectData.get_lifespans(strain_dfs[0])/24,lifespans,equal_var=False)))
        
        ax_h.set_ylabel('Percent Survival')
        ax_h.set_xlabel('Days Post-Maturity')
        ax_h.set_ylim([0,1.1])
        ax_h.set_xlim([0,max_life])
        ax_h.legend(plotting_tools.flatten_list(data_series),
            [('spe-9' if strain=='spe-9' else 'spe-9;'+strain) + ' (n={})'.format(len(strain_df.worms)) for strain,strain_df in zip(strains,strain_dfs)],
            frameon=False, bbox_to_anchor=(1.1,1.1))
        plotting_tools.clean_plot(ax_h, make_labels=make_labels,suppress_ticklabels=not make_labels)
    
    if len(out_dir)>0:
        survival_fig.savefig(out_dir+os.path.sep+'survival_curves.png')
        survival_fig.savefig(out_dir+os.path.sep+'survival_curves.svg')
    
    if plot_mode is 'cohorts':
        if not ax_provided: return (survival_fig, ax_h)
        else: return ax_h
    elif plot_mode in ['pop','population']:
        if not ax_provided: return (survival_fig, ax_h, data_series)
        else: return [ax_h, data_series]

def unique_items(seq):  # Credit to Peterbe
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

def test_lifespan_replicates(strain_df):
    rep_labels = unique_items([(' '.join(worm_label.split(' ')[:-1])) for worm_label in strain_df.worms])
    rep_labels = unique_items([rep_label if rep_label[-1].isnumeric() else rep_label[:-1] for rep_label in rep_labels])
    print(rep_labels)
    worm_assignments = np.array(
        [num for worm_label in strain_df.worms for num,rep_label in enumerate(rep_labels) if rep_label in ' '.join(worm_label.split(' ')[:-1])])
    lifespans = analyzeHealth.selectData.get_lifespans(strain_df)/24
    
    survival_fig, ax_h = plt.subplots(1,1)
    max_life = max(lifespans)//1 + 1
    data_series = []
    stats = []
    
    for num,rep_label in enumerate(rep_labels):
        sorted_ls = sorted(lifespans[worm_assignments==num])
        prop_alive = 1 - np.arange(start=0,stop=np.size(sorted_ls))/np.size(sorted_ls)
        data_series.append(
            ax_h.plot(np.append([0],sorted_ls),np.append([1],prop_alive),linewidth=2.0))
        stats.append([np.mean(sorted_ls),np.std(sorted_ls)])
        print('Stats for rep '+rep_label+': ({:.2f}+/-{:.2f})'.format(stats[-1][0],stats[-1][1]))
    
    ax_h.set_ylabel('Percent Survival')
    ax_h.set_xlabel('Days Post-Maturity')
    ax_h.set_ylim([0,1.1])
    ax_h.set_xlim([0,max_life])
    ax_h.legend(plotting_tools.flatten_list(data_series),
        [rep_label + ' (n={})'.format(np.count_nonzero(worm_assignments==num)) for num,rep_label in enumerate(rep_labels)],
        frameon=False)
    print('Intertrial variability across replicates for tested strain: {:.2f} +/- {:.2f}'.format(
        np.mean([item[0] for item in stats]),np.std([item[0] for item in stats])))
    
    return (survival_fig,ax_h)
        

def test_adultspans(strain_dfs,strain_labels):
    figs, axs = [], []
    
    for strain_df,strain_label in zip(strain_dfs,strain_labels):
        adultspans = analyzeHealth.selectData.get_adultspans(strain_df)/24
        lifespans = analyzeHealth.selectData.get_lifespans(strain_df)/24
        devspans = lifespans-adultspans
        
        fig_h,ax_h = plt.subplots(1,1)
        ax_h.scatter(adultspans,lifespans)
        ax_h.set_title(strain_label+'\nr={:.2f}\np={:.2f}'.format(
            *stats.pearsonr(adultspans,lifespans)))
        
        figs.append(fig_h)
        axs.append(ax_h)
        #~ fig_h,ax_h = plt.subplots(1,1)
        #~ ax_h.scatter(devspans,lifespans)
        #~ ax_h.set_title('r={:.2f}\np={:.2f}'.format(
            #~ *stats.pearsonr(devspans,lifespans)))
    return(figs,axs)

def generate_regression_data(strain_dfs,strains=None,out_dir='',make_labels=True,verbose=True):
    if strains is None:
        strains = [{}.format(i) for i,strain in enumerate(strains)]
    
    health_measure_keys = ['autofluorescence', 'size', 'eggs','texture','movement']
    health_measure_values = [
        ['intensity_80'], 
        ['adjusted_size','adjusted_size_rate'], 
        ['cumulative_eggs', 'cumulative_eggs_rate'], 
        ['life_texture'], 
        ['bulk_movement', 'stimulated_rate_a','stimulated_rate_b','unstimulated_rate']]
    health_measures = collections.OrderedDict([(k,v) for k,v in zip(health_measure_keys, health_measure_values)])
    all_physiology = plotting_tools.flatten_list([health_vals for k,health_vals in health_measures.items()])
    
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
        rawbiomarker_fig_h.savefig(out_dir+os.path.sep+'reg_biomarkers_raw.png')
    
    predictbiomarker_fig_h, predictbiomarker_ax_h = plt.subplots(len(health_measure_keys),len(strains))
    for strain_biomarker_data,strain_axs in zip(biomarker_predict_data, predictbiomarker_ax_h.T):
        for (predictor_data, b_ax, b_name) in zip(strain_biomarker_data[0].T,strain_axs,health_measure_keys):
            b_ax.scatter(predictor_data, strain_biomarker_data[1])
            b_ax.set_ylabel(b_name+'prognostic',rotation=0)
    if len(out_dir)>0:
        [clean_plot(my_plot, make_labels=make_labels,suppress_ticklabels=True) for my_plot in predictbiomarker_ax_h.flatten()]
        predictbiomarker_fig_h.savefig(out_dir+os.path.sep+'reg_biomarkers_prognostic.png')
        

def add_strain_legend(ax_h,strains,data_series=None):
    if data_series is None:
        last_idx = [idx for idx,child in enumerate(ax_h.get_children()) if type(child) is not matplotlib.lines.Line2D][0]
        data_series = ax_h.get_children()[last_idx-len(strains):last_idx]
        ax_h.legend(plotting_tools.flatten_list(data_series),
            #[('spe-9;'+strain if strain != 'spe-9' else strain) for strain in strains],
            strains,
            frameon=False, bbox_to_anchor=(1.1,1.1))

def plot_strain_health(strain_dfs,strains=None,group_mode='population',make_labels=True,out_dir='',collapse_mutant=False,var_plot_mode='combined'):
    default_fnames = {
        'population':'pop',
        'spe9_allcohorts':'allcohorts',
        'spe9_3cohorts':'3cohorts',
        'equalbins_7':'7equalbins',
        'equalbins_5':'5equalbins',
    }
    custom_fname = default_fnames[group_mode]
    
    if group_mode is 'population':
        animal_bins = np.array([100])
        cohorts_to_use=[]
    elif group_mode is 'spe9_allcohorts':
        animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
        animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
        cohorts_to_use = []
    elif group_mode is 'spe9_3cohorts':
        animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
        animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
        cohorts_to_use = [0,int(len(animal_bins)/2)+1, -1]
    elif group_mode is 'equalbins_7':
        animal_bins = np.linspace(100/7,100,7)
        cohorts_to_use = []
    elif group_mode is 'equalbins_5':
        animal_bins = np.linspace(100/5,100,5)
        cohorts_to_use = []
    
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    if var_plot_mode is 'combined':
        var_fig, var_ax = plt.subplots(len(health_vars),1)
        for var_plot,var in zip(var_ax.T,health_vars):
            graphingFigures.cannedFigures.cohort_traces(var_plot,var,strain_dfs[0],make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style='-',cohorts_to_use=cohorts_to_use)
            for strain_health,line_style in zip(strain_dfs[1:],['--','-.',':']):
                if collapse_mutant:
                    graphingFigures.cannedFigures.cohort_traces(var_plot,var,strain_health,make_labels=make_labels,bin_width_days=np.array([100]),bin_mode='percentile',line_style=line_style,cohorts_to_use=cohorts_to_use,stop_with_death=False)
                else:
                    graphingFigures.cannedFigures.cohort_traces(var_plot,var,strain_health,make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style=line_style,cohorts_to_use=cohorts_to_use)
        if len(out_dir)>0:
            [plotting_tools.clean_plot(var_plot,make_labels,suppress_ticklabels=not make_labels) for var_plot in var_ax] 
            if collapse_mutant is not True: 
                var_fig.savefig(out_dir+os.path.sep+'health_vars_'+custom_fname+'.svg')
                var_fig.savefig(out_dir+os.path.sep+'health_vars_'+custom_fname+'.png')
            else:
                var_fig.savefig(out_dir+os.path.sep+'health_vars_'+custom_fname+'_mutantvsWT.svg')
                var_fig.savefig(out_dir+os.path.sep+'health_vars_'+custom_fname+'_mutantvsWT.png')
    elif var_plot_mode is 'separate':
        for var in health_vars:
            var_fig, var_ax = plt.subplots(1,1)
            graphingFigures.cannedFigures.cohort_traces(var_ax,var,strain_dfs[0],make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style='-',cohorts_to_use=cohorts_to_use)
            for strain_health,line_style in zip(strain_dfs[1:],['--','-.',':']):
                if collapse_mutant:
                    graphingFigures.cannedFigures.cohort_traces(var_ax,var,strain_health,make_labels=make_labels,bin_width_days=np.array([100]),bin_mode='percentile',line_style=line_style,cohorts_to_use=cohorts_to_use,stop_with_death=False)
                else:
                    graphingFigures.cannedFigures.cohort_traces(var_ax,var,strain_health,make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style=line_style,cohorts_to_use=cohorts_to_use)
            if len(out_dir)>0:
                plotting_tools.clean_plot(var_ax,make_labels,suppress_ticklabels=not make_labels)
                if collapse_mutant is not True: 
                    var_fig.savefig(out_dir+os.path.sep+'health_var_'+var+'_'+custom_fname+'.svg')
                    var_fig.savefig(out_dir+os.path.sep+'health_var_'+var+'_'+custom_fname+'.png')
                else:
                    var_fig.savefig(out_dir+os.path.sep+'health_var_'+var+'_'+custom_fname+'_mutantvsWT.svg')
                    var_fig.savefig(out_dir+os.path.sep+'health_var_'+var+'_'+custom_fname+'_mutantvsWT.png')
    
    health_fig, health_ax = plt.subplots(1,1)
    graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_dfs[0],make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style='-')
    #~ for strain_health,line_style in zip(strain_dfs[1:],['--','-.',':']):
        #~ if collapse_mutant:
            #~ graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_health,make_labels=make_labels,bin_width_days=np.array([100]),bin_mode='percentile',line_style=line_style,stop_with_death=False)
        #~ else:
            #~ graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_health,make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style=line_style,cohorts_to_use=cohorts_to_use)
    if collapse_mutant: # Use qualitative color map to compare populations
        for strain_health,line_color in zip(strain_dfs[1:],plotting_tools.qual_colors[1:len(strain_dfs)]):
            graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_health,make_labels=make_labels,bin_width_days=np.array([100]),bin_mode='percentile',line_color=line_color,stop_with_death=False)
        add_strain_legend(health_ax,strains[1:])
    else:
        for strain_health,line_style in zip(strain_dfs[1:],['--','-.',':']):
            graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_health,make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style=line_style,cohorts_to_use=cohorts_to_use)

    if len(out_dir)>0: 
        plotting_tools.clean_plot(health_ax,make_labels,suppress_ticklabels=not make_labels)
        #~ plotting_tools.clean_plot(health_ax,make_labels,suppress_ticklabels=False)
        if collapse_mutant is not True: 
            health_fig.savefig(out_dir+os.path.sep+'health_overall_'+custom_fname+'.svg')
            health_fig.savefig(out_dir+os.path.sep+'health_overall_'+custom_fname+'.png')
        else:
            health_fig.savefig(out_dir+os.path.sep+'health_overall_'+custom_fname+'_mutantvsWT.svg')
            health_fig.savefig(out_dir+os.path.sep+'health_overall_'+custom_fname+'_mutantvsWT.png')
    
    return [[health_fig,health_ax],[var_fig,var_ax]]
   
def normalize_curve(f_tonorm, direction='decreasing'):
    if direction is 'decreasing':
        min_val = np.mean(f_tonorm[int(len(f_tonorm)*0.95):])
        max_val = np.mean(f_tonorm[:np.ceil(len(f_tonorm)*0.05)])
    else:
        max_val = np.mean(f_tonorm[int(len(f_tonorm)*0.95):])
        min_val = np.mean(f_tonorm[:np.ceil(len(f_tonorm)*0.05)])

    return (f_tonorm-min_val)/(max_val-min_val)

def find_scalefactor(t1,f1,t2,f2):
    def err_fun(t1,f1,t2,f2,time_scale,range_scale,offset):
        int_f1 = scipy.interpolate.interp1d(t1,f1)
        int_f2 = scipy.interpolate.interp1d(time_scale*t2,range_scale*f2+offset)
        
        num_timepts = min(len(t1),len(t2))
        rescaled_t1 = np.linspace(t1.min(),t1.max(),num_timepts)
        rescaled_t2 = time_scale*np.linspace(t2.min(),t2.max(),num_timepts)
        
        return np.sum((int_f2(rescaled_t2)-int_f1(rescaled_t1))**2)
    
    # Make error function and optimize on it
    result = scipy.optimize.minimize(lambda scale_factors: err_fun(t1,f1,t2,f2, scale_factors[0],scale_factors[1],scale_factors[2]),
        [0.9*t1.max()/t2.max(),f1.max()/f2.max(),-1*f2.min()],
        method='L-BFGS-B',
        bounds=[(0,None), (0,None),(None,0)])
    return result

def derive_bestscaling(t_ref,f_ref, t_tofit, f_tofit,num_timepts=None):
    # Pass back equally-spaced data for a reference and new functon associated with its corresponding range
    scale_results = find_scalefactor(t_ref,f_ref,t_tofit,f_tofit)
    #print(scale_results.x)
        
    #num_timepts = min(len(t_ref),len(t_tofit))
    num_timepts = len(t_ref) if num_timepts == None else num_timepts
    t_ref_resampled = np.linspace(t_ref.min(),t_ref.max(),num_timepts)
    t_tofit_resampled = np.linspace(t_tofit.min(),t_tofit.max(),num_timepts)*scale_results.x[0]
    #int_fref = scipy.interpolate.interp1d(t_ref,(f_ref-f_ref.min())/(f_ref.max()-f_ref.min()))
    int_fref = scipy.interpolate.interp1d(t_ref,normalize_curve(f_ref))
    #int_ffit = scipy.interpolate.interp1d(t_tofit*scale_results.x[0],
        #(f_tofit*scale_results.x[1]+scale_results.x[2])/((f_tofit.max()*scale_results.x[1]+scale_results.x[2])-(f_tofit.min()*scale_results.x[1]+scale_results.x[2])))   # Scale factor stored in 'x'
    int_ffit = scipy.interpolate.interp1d(t_tofit*scale_results.x[0],
        normalize_curve(f_tofit*scale_results.x[1]+scale_results.x[2]))
    
    
    return(np.linspace(0,1,num_timepts),
        int_fref(t_ref_resampled),
        int_ffit(t_tofit_resampled))
        
def plot_strain_rescaling(strain_dfs,strains,out_dir='', make_labels=True, do_bootstrap = True, percent_range = [0,100]):
    '''
        percent_range - range of lifespans to specify which worms to grab from the df
    '''
    animal_bins = np.array([percent_range[1]])
    if percent_range[1] != 100:
        animal_bins = np.append(animal_bins,100)
    if percent_range[0] != 0:
        animal_bins = np.append([percent_range[0]],animal_bins)
    bin_id = np.where(np.array(animal_bins)==percent_range[1])[0][0]
    #print(bin_id)
    # This should be [#,#?,100]
    if (animal_bins==np.array([100])).all(): stop_with_death=False
    else: stop_with_death=True
            
    # Experimental fitting
    adult_cohort_bin_data = [analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days=animal_bins,bin_mode='percentile') for strain_df in strain_dfs]
    compiled_cohort_data = [analyzeHealth.selectData.get_cohort_data(strain_df, cohort_assignments, stop_with_death=stop_with_death) for strain_df, cohort_assignments in zip(strain_dfs, [bin_data[0] for bin_data in adult_cohort_bin_data])] #[strain, ..., bin]
    res_scaling = [derive_bestscaling(np.array(compiled_cohort_data[0][1][bin_id]), 
            np.array(compiled_cohort_data[0][0][bin_id]),
            np.array(strain_cohort_data[1][bin_id]), 
            np.array(strain_cohort_data[0][bin_id])) for strain_cohort_data in compiled_cohort_data[1:]]
    
    if do_bootstrap:
        num_iterations = 1000
        
        resampled_traj = []
        worms_to_resample = [[np.array(strain_df.worms)[np.random.randint(low=0,high=len(strain_df.worms),size=len(strain_df.worms))] for i in range(num_iterations)] for strain_df in strain_dfs]
        WT_rep_bins = [analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = worm_set, bin_width_days=animal_bins,bin_mode='percentile') for worm_set in worms_to_resample[0]]
        WT_rep_data = [analyzeHealth.selectData.get_cohort_data(strain_dfs[0], cohort_assignments=replicate_bin_assignments,my_worms=worm_set, stop_with_death=stop_with_death) for worm_set, replicate_bin_assignments in zip(worms_to_resample[0], [bin_data[0] for bin_data in WT_rep_bins])]
        resampled_traj.append([scipy.interpolate.interp1d(rep_data[1][bin_id],normalize_curve(rep_data[0][bin_id]))(np.linspace(min(rep_data[1][bin_id]),max(rep_data[1][bin_id]),len(compiled_cohort_data[0][1][bin_id]))) 
            for rep_data in WT_rep_data])
        
        for strain_df, strain_worm_replicates in zip(strain_dfs[1:],worms_to_resample[1:]):
            compiled_rep_bins = [analyzeHealth.selectData.adult_cohort_bins(strain_df, my_worms = worm_set, bin_width_days=animal_bins,bin_mode='percentile') for worm_set in strain_worm_replicates]
            compiled_rep_data = [analyzeHealth.selectData.get_cohort_data(strain_df, cohort_assignments=replicate_bin_assignments,my_worms=worm_set, stop_with_death=stop_with_death) for worm_set, replicate_bin_assignments in zip(strain_worm_replicates, [bin_data[0] for bin_data in compiled_rep_bins])]
            
            compiled_scaling = [derive_bestscaling(np.array(WT_rep[1][bin_id]),np.array(WT_rep[0][bin_id]),
                np.array(mutant_rep[1][bin_id]),np.array(mutant_rep[0][bin_id]),num_timepts=len(compiled_cohort_data[0][1][bin_id])) for WT_rep,mutant_rep in zip(WT_rep_data,compiled_rep_data)]
            resampled_traj.append(np.array([scaling[2] for scaling in compiled_scaling]))

    fig_data = []
    for strain_num,strain_scaling in enumerate(res_scaling):
        fig_h, ax_h = plt.subplots(1,1)
        plot_data = [(ax_h.plot(strain_scaling[0], scaled_data,linewidth=2))[0] for scaled_data in strain_scaling[1:]]
        if do_bootstrap:
            [ax_h.fill_between(np.linspace(0,1,len(compiled_cohort_data[0][1][bin_id])),*np.percentile(strain_curves,[2.5,97.5],axis=0),color=(strain_color+1)/2) \
                for (strain_curves, strain_color) in zip([resampled_traj[0],resampled_traj[strain_num+1]],np.array([[0,0,1],[0,1,0]]))]

        if make_labels: ax_h.legend(plot_data, [strain[0],strain[strain_num]])
        ax_h.set_xlabel('Normalized Time in Adulthood (rel. to max lifespan)')
        ax_h.set_ylabel('Normalized Prognosis')
        fig_data.append([fig_h,ax_h])
        
        if len(out_dir)>0:
            plotting_tools.clean_plot(ax_h, make_labels=make_labels,suppress_ticklabels=not make_labels)
            #~ plotting_tools.clean_plot(ax_h, make_labels=make_labels,suppress_ticklabels=False)
            fig_h.savefig(out_dir+os.path.sep+'health_overall_scaling.png')
    return fig_data

# Parameters scatter (start vs. rate vs. death)
def parameter_analysis(strain_dfs, cohort_info=None, out_dir='',make_labels=None,strain_colors=None,plot_trenddata=True,orient='horiz'):
    if strain_colors == None: strain_colors = 4*[None]
    
    # Use spe-9 percentile bins
    animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins)
    
    if orient=='horiz':
        par_fig, ax_h = plt.subplots(3,len(strain_dfs))
        ax_h = ax_h.T
    elif orient =='vert':
        par_fig, ax_h = plt.subplots(len(strain_dfs),3)
    for strain_health, strain_axs,strain_color in zip(strain_dfs, ax_h,strain_colors):
        my_adultspans = analyzeHealth.selectData.get_adultspans(strain_health)/24  
        
        geometry_dict = analyzeHealth.computeStatistics.one_d_geometries(strain_health, 'health')
        start_data = geometry_dict['start']
        end_data = geometry_dict['end']
        rate_data = (start_data - end_data)/my_adultspans
        
        # Make scatters
        graphingFigures.cannedFigures.cohort_scatters(strain_axs[0], my_adultspans, start_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile',the_title = 'Start', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Starting Prognosis (Remaining Days)', label_coordinates = (4, 5),no_cohorts_color=strain_color,s=1**2,plot_trenddata=plot_trenddata)
        graphingFigures.cannedFigures.cohort_scatters(strain_axs[1], my_adultspans, rate_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Rate', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Aging Rate (Dimensionless)', label_coordinates = (10, 1.5), polyfit_degree = 2,no_cohorts_color=strain_color,s=1**2,plot_trenddata=plot_trenddata)
        graphingFigures.cannedFigures.cohort_scatters(strain_axs[2], my_adultspans, end_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'End', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Ending Prognosis (Remaining Days)', label_coordinates = (4, 6), polyfit_degree = 2,no_cohorts_color=strain_color,s=1**2,plot_trenddata=plot_trenddata)
    if orient == 'horiz': [plotting_tools.force_same_plot_attributes(axs,'ylim') for axs in ax_h.T]
    elif orient == 'vert': [plotting_tools.force_same_plot_attributes(axs,'ylim') for axs in ax_h.T]
    
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        par_fig.savefig(out_dir+os.path.sep+'par_analysis.png')
    return (par_fig,ax_h)

def deviation_analysis(strain_dfs, anal_mode='absolute', cohort_info=None,out_dir='',make_labels=True):
    strain_colors=[[0,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0.5]
    ]
    
    # Use percentile bins
    animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
    
    dev_fig, ax_h = plt.subplots(2, len(strain_dfs))
    
    for strain_num, (strain_health, strain_axs) in enumerate(zip(strain_dfs, ax_h.T)):
        my_adultspans = analyzeHealth.selectData.get_adultspans(strain_health)/24  

        # Prepare my "inflection" data. 
        geometry_dict = analyzeHealth.computeStatistics.one_d_geometries(strain_health, 'health')
        start_data = geometry_dict['start']
        mean_start = np.mean(start_data)            
        inflection_data = geometry_dict['absolute_inflection']
        relative_inflection = geometry_dict['self_inflection']
        
        if anal_mode is 'absolute':
            # Plot the traces and scatter for absolute inflection.
            graphingFigures.cannedFigures.cohort_traces(strain_axs[0], 'health', strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Prognosis Over Normalized Time', the_xlabel = 'Fractional Adult Lifespan', the_ylabel = 'Prognosis (Remaining Days)', x_normed = True, make_labels=make_labels)
            strain_axs[0].set_ylim([0, 1.1*mean_start])
            graphingFigures.cannedFigures.cohort_scatters(strain_axs[1], my_adultspans, inflection_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Absolute Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Days)', label_coordinates = (12, 2.5), make_labels=make_labels)
        elif anal_mode is 'relative':
            # Plot the traces and scatter for relative inflection.
            graphingFigures.cannedFigures.cohort_traces(strain_axs[0], 'health', strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Relative Prognosis Over Normalized Time', the_xlabel = 'Fractional Adult Lifespan', the_ylabel = 'Relative Prognosis (Fractional Remaining Life)', x_normed = True, y_normed = True, make_labels=make_labels)
            strain_axs[0].set_ylim([-0.1, 1.1])
            graphingFigures.cannedFigures.cohort_scatters(strain_axs[1], my_adultspans, relative_inflection, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Relative Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Relative Prognosis)', label_coordinates = (4, -0.4),make_labels=make_labels)
    plotting_tools.force_same_plot_attributes(ax_h[1,:],'ylim')
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        dev_fig.savefig(out_dir+os.path.sep+anal_mode+'_deviation.png')
    return (dev_fig, ax_h)
    
def deviation_rescaling(strain_dfs, anal_mode='absolute', cohort_info=None,out_dir='',make_labels=True,color_mode='cohort'):
   
    # Use percentile bins
    animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
    
    dev_fig, ax_h = plt.subplots(1,1)
    
    strain_adultspans = []
    for strain_num, strain_health in enumerate(strain_dfs):
        my_adultspans = analyzeHealth.selectData.get_adultspans(strain_health)/24  

        # Prepare my "inflection" data. 
        geometry_dict = analyzeHealth.computeStatistics.one_d_geometries(strain_health, 'health')
        start_data = geometry_dict['start']
        mean_start = np.mean(start_data)            
        inflection_data = geometry_dict['absolute_inflection']
        relative_inflection = geometry_dict['self_inflection']
        
        if anal_mode is 'absolute':
            # Plot the traces and scatter for absolute inflection.
            #~ graphingFigures.cannedFigures.cohort_scatters(ax_h, my_adultspans/my_adultspans.max(), inflection_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Absolute Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Days)', label_coordinates = (12, 2.5), make_labels=make_labels,no_cohorts_color=plotting_tools.qual_colors[strain_num],plot_trenddata=False)
            graphingFigures.cannedFigures.cohort_scatters(ax_h, my_adultspans, relative_inflection, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Absolute Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Relative Deviation', label_coordinates = (12, 2.5), make_labels=make_labels,no_cohorts_color=plotting_tools.qual_colors[strain_num],plot_trenddata=False)
        elif anal_mode is 'relative':
            # Plot the traces and scatter for relative inflection.
            graphingFigures.cannedFigures.cohort_scatters(ax_h, my_adultspans/np.median(my_adultspans), relative_inflection, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Relative Deviation', the_xlabel = 'Adultspan/Pop Median', the_ylabel = 'Average Deviation (Relative Prognosis)', label_coordinates = (4, -0.4),make_labels=make_labels,no_cohorts_color=plotting_tools.qual_colors[strain_num],plot_trenddata=False)
        strain_adultspans.append(my_adultspans/my_adultspans.max())
    #ax_h.set_title('r^2:{:3f}\np:{:3f}'.format(*scipy.stats.pearsonr(*strain_adultspans)))
    
    if len(out_dir)>0:
        plotting_tools.clean_plot(ax_h, make_labels=make_labels,suppress_ticklabels=not make_labels)
        dev_fig.savefig(out_dir+os.path.sep+anal_mode+'_deviationrescaling.png')
    return (dev_fig, ax_h)

def span_analysis(strain_dfs, strains, out_dir='', make_labels=True, cutoff_type=None,plot_mode='separate', bin_info={'bin_mode':'WT'}):
    '''
        Look at the health/gerospans of each strain based on Willie's old span finding code
    '''
    
    # Use percentile bins
    if bin_info['bin_mode'] == 'WT':
        animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
        animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
    elif bin_info['bin_mode'] in ['percent','percentile']:
        ntiles=bin_info['ntiles']
        animal_bins = [(i+1)*100/ntiles for i in range(ntiles)]

    
    if cutoff_type is 'global':
        flat_data = np.concatenate(np.array([np.ndarray.flatten(strain_df.mloc(strain_df.worms,['health'])) for strain_df in strain_dfs])).ravel()
        flat_data = flat_data[~np.isnan(flat_data)]
        my_cutoff = np.percentile(flat_data, 0.5*100)
        my_cutoff = strain_dfs[0].display_variables(my_cutoff, 'health')[0]
    elif cutoff_type is 'WT':
        flat_data = np.ndarray.flatten(strain_dfs[0].mloc(strain_dfs[0].worms, ['health']))
        flat_data = flat_data[~np.isnan(flat_data)]
        my_cutoff = np.percentile(flat_data, 0.5*100)
        my_cutoff = strain_dfs[0].display_variables(my_cutoff, 'health')[0]
    elif cutoff_type is None:
        my_cutoff = None
    else:
        raise Exception('(span_analysis) Bad cutofftype')
        
        
    span_figs = []
    span_axs = []
    if plot_mode is 'separate':
        for strain, strain_health in zip(strains,strain_dfs):
            span_fig, ax_h = plt.subplots(2,2)
            graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
                strain_health,make_labels=make_labels,
                cohort_info=analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = animal_bins,bin_mode='percentile'),
                my_cutoff=my_cutoff,
                bar_plot_mode='uniform')
            span_figs.append(span_fig)
            span_axs.append(ax_h)
            if len(out_dir)>0:
                [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
                span_fig.savefig(out_dir+os.path.sep+strain+'_healthspans.png')
    elif plot_mode is 'comb_collapsemutant':
        span_fig, ax_h = plt.subplots(2,2)
        for strain_num,(strain, strain_health,strain_color) in enumerate(zip(strains,strain_dfs,plotting_tools.qual_colors[:len(strains)])):
            if strain_num is 0:
                graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
                    strain_health,make_labels=make_labels,
                    cohort_info=analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = animal_bins,bin_mode='percentile'),
                    my_cutoff=my_cutoff,
                    bar_plot_mode='uniform')
            else:
                cohort_info = analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = np.array([100]),bin_mode='percentile')
                #cohort_info[3] = np.array([[42/255,160/255,75/255]])
                cohort_info[3] = np.array([strain_color])
                graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
                    strain_health,make_labels=make_labels,
                    cohort_info=cohort_info,
                    my_cutoff=my_cutoff,
                    bar_plot_mode='uniform',
                    stop_with_death=False)
        span_axs.append(ax_h)
        span_figs.append(span_fig)
        if len(out_dir)>0:
            [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
            #~ [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=False) for my_ax in ax_h.flatten()]
            span_fig.savefig(out_dir+os.path.sep+'comb_collapsemutant_healthspans.png')
            span_fig.savefig(out_dir+os.path.sep+'comb_collapsemutant_healthspans.svg')

    return (span_figs,span_axs)

def span_comparison(strain_dfs,cutoff_type=None,plot_layout='combined',time_mode='absolute',plot_var = 'gerospan', relative_time = 0.5):    
    '''
        Generates a scatter of health- or gerospans against adultspan (of individuals) in each supplied populations; also generates a regression curve using lowess
    '''
    
    if cutoff_type is 'global':
        flat_data = np.concatenate(np.array([np.ndarray.flatten(strain_df.mloc(strain_df.worms,['health'])) for strain_df in strain_dfs])).ravel()
        flat_data = flat_data[~np.isnan(flat_data)]
        my_cutoff = np.percentile(flat_data, relative_time*100)
        my_cutoff = strain_dfs[0].display_variables(my_cutoff, 'health')[0]
    elif cutoff_type is 'WT':
        flat_data = np.ndarray.flatten(      # Either one is the same but units and sign are different
            -1*strain_dfs[0].mloc(measures=['health'])[:,0,:])
        flat_data = flat_data[~np.isnan(flat_data)]
        my_cutoff = np.percentile(flat_data, (1-relative_time)*100)
    elif cutoff_type is None:
        my_cutoff = None
    else:
        raise Exception('(span_analysis) Bad cutoff type')
    
    def plot_helper(strain_df,strain_color, plot_layout, time_mode,plot_var,relative_time, my_cutoff,ax_h):
        adultspans = analyzeHealth.selectData.get_adultspans(strain_df)/24
        healthspans = analyzeHealth.computeStatistics.get_spans(strain_df,'health', 'overall_time',fraction=0.5, cutoff_value = my_cutoff,reverse_direction=True)/24
        
        if plot_var == 'gerospan':
            y_vals = (adultspans-healthspans)/adultspans
        elif plot_var == 'healthspan':
            y_vals = healthspans/adultspans
        
        if time_mode == 'absolute':
            t_vals = adultspans
        elif time_mode == 'relative':
            t_vals = adultspans/np.median(adultspans)
            
        ax_h.scatter(t_vals,y_vals,color=strain_color)
        reg_curve = smoothing.lowess(t_vals,y_vals)

        sorted_adultspan_idx = sorted(range(len(t_vals)), key = lambda k:t_vals[k])
        ax_h.plot(t_vals[sorted_adultspan_idx],reg_curve[sorted_adultspan_idx],color=strain_color,linewidth=2.5)
    
        if plot_var == 'gerospan':
            ax_h.set_ylabel('Gerospan (Fraction of Life)')
        elif plot_var == 'healthspan':
            ax_h.set_ylabel('Healthspan (Fraction of Life)')
        ax_h.set_ylim([0,1.4])
        
        if time_mode == 'absolute': ax_h.set_xlabel('Adultspan')
        elif time_mode == 'relative': ax_h.set_xlabel('Adultspan rel. to population median')
        else: raise Exception('(span_comparison) Bad time mode')
        
        return ax_h
        
    
    if plot_layout is 'combined':
        health_scatter, ax_h = plt.subplots(1,1)
        for strain_df,strain_color in zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)]):
            ax_h = plot_helper(strain_df,strain_color, plot_layout, time_mode,plot_var,relative_time,my_cutoff,ax_h)
        return health_scatter, ax_h
    elif plot_layout is 'separate':
        #health_scatter,axs = plt.subplots(1,len(strain_dfs))
        figs, axs = [], []
        for strain_df,strain_color in zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)]):
            health_scatter, ax_h = plt.subplots(1,1)
            ax_h = plot_helper(strain_df,strain_color, plot_layout, time_mode,plot_var,relative_time,my_cutoff,ax_h)
            figs.append(health_scatter)
            axs.append(ax_h)
        return figs, axs

def avghealth_comparison(strain_dfs, plot_layout='combined',time_mode='absolute'):
    '''
        Generates a scatter of average health across life against adultspan (of individuals) in each supplied population; also generates a regression curve using lowess
    '''
    
    def plot_helper(strain_df, strain_color, time_mode, plot_layout,ax_h):
        adultspans = analyzeHealth.selectData.get_adultspans(strain_df)/24
        
        pop_data = strain_df.mloc(strain_df.worms,['health'])[:,0,:]
        avg_health = np.nanmean(pop_data,axis=1)
        (avg_health,trash,trash) = strain_df.display_variables(avg_health, 'health')
        
        if time_mode == 'absolute':
            t_vals = adultspans
        elif time_mode == 'relative':
            t_vals = adultspans/np.median(adultspans)
        
        ax_h.scatter(t_vals,avg_health,color=strain_color)
        reg_curve = smoothing.lowess(t_vals,avg_health)
        
        sorted_adultspan_idx = sorted(range(len(t_vals)), key = lambda k:t_vals[k])
        ax_h.plot(t_vals[sorted_adultspan_idx],reg_curve[sorted_adultspan_idx],color=strain_color,linewidth=2.5)
        ax_h.set_ylabel('Average prognosis (days remaining)')
        ax_h.set_ylim([0,10])
        
        if time_mode == 'absolute': ax_h.set_xlabel('Adultspan')
        elif time_mode == 'relative': ax_h.set_xlabel('Adultspan rel. to population median')
        else: raise Exception('(span_comparison) Bad time mode')
        
        return ax_h
    
    if plot_layout is 'combined':
        health_scatter, ax_h = plt.subplots(1,1)
        for strain_df,strain_color in zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)]):
            ax_h = plot_helper(strain_df,strain_color, time_mode, plot_layout,ax_h)
        return health_scatter, ax_h
    elif plot_layout is 'separate':
        #health_scatter,axs = plt.subplots(1,len(strain_dfs))
        figs, axs = [], []
        for strain_df,strain_color in zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)]):
            health_scatter, ax_h = plt.subplots(1,1)
            ax_h = plot_helper(strain_df,strain_color, time_mode, plot_layout,ax_h)
            figs.append(health_scatter)
            axs.append(ax_h)
        return figs, axs
        
def healthagainstdev_comparison(strain_dfs, time_mode = 'absolute'):
    '''
        Generates a scatter of average health across life against deviation in each supplied population (as a sanity check - should be correlated decently well....
    '''
    
    health_scatter, ax_h = plt.subplots(1,1)
    for strain_df,strain_color in zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)]):
        pop_data = strain_df.mloc(strain_df.worms,['health'])[:,0,:]
        avg_health = np.nanmean(pop_data,axis=1)
        (avg_health,trash,trash) = strain_df.display_variables(avg_health, 'health')
        
        geometry_dict = analyzeHealth.computeStatistics.one_d_geometries(strain_df, 'health')
        start_data = geometry_dict['start']
        mean_start = np.mean(start_data)            
        inflection_data = geometry_dict['absolute_inflection']
        relative_inflection = geometry_dict['self_inflection']
        
        if time_mode == 'absolute':
            ax_h.scatter(avg_health,inflection_data,color=strain_color)
        elif time_mode == 'relative':
            ax_h.scatter(avg_health,relative_inflection,color=strain_color)
    ax_h.set_xlabel('Average Health')
    if time_mode == 'absolute':
        ax_h.set_ylabel('Absolute Deviation')
    elif time_mode == 'relative':
        ax_h.set_ylabel('Normalized Deviation')
    return health_scatter, ax_h

def cohort_percentiles(strain_dfs,time_mode='absolute',plot_var='gerospan',ntiles=5,bin_mode='percent'):    
    '''
        Generates a scatter of gerospan or average health across life (for percentile cohorts) against adultspan in each supplied population
    '''
    
    fig_h, ax_h = plt.subplots(1,1)
    
    for strain_num, (strain_df,strain_color) in enumerate(zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)])):
        adultspans = analyzeHealth.selectData.get_adultspans(strain_df)/24
        if bin_mode == 'percent':
            animal_bins = analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days = [(i+1)*100/ntiles for i in range(ntiles)],bin_mode='percentile')[0]
        elif bin_mode == 'WT':
            WT_binsizes = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
            WT_cumpercent = 100*np.cumsum(WT_binsizes)/sum(WT_binsizes) # Make percentile bins
            animal_bins = analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days = WT_cumpercent, bin_mode='percentile')[0]
            
        
        if plot_var == 'gerospan':
            flat_data = np.ndarray.flatten(
                -1*strain_dfs[0].mloc(measures=['health'])[:,0,:])
            flat_data = flat_data[~np.isnan(flat_data)]
            my_cutoff = np.percentile(flat_data, (1-0.5)*100)
            
            healthspans = analyzeHealth.computeStatistics.get_spans(strain_df,'health', 'overall_time',fraction=0.5, cutoff_value = my_cutoff,reverse_direction=True)/24
            gerospans = (adultspans-healthspans)
            
            if strain_num == 0:     
                ax_h.scatter(adultspans,gerospans/adultspans,color=(strain_color+5)/6)
            [ax_h.scatter(np.mean(adultspans[cohort_bins]), np.mean(gerospans[cohort_bins]/adultspans[cohort_bins]),color = strain_color) for cohort_bins in animal_bins]
            ax_h.scatter(np.mean(adultspans), np.mean(gerospans/adultspans), color=strain_color,marker='x',s=40)
            #~ cohort_styles = ['.','o','s','D','x']
            #~ [ax_h.scatter(np.mean(adultspans[cohort_bins]), np.mean(gerospans[cohort_bins]/adultspans[cohort_bins]),color = strain_color,marker=cohort_style,s=40) for cohort_style,cohort_bins in zip(cohort_styles,animal_bins)]
        elif plot_var == 'avg_health':
            pop_data = strain_df.mloc(strain_df.worms,['health'])[:,0,:]
            avg_health = np.nanmean(pop_data,axis=1)
            (avg_health,trash,trash) = strain_df.display_variables(avg_health, 'health')
            
            if strain_num == 0:     
                ax_h.scatter(adultspans,avg_health,color=(strain_color+10)/11)
            [ax_h.scatter(np.mean(adultspans[cohort_bins]), np.mean(avg_health[cohort_bins]), color=strain_color) for cohort_bins in animal_bins]
            ax_h.scatter(np.mean(adultspans), np.mean(avg_health), color=strain_color,marker='x',s=40)

        
    ax_h.set_xlabel('Adultspan (d)')
    if plot_var == 'gerospan': 
        ax_h.set_ylabel('Gerospan')
        ax_h.set_ylim([0,1])
    elif plot_var == 'avg_health': ax_h.set_ylabel('Average Health')
    
    
    return (fig_h, ax_h)

def cohort_percentiles_oldadapted(strain_dfs,ntiles=5,mean_analysis='ind'):    
    '''
        Generates a scatter of gerospan for percentile cohorts against adultspan using a Willie-ish apporach to calculating gerospan in each supplied population
        Helpful for debugging new analyses against Willie's old approach
    '''
    
    fig_h, ax_h = plt.subplots(1,1)
    
    flat_data = np.ndarray.flatten(
        -1*strain_dfs[0].mloc(measures=['health'])[:,0,:])
    flat_data = flat_data[~np.isnan(flat_data)]
    my_cutoff_raw = np.percentile(flat_data, (1-0.5)*100)
    my_cutoff = strain_dfs[0].display_variables(my_cutoff_raw,'health')[0]
    
    strain_transitions = []
    strain_min_life = []
    for strain_num, (strain_df,strain_color) in enumerate(zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)])):
        adultspans = analyzeHealth.selectData.get_adultspans(strain_df)/24
        
        animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
        animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
        animal_bins = analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days = animal_bins,bin_mode='percentile')[0]
        
        #~ animal_bins = analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days = [(i+1)*100/ntiles for i in range(ntiles)],bin_mode='percentile')[0]
        
        healthspans_ind = analyzeHealth.computeStatistics.get_spans(strain_df,'health', 'overall_time',fraction=0.5, cutoff_value = my_cutoff_raw,reverse_direction=True)/24
        gerospans_ind = (adultspans-healthspans_ind)
        
        cohort_transitions = []
        minimum_life_cohorts = []
        for my_cohort in animal_bins:
            if len(my_cohort) > 0:
                cohort_data = strain_df.mloc(strain_df.worms, ['health'])[my_cohort, 0, :]
                cohort_data = cohort_data[~np.isnan(cohort_data).all(axis = 1)]
                cohort_data = np.mean(cohort_data, axis = 0)
                (cohort_data, my_unit, fancy_name) = strain_df.display_variables(cohort_data, 'health')
                
                cohort_ages = np.array(strain_df.ages[:cohort_data.shape[0]])            

                healthy_mask = (cohort_data > my_cutoff)
                first_unhealthy_index = healthy_mask.argmin()
                unhealthy_mask = (cohort_data < my_cutoff)          
                unhealthy_mask[first_unhealthy_index - 1] = True
                minimum_life_cohorts.append(cohort_ages[unhealthy_mask][-1])

                if healthy_mask.all():
                    cohort_transitions.append(cohort_ages[-1])                  
                else:
                    cohort_transitions.append(cohort_ages[first_unhealthy_index])
        
        
        if strain_num == 0:     
            ax_h.scatter(adultspans,gerospans_ind/adultspans,color=(strain_color+10)/11)    # This is potentially misleading though.... maybe.....
        
        [ax_h.scatter(np.mean(adultspans[cohort_bins]), 1-(cohort_hs/cohort_ml),color = strain_color) for cohort_bins, cohort_hs,cohort_ml in zip(animal_bins,cohort_transitions, minimum_life_cohorts)]    # Individual-level
        #print([cohort_hs/cohort_ml for cohort_bins, cohort_hs,cohort_ml in zip(animal_bins,cohort_transitions, minimum_life_cohorts)])
        
        if mean_analysis == 'ind':
            # Do the grand mean of the individual data
            ax_h.scatter(np.mean(adultspans), np.mean(gerospans_ind/adultspans), color=strain_color,marker='x',s=40)
        elif mean_analysis in ['pop','population']:   
            cohort_data = strain_df.mloc(strain_df.worms, ['health'])[:,0,:]
            cohort_data = cohort_data[~np.isnan(cohort_data).all(axis = 1)]
            cohort_data = np.nanmean(cohort_data, axis = 0)
            (cohort_data, my_unit, fancy_name) = strain_df.display_variables(cohort_data, 'health')
                
            cohort_ages = np.array(strain_df.ages[:cohort_data.shape[0]])  
            
            healthy_mask = (cohort_data > my_cutoff)
            first_unhealthy_index = healthy_mask.argmin()
            unhealthy_mask = (cohort_data < my_cutoff)          
            unhealthy_mask[first_unhealthy_index - 1] = True
            minimum_life_pop = cohort_ages[unhealthy_mask][-1]
            
            if healthy_mask.all():
                pop_transition = cohort_ages[-1]
            else:
                pop_transition = cohort_ages[first_unhealthy_index]
            
            ax_h.scatter(np.mean(adultspans), 1-(pop_transition/minimum_life_pop), color=strain_color, marker='x',s=40)
        
    ax_h.set_xlabel('Adultspan (d)')
    ax_h.set_ylabel('Fractional Gerospan')
    ax_h.set_ylim([0,1])
    
    return (fig_h, ax_h)

def get_healthspans(adult_df, a_variable='health',cutoff_value=None,return_crossings=False,temp=None):
    '''
        Get healthspans based on dwell time under threshold (more robust than Willie spans, part. to end of life noise)
        Cutoff (i.e. everything) should be raw!! (i.e. not adjusted with CompleteDF.display_variables)
    '''
    
    unit_multipliers = {    #
            'intensity_90': None, 
            'intensity_80': None, 
            'cumulative_eggs': 1,
            'cumulative_eggs_rate': 1/3,
            'cumulative_area': (1.304/1000)**2,
            'visible_eggs': 1,
            'total_size': (1.304/1000)**2, 
            'age_texture': 1, 
            'bulk_movement': (1.304/1000)/3,
            'stimulated_rate_a': (1.304/1000),
            'stimulated_rate_b': (1.304/1000),
            'unstimulated_rate': (1.304/1000),
            'area': 1/100000,
            'life_texture': -1,
            'adjusted_size': (1.304/1000)**2,
            'adjusted_size_rate': ((1.304/1000)**2)/3,
            'great_lawn_area': (1.304/1000)**2, 
            'texture': (-1/24),
            'eggs': (-1/24),
            'autofluorescence': (-1/24),
            'movement': (-1/24),
            'size': (-1/24),
            'health': (-1/24),
    }
    measures_to_negate = ['intensity_90','intensity_80', 'life_texture', 'autofluorescence','movement','size','health','eggs']
    # TODO - think about how to span size.... Need to change anything?
    
    data_values = adult_df.mloc(measures=[a_variable])[:,0,:]
    if a_variable in measures_to_negate:
        data_values = data_values*-1 # Reverse direction (unit_multiplier for these is negative)
    if cutoff_value is None:
        all_data = np.ndarray.flatten(data_values*-1)
        if a_variable in measures_to_negate:
            all_data = all_data*-1
        all_data = all_data[~np.isnan(all_data)]
        cutoff_value = np.percentile(all_data, 0.5*100)
    adultspans = analyzeHealth.selectData.get_adultspans(adult_df)
    ghost_ages = adult_df.mloc(measures=['ghost_age'])[:,0,:]
    
    healthspans = []
    crossing_idxs = []
    for my_worm,worm_data,adultspan,ghost_age in zip(adult_df.worms,data_values,adultspans,ghost_ages):
        # Get all crossings and isolate those that are going down
        adj_data = worm_data-cutoff_value
        adj_data = adj_data[~np.isnan(adj_data)]
        crossings = np.where((adj_data[:-1]>0) & (adj_data[1:]<0))[0]
        
        adultspan_len = len(ghost_age[~np.isnan(ghost_age)])
        
        # Handle no crossings
        if len(crossings) == 0:
            #print('no crossings found')
            if adj_data[0]>0:
                healthspans.append(adultspan)
                crossing_idxs.append(adultspan_len)
            elif adj_data[0]<=0:
                healthspans.append(0)
                crossing_idxs.append(0)
        else:
            # Get the first crossing that lingers below the cutoff for more than 10% of lifetime
            found_crossing = False
            for crossing_idx in crossings:
                if (adj_data[crossing_idx+1:np.floor(crossing_idx+0.1*adultspan_len)+1]<0).all():
                    healthspans.append(adult_df.ages[crossing_idx]*24)
                    crossing_idxs.append(crossing_idx)
                    found_crossing=True
                    break
            if not found_crossing:
                #print('couldnt find crossing for '+my_worm)
                healthspans.append(adult_df.ages[crossing_idx]*24) # Default to take the last crossing in a bad situation
                crossing_idxs.append(crossing_idx)
    
    
    if not return_crossings:
        return np.array(healthspans)
    else:
        return np.array(healthspans),np.array(crossing_idxs)

def get_strain_dist(strain_df, plot_var='gerospan', health_var='health',**kw_data):
    '''
        Build kde estimates of various span constructs
    '''
    adultspans = analyzeHealth.selectData.get_adultspans(strain_df)/24
    healthspans = get_healthspans(strain_df,health_var, cutoff_value = kw_data['cutoff_value'])/24
    
    if plot_var == 'gerospan':
        gerospans = (adultspans-healthspans)
        gs_support, gs_density, kde_obj = zplib.scalar_stats.kde.kd_distribution(gerospans)
        return [gs_support, gs_density]
    elif plot_var == 'fgerospan':
        fgerospans = (adultspans-healthspans)/adultspans
        hs_support, hs_density, kde_obj = zplib.scalar_stats.kde.kd_distribution(fgerospans)
        return [hs_support, hs_density]
    elif plot_var == 'healthspan':
        hs_support, hs_density, kde_obj = zplib.scalar_stats.kde.kd_distribution(healthspans)
        return [hs_support, hs_density]
    elif plot_var == 'fhealthspan':
        hs_support, hs_density, kde_obj = zplib.scalar_stats.kde.kd_distribution(healthspans/adultspans)
        return [hs_support, hs_density]
    elif plot_var == 'avg_health':
        pop_data = strain_df.mloc(strain_df.worms,['health'])[:,0,:]
        avg_health = np.nanmean(pop_data,axis=1)
        (avg_health,trash,trash) = strain_df.display_variables(avg_health, 'health')
        
        health_support, health_density, kde_obj = zplib.scalar_stats.kde.kd_distribution(average_health)
        return [health_support, health_density]
        
def plot_health_dist(strain_dfs, plot_var,health_var='health', my_cutoff = None, ax_h = None):
    if my_cutoff is None: # Use WT cutoff
        flat_data = np.ndarray.flatten(
            -1*strain_dfs[0].mloc(measures=[health_var])[:,0,:])
        flat_data = flat_data[~np.isnan(flat_data)]
        my_cutoff = np.percentile(flat_data, (1-0.5)*100)
        
    if ax_h is None:
        fig_h, ax_h = plt.subplots(1,1)
        ax_provided = False
    else: ax_provided = True

    for strain_df,strain_color in zip(strain_dfs,plotting_tools.qual_colors[:len(strain_dfs)]):
        dist_x, dist_y = get_strain_dist(strain_df,plot_var = plot_var,health_var=health_var,cutoff_value = my_cutoff)
        ax_h.plot(dist_x,dist_y,color=strain_color,linewidth=2)

    if plot_var == 'gerospan': 
        ax_h.set_xlabel('Gerospan (days)')
        ax_h.set_xlim([0,ax_h.get_xlim()[1]])
    elif plot_var == 'fgerospan': 
        ax_h.set_xlabel('Fractional Gerospan')
        ax_h.set_xlim([0,1])
    elif plot_var =='healthspan':
        ax_h.set_xlabel('Healthspan (days)')
    elif plot_var == 'fhealthspan':
        ax_h.set_xlabel('Fractional Healthspan')
        ax_h.set_xlim([0,1])
    elif plot_var == 'avg_health':
        ax_h.set_xlabel('Average Health')
    
    ax_h.set_ylabel('Density')
    
    if ax_provided:
        return ax_h
    else:
        return (fig_h,ax_h)

if __name__ is "__main__":
    make_labels= False
    #strains = ['spe-9','age-1']
    #out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/processed_data/age-1_cohorts+regression_20160818/'
    #out_dir = ''
    #strains = ['spe-9percombinedweightedSVM','age-1percombinedweightedSVM']
    #~ strains = ['spe-9','age-1specific']
    strains = ['spe-9','age-1','ife-2','clk-1']
    plot_data = False
    out_dir='/media/Data/Documents/Presentations/LabMeeting_20170126/Analysis_Figures_WTclassifier/'
    out_dir = ''
    if make_labels and out_dir is not '': 
        out_dir = out_dir+os.path.sep+'labeled'+os.path.sep
        if not os.path.isdir(out_dir): os.mkdir(out_dir)
    
    #~ strain_dfs = load_strain_data(strains)
    strain_data = load_strain_data(strains,'/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/dfs/WT_classifier/')
    strain_dfs = [data['adult_df'] for data in strain_data]
    if plot_data: 
        plot_survival(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels)
        
        #plot_strain_health(strain_dfs,group_mode='population',make_labels=make_labels,out_dir=out_dir)
        plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir) # Plot all spe-9 cohorts with all cohorts from all mutants superimposed
        plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir,var_plot_mode='separate') #Plot all spe-9 cohorts and separate health variable plots by variable
        plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir,collapse_mutant=True) #  Plot all spe-9 cohorts with mutant population collapsed over it
        plot_strain_health(strain_dfs,group_mode='spe9_3cohorts',make_labels=make_labels,out_dir=out_dir)   # Plot 3 cohorts of based on spe-9 cohorts
        plot_strain_health(strain_dfs,group_mode='equalbins_7',make_labels=make_labels,out_dir=out_dir) # Plot 7 cohorts based on equal binning on eventual lifespan
        
        parameter_analysis(strain_dfs,
            out_dir=out_dir,make_labels=make_labels)
        
        [abs_fig, abs_ax] = deviation_analysis(strain_dfs,'absolute',
            out_dir=out_dir,make_labels=make_labels)
        [rel_fig, rel_ax] = deviation_analysis(strain_dfs,'relative',
            out_dir=out_dir,make_labels=make_labels)
        
        span_analysis(strain_dfs, strains, out_dir=out_dir,make_labels=make_labels,cutoff_type='WT')
        span_analysis(strain_dfs, strains, out_dir=out_dir,make_labels=make_labels,plot_mode='comb_collapsemutant',cutoff_type='WT') # Put all the mutants onto one plot

        
        #~ [abs_fig, abs_ax] = process_mutant_data.deviation_rescaling(strain_dfs,'absolute',      
            #~ out_dir=out_dir,make_labels=make_labels)
        #~ [abs_fig, abs_ax] = process_mutant_data.deviation_rescaling(strain_dfs,'relative',      
            #~ out_dir=out_dir,make_labels=make_labels)
        [abs_fig, abs_ax] = process_mutant_data.deviation_rescaling(strain_dfs,'absolute',      
            out_dir='',make_labels=True)
        [rel_fig, rel_ax] = process_mutant_data.deviation_rescaling(strain_dfs,'relative',      
            out_dir='',make_labels=True)
        if len(out_dir) > 0:
            [plotting_tools.clean_plot(ax_h,cleaning_mode = 'PPT' if not make_labels else 'verbose') for ax_h in [abs_ax, rel_ax]]
            [fig_h.savefig(out_dir+os.path.sep+'deviation_compare_'+time_label+'.png') for fig_h,time_label in zip([abs_fig,rel_fig],['abs','rel'])]
            [fig_h.savefig(out_dir+os.path.sep+'deviation_compare_'+time_label+'.svg') for fig_h,time_label in zip([abs_fig,rel_fig],['abs','rel'])]
        
        [span_figs, span_axs] = process_mutant_data.span_comparison(strain_dfs,cutoff_type='WT',plot_var='gerospan',time_mode='absolute',plot_layout='separate')
        if len(out_dir) > 0:
            [plotting_tools.clean_plot(ax_h,cleaning_mode = 'PPT' if not make_labels else 'verbose') for ax_h in span_axs]
            [fig_h.savefig(out_dir+os.path.sep+'span_compare_'+strain+'_abs.png') for fig_h, strain in zip(span_figs,strains)]
            [fig_h.savefig(out_dir+os.path.sep+'span_compare_'+strain+'_abs.svg') for fig_h, strain in zip(span_figs,strains)]
        [fig_h, ax_h] = process_mutant_data.span_comparison(strain_dfs,cutoff_type='WT',plot_var='gerospan',time_mode='absolute',plot_layout='combined')
        if len(out_dir) > 0:
            plotting_tools.clean_plot(ax_h,cleaning_mode = 'PPT' if not make_labels else 'verbose')
            fig_h.savefig(out_dir+os.path.sep+'span_compare_combined_abs.png')
            fig_h.savefig(out_dir+os.path.sep+'span_compare_combined_abs.svg')
        [span_figs, span_axs] = process_mutant_data.span_comparison(strain_dfs,cutoff_type='WT',plot_var='gerospan',time_mode='relative',plot_layout='separate')
        if len(out_dir) > 0:
            [plotting_tools.clean_plot(ax_h,cleaning_mode = 'PPT' if not make_labels else 'verbose') for ax_h in span_axs]
            [fig_h.savefig(out_dir+os.path.sep+'span_compare_'+strain+'_rel.png') for fig_h, strain in zip(span_figs,strains)]
            [fig_h.savefig(out_dir+os.path.sep+'span_compare_'+strain+'_rel.svg') for fig_h, strain in zip(span_figs,strains)]
        [fig_h, ax_h] = process_mutant_data.span_comparison(strain_dfs,cutoff_type='WT',plot_var='gerospan',time_mode='relative',plot_layout='combined')
        if len(out_dir) > 0:
            plotting_tools.clean_plot(ax_h,cleaning_mode = 'PPT' if not make_labels else 'verbose')
            fig_h.savefig(out_dir+os.path.sep+'span_compare_combined_rel.png')
            fig_h.savefig(out_dir+os.path.sep+'span_compare_combined_rel.svg')
        
        fig_abs, ax_abs = process_mutant_data.avghealth_comparison(strain_dfs,plot_layout='combined', time_mode='absolute')
        fig_rel, ax_rel = process_mutant_data.avghealth_comparison(strain_dfs,plot_layout='combined', time_mode='relative')
        fig_abs_sep, ax_abs_sep = process_mutant_data.avghealth_comparison(strain_dfs, plot_layout='separate', time_mode='absolute')
        fig_rel_sep, ax_rel_sep = process_mutant_data.avghealth_comparison(strain_dfs, plot_layout='separate', time_mode='relative')
        if len(out_dir) > 0:
            [plotting_tools.clean_plot(my_ax,cleaning_mode='PPT' if not make_labels else 'verbose') for my_ax in [ax_abs,ax_rel]]
            [plotting_tools.clean_plot(my_ax,cleaning_mode='PPT' if not make_labels else 'verbose') for my_ax in ax_abs_sep]
            [plotting_tools.clean_plot(my_ax,cleaning_mode='PPT' if not make_labels else 'verbose') for my_ax in ax_rel_sep]
            fig_abs.savefig(out_dir+os.path.sep+'avghealth_abs_comb.png')
            fig_rel.savefig(out_dir+os.path.sep+'avghealth_rel_comb.png')
            fig_abs.savefig(out_dir+os.path.sep+'avghealth_abs_comb.svg')
            fig_rel.savefig(out_dir+os.path.sep+'avghealth_rel_comb.svg')
            [fig_h.savefig(out_dir+os.path.sep+'avghealth_abs_'+strain+'.png') for fig_h,strain in zip(fig_abs_sep,strains)]
            [fig_h.savefig(out_dir+os.path.sep+'avghealth_rel_'+strain+'.png') for fig_h,strain in zip(fig_rel_sep,strains)]
            [fig_h.savefig(out_dir+os.path.sep+'avghealth_abs_'+strain+'.svg') for fig_h,strain in zip(fig_abs_sep,strains)]
            [fig_h.savefig(out_dir+os.path.sep+'avghealth_rel_'+strain+'.svg') for fig_h,strain in zip(fig_rel_sep,strains)]
        
        abs_fig, abs_h = process_mutant_data.healthagainstdev_comparison(strain_dfs)
        rel_fig, rel_h = process_mutant_data.healthagainstdev_comparison(strain_dfs,time_mode = 'relative')
        if len(out_dir) > 0:
            plotting_tools.clean_plot(abs_h,cleaning_mode='PPT' if not make_labels else 'verbose')
            abs_fig.savefig(out_dir+os.path.sep+'devhealth_abs.png')
            abs_fig.savefig(out_dir+os.path.sep+'devhealth_abs.svg')
            
            plotting_tools.clean_plot(rel_h,cleaning_mode='PPT' if not make_labels else 'verbose')
            abs_fig.savefig(out_dir+os.path.sep+'devhealth_rel.png')
            abs_fig.savefig(out_dir+os.path.sep+'devhealth_rel.svg') 

        #~ plot_strain_rescaling(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels,do_bootstrap=True)
        #~ plot_strain_rescaling(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels,do_bootstrap=True,percent_range=[20,40])
        #~ plot_strain_rescaling(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels,do_bootstrap=True,percent_range=[60,80])
        #~ plt.close('all')
        
        #~ out_dir ='/media/Data/Work/ZPLab/Analysis/MutantHealth/ForZachGrant/'
        #~ process_mutant_data.plot_survival(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels)
        #~ process_mutant_data.plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir,collapse_mutant=True) #  Plot all spe-9 cohorts with mutant population collapsed over it
        #~ process_mutant_data.span_analysis(strain_dfs, strains, out_dir=out_dir,make_labels=make_labels,plot_mode='comb_collapsemutant',cutoff_type='WT') # Put all the mutants onto one plot
        
        #~ for strain,strain_df in zip(strains, strain_dfs):
            #~ fig_h, ax_h = test_lifespan_replicates(strain_df)
            #~ plotting_tools.clean_plot(ax_h, make_labels=make_labels,suppress_ticklabels=False)
            #~ fig_h.savefig(out_dir+os.path.sep+strain+'rep_survival.png')
    
        fig_con, ax_con = process_mutant_data.cohort_percentiles(strain_dfs)
        if len(out_dir) > 0:
            #plotting_tools.clean_plot(ax_con,cleaning_mode='PPT' if not make_labels else 'verbose')
            fig_con.savefig(out_dir+os.path.sep+'continuum_percentiles_alldata.png')
            fig_con.savefig(out_dir+os.path.sep+'continuum_percentiles_alldata.svg')
