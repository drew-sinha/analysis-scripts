import matplotlib.pyplot as plt
import pickle
import numpy as np
import collections
import os
import numpy as np
import scipy.interpolate
import scipy.optimize

import analyzeHealth
import graphingFigures
import plotting_tools

plt.ion()
plt.show()
plt.close('all')

def load_strain_data(strains):
    strain_dfs = []
    for strain in strains:
        #with open('/mnt/bulkdata/wzhang/human_dir/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
        with open('/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
            strain_dfs.append(pickle.load(my_file)['adult_df'])
        print(strain+": n = {}".format(len(strain_dfs[-1].worms)))
        print(strain+": CV={:.2f}".format(
            np.std(analyzeHealth.selectData.get_adultspans(strain_dfs[-1]))/np.mean(analyzeHealth.selectData.get_adultspans(strain_dfs[-1]))))
    return strain_dfs
    
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
        

def plot_strain_health(strain_dfs,group_mode='population',make_labels=True,out_dir='',custom_fname='',collapse_mutant=False,var_plot_mode='combined'):
    default_fnames = {
        'population':'pop',
        'spe9_allcohorts':'allcohorts',
        'spe9_3cohorts':'3cohorts',
        'equalbins_7':'7equalbins',
    }
    if custom_fname is '': custom_fname = default_fnames[group_mode]
    
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
    
    health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
    if var_plot_mode is 'combined':
        var_fig, var_ax = plt.subplots(len(health_vars),1)
        for var_plot,var in zip(var_ax.T,health_vars):
            graphingFigures.cannedFigures.cohort_traces(var_plot,var,strain_dfs[0],make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style='-',cohorts_to_use=cohorts_to_use)
            for strain_health,line_style in zip(strain_dfs[1:],['--']):
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
            for strain_health,line_style in zip(strain_dfs[1:],['--']):
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
    for strain_health,line_style in zip(strain_dfs[1:],['--']):
        if collapse_mutant:
            graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_health,make_labels=make_labels,bin_width_days=np.array([100]),bin_mode='percentile',line_style=line_style,stop_with_death=False)
        else:
            graphingFigures.cannedFigures.cohort_traces(health_ax,'health',strain_health,make_labels=make_labels,bin_width_days=animal_bins,bin_mode='percentile',line_style=line_style)
    if len(out_dir)>0: 
        plotting_tools.clean_plot(health_ax,make_labels,suppress_ticklabels=not make_labels)
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
    print(bin_id)
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
            fig_h.savefig(out_dir+os.path.sep+'health_overall_scaling.png')
    return fig_data
    

def plot_survival(strain_dfs, out_dir='', make_labels=True):
    # Plot survival curves and lifespan distributions for each strain
    
    # Use percentile bins
    animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
    
    survival_fig, ax_h = plt.subplots(len(strains),2)
    for strain_health, strain_ax in zip(strain_dfs,ax_h.T):
        graphingFigures.cannedFigures.survival_lifespan(strain_ax[0],strain_ax[1],strain_health, make_labels=make_labels,
            cohort_info=analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = animal_bins,bin_mode='percentile'))
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        survival_fig.savefig(out_dir+os.path.sep+'survival_curves.png')
    return (survival_fig, ax_h)

# Parameters scatter (start vs. rate vs. death)
def parameter_analysis(strain_dfs, cohort_info=None, out_dir='',make_labels=None,):
    # Use spe-9 percentile bins
    animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins)
    
    par_fig, ax_h = plt.subplots(3,len(strain_dfs))
    for strain_health, strain_axs in zip(strain_dfs, ax_h.T):
        my_adultspans = analyzeHealth.selectData.get_adultspans(strain_health)/24  
        
        geometry_dict = analyzeHealth.computeStatistics.one_d_geometries(strain_health, 'health')
        start_data = geometry_dict['start']
        end_data = geometry_dict['end']
        rate_data = (start_data - end_data)/my_adultspans
        
        # Make scatters
        graphingFigures.cannedFigures.cohort_scatters(strain_axs[0], my_adultspans, start_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile',the_title = 'Start', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Starting Prognosis (Remaining Days)', label_coordinates = (4, 7.5))
        graphingFigures.cannedFigures.cohort_scatters(strain_axs[1], my_adultspans, rate_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Rate', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Aging Rate (Dimensionless)', label_coordinates = (10, 1.5), polyfit_degree = 2)
        graphingFigures.cannedFigures.cohort_scatters(strain_axs[2], my_adultspans, end_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'End', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Ending Prognosis (Remaining Days)', label_coordinates = (0.5, 1), polyfit_degree = 2)
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
    strain_colors=[[0,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0.5]
    ]
    
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
            graphingFigures.cannedFigures.cohort_scatters(ax_h, my_adultspans/my_adultspans.max(), inflection_data, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Absolute Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Days)', label_coordinates = (12, 2.5), make_labels=make_labels,no_cohorts_color=strain_colors[strain_num],plot_trenddata=False)
        elif anal_mode is 'relative':
            # Plot the traces and scatter for relative inflection.
            graphingFigures.cannedFigures.cohort_scatters(ax_h, my_adultspans/my_adultspans.max(), relative_inflection, strain_health, bin_width_days=animal_bins,bin_mode='percentile', the_title = 'Relative Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Relative Prognosis)', label_coordinates = (4, -0.4),make_labels=make_labels,no_cohorts_color=strain_colors[strain_num],plot_trenddata=False)
        strain_adultspans.append(my_adultspans/my_adultspans.max())
    #ax_h.set_title('r^2:{:3f}\np:{:3f}'.format(*scipy.stats.pearsonr(*strain_adultspans)))
    
    if len(out_dir)>0:
        plotting_tools.clean_plot(ax_h, make_labels=make_labels,suppress_ticklabels=not make_labels)
        dev_fig.savefig(out_dir+os.path.sep+anal_mode+'_deviationrescaling.png')
    return (dev_fig, ax_h)

# Spans analysis
def span_analysis(strain_dfs, strains, out_dir='', make_labels=True, global_cutoff=True,plot_mode='separate'):
    # Use percentile bins
    animal_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_dfs[0], my_worms = strain_dfs[0].worms, bin_width_days = 2)[0]])
    animal_bins = 100*np.cumsum(animal_bins)/sum(animal_bins) # Make percentile bins
    
    if global_cutoff:
        flat_data = np.concatenate(np.array([np.ndarray.flatten(strain_df.mloc(strain_df.worms,['health'])) for strain_df in strain_dfs])).ravel()
        flat_data = flat_data[~np.isnan(flat_data)]
        my_cutoff = np.percentile(flat_data, 0.5*100)
        my_cutoff = strain_dfs[0].display_variables(my_cutoff, 'health')[0]
        
    span_figs, span_axs = [], []
    if plot_mode is 'separate':
        for strain, strain_health in zip(strains,strain_dfs):
            span_fig, ax_h = plt.subplots(2,2)
            graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
                strain_health,make_labels=make_labels,
                cohort_info=analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = animal_bins,bin_mode='percentile'),
                my_cutoff=None if not global_cutoff else my_cutoff,
                bar_plot_mode='uniform')
            span_figs.append(span_fig)
            span_axs.append(ax_h)
            if len(out_dir)>0:
                [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
                span_fig.savefig(out_dir+os.path.sep+strain+'_healthspans.png')
    elif plot_mode is 'comb_collapsemutant':
        span_fig, ax_h = plt.subplots(2,2)
        for strain_num,(strain, strain_health) in enumerate(zip(strains,strain_dfs)):
            if strain_num is 0:
                graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
                    strain_health,make_labels=make_labels,
                    cohort_info=analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = animal_bins,bin_mode='percentile'),
                    my_cutoff=None if not global_cutoff else my_cutoff,
                    bar_plot_mode='uniform')
            else:
                cohort_info = analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = np.array([100]),bin_mode='percentile')
                cohort_info[3] = np.array([[42/255,160/255,75/255]])
                graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
                    strain_health,make_labels=make_labels,
                    cohort_info=cohort_info,
                    my_cutoff=None if not global_cutoff else my_cutoff,
                    bar_plot_mode='uniform',
                    stop_with_death=False)
            span_figs.append(span_fig)
            span_axs.append(ax_h)
        if len(out_dir)>0:
            [plotting_tools.clean_plot(my_ax, make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
            span_fig.savefig(out_dir+os.path.sep+'comb_collapsemutant_healthspans.png')

    return (span_figs,span_axs)
    
if __name__ is "__main__":
    make_labels= False
    #strains = ['spe-9','age-1']
    #out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/processed_data/age-1_cohorts+regression_20160818/'
    #out_dir = ''
    #strains = ['spe-9percombinedweightedSVM','age-1percombinedweightedSVM']
    strains = ['spe-9','age-1specific']
    out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/processed_data/spe-9+age-1percombinedweightedSVM_20161006/'
    out_dir=''
    if make_labels and out_dir is not '': 
        out_dir = out_dir+os.path.sep+'labeled'+os.path.sep
        if not os.path.isdir(out_dir): os.mkdir(out_dir)
    
    strain_dfs = load_strain_data(strains)
    plot_strain_rescaling(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels,do_bootstrap=True)
    plot_strain_rescaling(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels,do_bootstrap=True,percent_range=[20,40])
    plot_strain_rescaling(strain_dfs,strains,out_dir=out_dir,make_labels=make_labels,do_bootstrap=True,percent_range=[60,80])
    bob
    
    plot_survival(strain_dfs,out_dir=out_dir,make_labels=make_labels)
    
    #plot_strain_health(strain_dfs,group_mode='population',make_labels=make_labels,out_dir=out_dir)
    plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir)
    plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir,var_plot_mode='separate')
    plot_strain_health(strain_dfs,group_mode='spe9_allcohorts',make_labels=make_labels,out_dir=out_dir,collapse_mutant=True)
    plot_strain_health(strain_dfs,group_mode='spe9_3cohorts',make_labels=make_labels,out_dir=out_dir)
    plot_strain_health(strain_dfs,group_mode='equalbins_7',make_labels=make_labels,out_dir=out_dir)
    
    parameter_analysis(strain_dfs,
        out_dir=out_dir,make_labels=make_labels)
    
    [abs_fig, abs_ax] = deviation_analysis(strain_dfs,'absolute',
        out_dir=out_dir,make_labels=make_labels)
    [rel_fig, rel_ax] = deviation_analysis(strain_dfs,'relative',
        out_dir=out_dir,make_labels=make_labels)
    
    span_analysis(strain_dfs, strains, out_dir=out_dir,make_labels=make_labels)
    span_analysis(strain_dfs, strains, out_dir=out_dir,make_labels=make_labels,plot_mode='comb_collapsemutant')

    
    [abs_fig, abs_ax] = deviation_rescaling(strain_dfs,'absolute',
        out_dir=out_dir,make_labels=make_labels)
    [abs_fig, abs_ax] = deviation_rescaling(strain_dfs,'relative',
        out_dir=out_dir,make_labels=make_labels)
