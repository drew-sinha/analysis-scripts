#from wzhang import analyzeHealth, graphingFigures
# wzhang files on path
import graphingFigures
import analyzeHealth

import matplotlib.pyplot as plt
import numpy as np
import pickle
import collections

plt.ion()
plt.show()
plt.close('all')

'''
Things to run
1. Survival curves DONE
2. Individual health variables and overall health DONE
3. Absolute and relative healthspans DONE
4. Relative health by deviation DONE
5. Redo analyses with defining cohorts so that percentiles match b/t spe-9 and age-1
6. Do some analysis to address scaling...
'''

out_dir = '/media/Data/Work/ZPLab/Analysis/MutantHealth/age-1_fixedbins/'
reload_data = True
make_labels = False

# Loads dfs
strains = ['spe-9','age-1']
#strain_health = collections.OrderedDict()
if reload_data:
    strain_dfs = []
    for strain in strains:
        with open('/mnt/bulkdata/wzhang/human_dir/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
            strain_dfs.append(pickle.load(my_file)['adult_df'])
        print(strain+": n = {}".format(len(strain_dfs[-1].worms)))

# Plot survival curves and lifespan distributions for each strain
survival_fig, ax_h = plt.subplots(len(strains),2)
for strain_health, strain_ax in zip(strain_dfs,ax_h.T):
    #survival_plot = plotFigures.consistent_subgrid_coordinates((6, 8), (1, 6), my_width = 2, my_height = 2)
    #lifespans_plot = plotFigures.consistent_subgrid_coordinates((6, 8), (3, 6), my_width = 2, my_height = 2)
    graphingFigures.cannedFigures.survival_lifespan(strain_ax[0],strain_ax[1],strain_health, make_labels=make_labels)
    
if out_dir is not '': survival_fig.savefig(out_dir+'survival_curves.svg')
if out_dir is not '': survival_fig.savefig(out_dir+'survival_curves.png')


# Plot health variables
health_vars = ['bulk_movement', 'intensity_80', 'life_texture', 'cumulative_eggs','adjusted_size']
health_fig, ax_h = plt.subplots(len(health_vars),len(strains))
for strain_health, strain_plots in zip(strain_dfs, ax_h.T):
    for var_plot, var in zip(strain_plots, health_vars):
        graphingFigures.cannedFigures.cohort_traces(var_plot,var, strain_health, make_labels=make_labels)
if out_dir is not '': health_fig.savefig(out_dir+'health_vars.svg')
if out_dir is not '': health_fig.savefig(out_dir+'health_vars.png')

# Plot overall prognosis
health_fig, ax_h = plt.subplots(1,len(strains))
for strain_health, health_plot in zip(strain_dfs, ax_h):
    graphingFigures.cannedFigures.cohort_traces(health_plot,'health', strain_health, make_labels=make_labels)
if out_dir is not '': health_fig.savefig(out_dir+'health_overall.svg')
if out_dir is not '': health_fig.savefig(out_dir+'health_overall.png')


# Deviation analysis
def deviation_analysis(strain_dfs, mode='absolute', make_labels=True):
    dev_fig, ax_h = plt.subplots(2, len(strain_dfs))
    
    for strain_num, (strain_health, strain_axs) in enumerate(zip(strain_dfs, ax_h.T)):
        # Make bins of lifespans.
        (life_cohorts, bin_lifes, my_bins, my_colors) = analyzeHealth.selectData.adult_cohort_bins(strain_health, my_worms = strain_health.worms, bin_width_days = 2)
        my_adultspans = analyzeHealth.selectData.get_adultspans(strain_health)/24  
        my_cohorts = life_cohorts

        # Prepare my "inflection" data. 
        geometry_dict = analyzeHealth.computeStatistics.one_d_geometries(strain_health, 'health')
        start_data = geometry_dict['start']
        mean_start = np.mean(start_data)            
        inflection_data = geometry_dict['absolute_inflection']
        relative_inflection = geometry_dict['self_inflection']
        
        if mode is 'absolute':
            # Plot the traces and scatter for absolute inflection.
            graphingFigures.cannedFigures.cohort_traces(strain_axs[0], 'health', strain_health, the_title = 'Prognosis Over Normalized Time', the_xlabel = 'Fractional Adult Lifespan', the_ylabel = 'Prognosis (Remaining Days)', x_normed = True, make_labels=make_labels)
            strain_axs[0].set_ylim([0, 1.1*mean_start])
            graphingFigures.cannedFigures.cohort_scatters(strain_axs[1], my_adultspans, inflection_data, strain_health, the_title = 'Absolute Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Days)', label_coordinates = (12, 2.5), make_labels=make_labels)
        elif mode is 'relative':
            # Plot the traces and scatter for relative inflection.
            graphingFigures.cannedFigures.cohort_traces(strain_axs[0], 'health', strain_health, the_title = 'Relative Prognosis Over Normalized Time', the_xlabel = 'Fractional Adult Lifespan', the_ylabel = 'Relative Prognosis (Fractional Remaining Life)', x_normed = True, y_normed = True, make_labels=make_labels)
            strain_axs[0].set_ylim([-0.1, 1.1])
            graphingFigures.cannedFigures.cohort_scatters(strain_axs[1], my_adultspans, relative_inflection, strain_health, the_title = 'Relative Deviation', the_xlabel = 'Days of Adult Lifespan', the_ylabel = 'Average Deviation (Relative Prognosis)', label_coordinates = (4, -0.4),make_labels=make_labels)
    return (dev_fig, ax_h)

[abs_fig, abs_ax] = deviation_analysis(strain_dfs,'absolute', make_labels=make_labels)
[rel_fig, rel_ax] = deviation_analysis(strain_dfs,'relative', make_labels=make_labels)
if out_dir is not '': abs_fig.savefig(out_dir+'abs_deviation.svg')
if out_dir is not '': abs_fig.savefig(out_dir+'abs_deviation.png')
if out_dir is not '': rel_fig.savefig(out_dir+'rel_deviation.svg')
if out_dir is not '': rel_fig.savefig(out_dir+'rel_deviation.png')

# Spans analysis
for strain, strain_health in zip(strains,strain_dfs):
    span_fig, ax_h = plt.subplots(2,2)
    graphingFigures.cannedFigures.show_spans(ax_h[0,0],ax_h[0,1],ax_h[1,0],ax_h[1,1],
        strain_health,make_labels=make_labels)
    if out_dir is not '': span_fig.savefig(out_dir+strain+'healthspans.svg')
    if out_dir is not '': span_fig.savefig(out_dir+strain+'healthspans.png')

# Remake Willie's original health figure
health_fig, ax_h = plt.subplots(1,1)
graphingFigures.cannedFigures.cohort_traces(ax_h,'health', strain_dfs[0], make_labels=True,the_xlabel='Days Post-Maturity')
if out_dir is not '': health_fig.savefig(out_dir+'Willie_health_overall.svg')
if out_dir is not '': health_fig.savefig(out_dir+'Willie_health_overall.png')
