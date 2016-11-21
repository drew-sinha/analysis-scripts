import analyzeHealth
import matplotlib.pyplot as plt
import pickle
import numpy as np
import collections
import os
import numpy as np
import scipy.stats

import plotting_tools
import zplib.scalar_stats.kde
import zplib.image.colorize

do_dynamics=False
do_residuals=True
do_demofigs=False

make_labels = False
out_dir = '/media/Data/Work/ZPLab/Analysis/ResidualAnalysis/'
#out_dir = ''
#out_dir=''
if make_labels and len(out_dir)>0:
    out_dir = out_dir + '/labeled/'
    if not os.path.isdir(out_dir): os.mkdir(out_dir)

plt.ion()
plt.close('all')
plt.show()

#def get_cohort_data(adult_df, cohort_assignments,variable_to_get='health', stop_with_death=True, skip_conversion=False):
    #'''
        #cohort_assignments - Output from adult_cohort_bins (lists of indices for worms in adult_df corresponding to each cohort)
        #variable_to_get - String corresponding to variable of interest
        #normalize_var - Normalize variables into the range [0,1]
    #'''
    
    #cohort_data = []
    #cohort_ages = []
    #for a_cohort in cohort_assignments:
        #worm_cohort_data = adult_df.mloc(adult_df.worms, [variable_to_get])[a_cohort, 0,:] # Data for individual worms in cohort
        #worm_cohort_data = worm_cohort_data[~np.isnan(worm_cohort_data).all(axis=1)]
        #if stop_with_death:
            #cohort_data.append(np.mean(worm_cohort_data,axis=0))
        #else:
            #cohort_data.append(np.nanmean(worm_cohort_data,axis=0))
        #if not skip_conversion:
            #(cohort_data[-1], my_unit, fancy_name) = adult_df.display_variables(cohort_data[-1], variable_to_get)
        #cohort_data[-1] = cohort_data[-1][~np.isnan(cohort_data[-1])]
        ##if normalize_var: cohort_data[-1] = (cohort_data[-1]-cohort_data[-1][0])/(cohort_data[-1][-1]-cohort_data[-1][0])
        #cohort_ages.append(adult_df.ages[:cohort_data[-1].shape[0]])
    #return (cohort_data, cohort_ages)

# Load data
strain='spe-9'
with open('/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/'+strain+'_health/df_'+strain+'.pickle','rb') as my_file:
    strain_df=pickle.load(my_file)['adult_df']
    print(strain+": n = {}".format(len(strain_df.worms)))
    print(strain+": CV={:.2f}".format(
        np.std(analyzeHealth.selectData.get_adultspans(strain_df))/np.mean(analyzeHealth.selectData.get_adultspans(strain_df))))

if do_dynamics:
    # Dynamical modeling
    # Start with mean health of population
    animal_bins = np.array([100])
    adult_cohort_bin_data = analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days=animal_bins,bin_mode='percentile')
    cohort_assignments = adult_cohort_bin_data[0]
    my_cohort_data = analyzeHealth.selectData.get_cohort_data(strain_df, cohort_assignments, stop_with_death=False)

    cohort_data_rate = np.diff(my_cohort_data[0][0])
    cohort_data_rate = np.append(cohort_data_rate[0],cohort_data_rate)
    cohort_data_doublerate = np.diff(cohort_data_rate)
    cohort_data_doublerate = np.append(cohort_data_doublerate[0], cohort_data_doublerate)

    # Plot data as time series
    fig_h, ax_h = plt.subplots(3,1)
    fig_h.suptitle('Mean Health of Population - Aging Across Time')
    [my_ax.plot(my_cohort_data[1][0],my_data) for my_data,my_ax in zip([my_cohort_data[0][0],cohort_data_rate,cohort_data_doublerate],ax_h)]
    [my_ax.set_xlabel('Time (d)') for my_ax in ax_h]
    [my_ax.set_ylabel(my_label) for my_label,my_ax in zip(['Health','d/dt(Health)','d2/dt2(Health)'],ax_h)]

    fig_h, ax_h = plt.subplots(1,3)
    fig_h.suptitle('Mean Health of Population - Phase Space Analysis')
    ax_h[0].plot(my_cohort_data[0][0],cohort_data_rate)
    ax_h[0].set_xlabel('Health')
    ax_h[0].set_ylabel('d/dt(Health)')
    ax_h[1].plot(my_cohort_data[0][0], cohort_data_doublerate)
    ax_h[1].set_xlabel('Health')
    ax_h[1].set_ylabel('d2/dt2(Health)')
    ax_h[2].plot(cohort_data_rate, cohort_data_doublerate)
    ax_h[2].set_xlabel('d/dt(Health)')
    ax_h[2].set_ylabel('d2/dt2(Health')


    # Individual trajectories
    max_age = 14

    indiv_health = strain_df.display_variables(strain_df.mloc(measures=['health']), 'health')[0][:,0,:]   # Get rid of 'measures' axis
    indiv_ages = strain_df.ages.copy()
    if max_age is not None:
        last_time_idx = np.where(np.array(strain_df.ages)<=max_age)[0][-1]+1
        indiv_health = indiv_health[:,:last_time_idx]
        indiv_ages = indiv_ages[:last_time_idx]

    indiv_health_vel = np.diff(indiv_health, axis=1)
    #indiv_health_vel = np.append(indiv_health_vel[:,0][:,np.newaxis], indiv_health_vel, axis=1)
    indiv_health_acc = np.diff(indiv_health,n=2, axis=1)
    #indiv_health_acc = np.append(
        #np.repeat(indiv_health_acc[:,0][:,np.newaxis],repeats=2,axis=1), 
        #indiv_health_acc, axis=1)
    indiv_health = indiv_health[:,2:]
    indiv_health_vel = indiv_health_vel[:,1:]
    indiv_ages = indiv_ages[2:]

    # Plot time series
    fig_h, ax_h = plt.subplots(3,1)
    fig_h.suptitle('Indiv Trajectories Across Time')
    [my_ax.plot(indiv_ages, worm_health_var) for worm_health_var, my_ax in zip([indiv_health.T,indiv_health_vel.T,indiv_health_acc.T], ax_h)]
    [my_ax.set_xlabel('Time (d)') for my_ax in ax_h]
    [my_ax.set_ylabel(my_label) for my_label,my_ax in zip(['Health','d/dt(Health)','d2/dt2(Health)'],ax_h)]

    # Phase-space analysis
    fig_h, ax_h = plt.subplots(1,3)
    fig_h.suptitle('Indiv Trajectories - Phase Space Analysis')
    ax_h[0].plot(indiv_health.T,indiv_health_vel.T)
    ax_h[0].set_xlabel('Health')
    ax_h[0].set_ylabel('d/dt(Health)')
    ax_h[1].plot(indiv_health.T, indiv_health_acc.T)
    ax_h[1].set_xlabel('Health')
    ax_h[1].set_ylabel('d2/dt2(Health)')
    ax_h[2].plot(indiv_health_vel.T,indiv_health_acc.T)
    ax_h[2].set_xlabel('d/dt(Health)')
    ax_h[2].set_ylabel('d2/dt2(Health)')


    # Cohorts
    spe9_bins = np.array([len(my_bin) for my_bin in analyzeHealth.selectData.adult_cohort_bins(strain_df, my_worms = strain_df.worms, bin_width_days = 2)[0]])
    spe9_bins = 100*np.cumsum(spe9_bins)/sum(spe9_bins) # Make percentile bins
    adult_cohort_bin_data = analyzeHealth.selectData.adult_cohort_bins(strain_df, bin_width_days=spe9_bins,bin_mode='percentile')
    compiled_cohort_data = analyzeHealth.selectData.get_cohort_data(strain_df, adult_cohort_bin_data[0])

    #def derive_dynamics(adult_df, worms=None, max_age=14,trim_front=True):
        #health_data = adult_df.display_variables(strain_df.mloc(worms=worms,measures=['health']), 'health')[0][:,0,:]
        #ages = adult_df.ages.copy()
        
        #if max_age is not None:
            #last_time_idx = np.where(np.array(ages)<=max_age)[0][-1]+1
            #health_data = health_data[:,:last_time_idx]
            #ages = ages[:last_time_idx]
        
        #health_vel = np.diff(health_data,axis=1)
        #health_acc = np.diff(health_data,n=2,axis=1)
        
        #if trim_front:
            #health_data = health_data[:,2:]
            #ages = ages[2:]
            #health_vel = health_vel[:,1:]
        
        #return [health_data, ages, health_vel, health_acc]
        
    def derive_dynamics(health_data, ages, max_age=14,trim_front=True): 
        if max_age is not None:
            last_time_idx = np.where(np.array(ages)<=max_age)[0][-1]+1
            health_data = health_data[:last_time_idx]
            ages = ages[:last_time_idx]
        
        health_vel = np.diff(health_data)   # Uses last axis! should be good.
        health_acc = np.diff(health_data,n=2)
        
        if trim_front:
            health_data = health_data[2:]
            ages = ages[2:]
            health_vel = health_vel[1:]
        
        return [health_data, ages, health_vel, health_acc]

    def plot_dynamics_timeseries(health_data, ages, health_vel, health_acc, fig_h=None):
        if fig_h is None: fig_h, ax_h = plt.subplots(1,3)
        else: ax_h = fig_h.get_axes()
        
        [my_ax.plot(ages, worm_health_var) for worm_health_var, my_ax in zip([health_data.T,health_vel.T,health_acc.T], ax_h)]
        [my_ax.set_xlabel('Time (d)') for my_ax in ax_h]
        [my_ax.set_ylabel(my_label) for my_label,my_ax in zip(['Health','d/dt(Health)','d2/dt2(Health)'],ax_h)]
        
        return [fig_h,ax_h]
        
    def plot_phaseplane(health_data, ages, health_vel, health_acc, fig_h= None):
        if fig_h is None: fig_h, ax_h = plt.subplots(1,3)
        else: ax_h = fig_h.get_axes()

        ax_h[0].plot(health_data.T,health_vel.T)
        ax_h[0].set_xlabel('Health')
        ax_h[0].set_ylabel('d/dt(Health)')
        ax_h[1].plot(health_data.T, health_acc.T)
        ax_h[1].set_xlabel('Health')
        ax_h[1].set_ylabel('d2/dt2(Health)')
        ax_h[2].plot(health_vel.T,health_acc.T)
        ax_h[2].set_xlabel('d/dt(Health)')
        ax_h[2].set_ylabel('d2/dt2(Health)')
        
        return [fig_h, ax_h]


    compiled_traj = []
    #for cohort_assignments in adult_cohort_bin_data[0]:
        #compiled_traj.append(derive_trajectories(strain_df, worms=np.array(strain_df.worms)[cohort_assignments]))

    for [health_data, ages] in zip(*compiled_cohort_data):
        compiled_traj.append(derive_cohort_dynamics(health_data, ages))

    # Plot by cohort
    # Get colors
    #my_colors = zplib_image_colorize.color_map(bin_lifes/np.max(bin_lifes[~np.isnan(bin_lifes)]))
    #my_colors = my_colors/255

    fig_h, ax_h = plt.subplots(3,1)
    [plot_dynamics_timeseries(health_data,ages, health_vel, health_acc, fig_h) for [health_data,ages,health_vel,health_acc] in compiled_traj]

    fig_h, ax_h = plt.subplots(3,1)
    [plot_phaseplane(health_data,ages, health_vel, health_acc, fig_h) for [health_data,ages,health_vel,health_acc] in compiled_traj]

if do_residuals:
    def calc_diff_residuals(health_data, ages, bin_mode='time',bin_interval=2,flatten=True,return_stats=False,normalize_mode='zscore'):
        inc_healthchg = np.diff(health_data,axis=1)
        mean_traj = np.nanmean(inc_healthchg,axis=0)
        var_traj = np.nanvar(inc_healthchg,axis=0)
        
        resid=[]
        for epoch in range(int(np.ceil(ages[-1]/bin_interval))):
            epoch_timepts = (ages[1:]>=epoch*bin_interval) & (ages[1:]<(epoch+1)*bin_interval) # Do 1: since lose one point from the differenching
            temp_resid = inc_healthchg[:,epoch_timepts]
            if normalize_mode=='zscore':
                temp_resid = (temp_resid-mean_traj[np.newaxis,epoch_timepts])/(np.sqrt(var_traj))[np.newaxis,epoch_timepts]
            elif normalize_mode is not 'none':
                raise Exception('normalize_mode not implemented')
            if flatten:
                resid.append(temp_resid[~np.isnan(temp_resid)].flatten())
            else:
                resid.append(temp_resid)

        if not return_stats:
            return resid
        else:
            return (resid, mean_traj,var_traj)
    def plot_residuals(resid,dist_mode='both',plot_gaussian=True,display_info=True,**kws):
        '''
            resid - M x N list containing n arrays of residuals
            dist_mode - Mode for displaying the distribution
            plot_gaussian - Flag for superimposing the corresponding Gaussian (i.e. of same mean and variance)
        '''
        if 'x_min' not in kws:
            x_min = -5
        else: x_min = kws['x_min']
        if 'x_max' not in kws:
            x_max = 5
        else: x_max = kws['x_max']
        if 'draw_color' not in kws:
            draw_color='b'
        else: draw_color = kws['draw_color']
        
        if 'ax_h' in kws:
            ax_h = kws['ax_h']
            fig_h = ax_h[0][0].figure
        else:
            resid_shape = [len(resid), len(resid[0])]
            fig_h, ax_h = plt.subplots(*resid_shape,squeeze=False)
        if 'orientation' in kws:
            orientation = kws['orientation']
        else: orientation='vertical'
        
        for interval_resid, interval_ax in zip(plotting_tools.flatten_list(resid,to_level=1),ax_h.flatten()):
            #print(interval_resid)
            if interval_resid.size is not 0:
                if dist_mode is 'bar' or dist_mode is 'both':
                    interval_ax.hist(interval_resid,bins=20,range=[interval_resid.min(),interval_resid.max()],normed=True,color='w',edgecolor=draw_color,orientation=orientation)
                if dist_mode is 'kde' or dist_mode is 'both':
                    interval_supp,interval_density,interval_obj = zplib.scalar_stats.kde.kd_distribution(interval_resid)
                    if orientation is 'vertical': interval_ax.plot(interval_supp,interval_density,color=draw_color,linewidth=2)
                    else: interval_ax.plot(interval_density,interval_supp,color=draw_color,linewidth=2)
                if plot_gaussian:
                    [data_mean,data_var] = [np.mean(interval_resid), np.var(interval_resid)]
                    #print(data_mean, data_var)
                    gauss_supp = np.linspace(interval_resid.min(),interval_resid.max(),20 if dist_mode is 'bar' else len(interval_supp))
                    if orientation is 'vertical':
                        interval_ax.plot(gauss_supp, 1/np.sqrt(2*np.pi*data_var)*np.exp(-(gauss_supp-data_mean)**2/(2*data_var)),color='g',linewidth=2)
                    else:
                        interval_ax.plot(1/np.sqrt(2*np.pi*data_var)*np.exp(-(gauss_supp-data_mean)**2/(2*data_var)),gauss_supp,color='g',linewidth=2)
                    
                if display_info: interval_ax.text(0.05,0.95,
                    'n={}\nD={:.3f}\np={:.3f}'.format(interval_resid.size,*(scipy.stats.kstest(interval_resid,'norm'))),
                    fontsize=10,transform=interval_ax.transAxes)
            #interval_ax.set_ylim([0,1])
            interval_ax.set_xlim([x_min,x_max])
        return fig_h,ax_h
    
    # Individual residuals binned over diff. times
    indiv_health = strain_df.display_variables(strain_df.mloc(measures=['health']), 'health')[0][:,0,:]   # Get rid of 'measures' axis
    indiv_ages = np.array(strain_df.ages.copy())
    (residuals, indiv_mean_traj, indiv_var_traj) = calc_diff_residuals(indiv_health,indiv_ages, return_stats=True)
    
    # Plot mean and variance of the incremental residuals over all intervals
    fig_h, ax_h = plt.subplots(2,1)
    ax_h[0].plot(indiv_ages[1:], indiv_mean_traj)
    ax_h[0].set_xlabel('Days post-adulthood')
    ax_h[0].set_ylabel('Mean incremental residual of prognosis (d)')
    ax_h[1].plot(indiv_ages[1:],indiv_var_traj)
    ax_h[1].set_xlabel('Days post-adulthood')
    ax_h[1].set_ylabel('Variance in incremental residual of prognosis (d^2)')
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h]
        fig_h.savefig(out_dir+'meanvar_pooled_acrosstime.png')
    
    # Plot first interval pooled
    fig_h, ax_h = plot_residuals([[residuals[0]]])
    ax_h[0,0].set_title('[0,2]d pooled')
    fig_h.canvas.set_window_title('First Interval (Pooled)')
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'firstinterval_pooled.png')    
    
    # Plot all intervals (with pooled subintervals) superimposed on each other
    fig_h,ax_h=plt.subplots(1,1,squeeze=False)
    [plot_residuals([[interval_resid]],dist_mode='kde',ax_h=ax_h,display_info=False,draw_color=my_color) for interval_resid,my_color in zip(residuals,plt.cm.Reds(np.linspace(0.1,0.9,len(residuals))))]
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'allintervals_pooled_superimposed.png')
    
    # Plot all intervals with each interval pooling all subintervals' residuals
    fig_h, ax_h = plot_residuals(np.atleast_2d(residuals))
    [my_ax.set_title('[{},{}]d'.format(plot_num*2,(plot_num+1)*2)) for plot_num,my_ax in enumerate(ax_h.flatten())]
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'allintervals_pooled.png')
    
    # By cohort
    adult_cohort_bin_data = analyzeHealth.selectData.adult_cohort_bins(strain_df, my_worms = strain_df.worms, bin_width_days = 2)
    compiled_resid = []
    compiled_mean_traj = []
    compiled_var_traj = []
    for cohort_bins in adult_cohort_bin_data[0]:
        cohort_data = calc_diff_residuals(indiv_health[cohort_bins,:],np.array(indiv_ages),return_stats = True)
        compiled_resid.append(cohort_data[0])
        compiled_mean_traj.append(cohort_data[1])
        compiled_var_traj.append(cohort_data[2])
    
    # Plot mean and variance of cohort health decline over time
    fig_h,ax_h = plt.subplots(2,1)
    [ax_h[0].plot(indiv_ages[1:], cohort_mean_traj,color=cohort_color) for cohort_mean_traj,cohort_color in zip(compiled_mean_traj,adult_cohort_bin_data[3])]
    ax_h[0].set_xlabel('Days post-adulthood')
    ax_h[0].set_ylabel('Mean incremental residual of prognosis (d)')
    [ax_h[1].plot(indiv_ages[1:], cohort_var_traj,color=cohort_color) for cohort_var_traj,cohort_color in zip(compiled_var_traj,adult_cohort_bin_data[3])]
    ax_h[1].set_xlabel('Days post-adulthood')
    ax_h[1].set_ylabel('Variance in incremental residual of prognosis (d^2)')
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h]
        fig_h.savefig(out_dir+'meanvar_stratified_acrosstime.png')
    
    # Plot first interval by cohort
    fig_h,ax_h=plot_residuals(np.atleast_2d([cohort_resid[0] for cohort_resid in compiled_resid]))
    fig_h.canvas.set_window_title('First Interval (Stratified)')
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'firstinterval_cohortstratified.png')
    
    # Plot first interval with cohorts superimposed on same figure
    first_residuals=[cohort_resid[0] for cohort_resid in compiled_resid]
    fig_h,ax_h=plt.subplots(1,1,squeeze=False)
    [plot_residuals([[cohort_resid]],dist_mode='kde',ax_h=ax_h,display_info=False,draw_color=cohort_color) for cohort_resid,cohort_color in zip(first_residuals,adult_cohort_bin_data[3])]
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'firstinterval_cohorts_superimposed.png')
    
    # Plot all intervals by cohort
    fig_h, ax_h = plot_residuals(np.array(compiled_resid))
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'allintervals_cohortstratified.png')
    
    # Correlations in individual residuals in absolute time
    bin_interval=2
    ordered_residuals = calc_diff_residuals(indiv_health,indiv_ages,bin_interval=bin_interval,flatten=False,normalize_mode='none')
    
    # First pool all data together
    
    start_resid, end_resid = np.array([]),np.array([])
    for start_interval,end_interval in zip(ordered_residuals[:-1],ordered_residuals[1:]):
        alive_individuals = (~np.isnan(start_interval).any(axis=1)) & (~np.isnan(end_interval).any(axis=1))
        start_resid = np.append(start_resid, np.sum(start_interval,axis=1)[alive_individuals])
        end_resid = np.append(end_resid, np.sum(end_interval,axis=1)[alive_individuals])
    corr_data = scipy.stats.pearsonr(start_resid,end_resid)   
    fig_h, ax_h = plt.subplots(1,1)
    ax_h.scatter(start_resid, end_resid)
    ax_h.set_title('r^2={:.3f}\np = {:.3f}'.format(
        corr_data[0]**2,corr_data[1]))
    if len(out_dir) > 0:
        plotting_tools.clean_plot(ax_h,make_labels=make_labels,suppress_ticklabels=not make_labels)
        fig_h.savefig(out_dir+'intervalresidual_correlations_pooled.png')
    
    # Then, do over time
    fig_h, ax_h = plt.subplots(1,len(ordered_residuals)-1)
    for ax_num,(my_ax,start_interval,end_interval) in enumerate(zip(ax_h,ordered_residuals[:-1],ordered_residuals[1:])):            
        alive_individuals = (~np.isnan(start_interval).any(axis=1)) & (~np.isnan(end_interval).any(axis=1))
        total_chg = np.array([np.sum(start_interval,axis=1),np.sum(end_interval,axis=1)]).T
        my_ax.scatter(total_chg[alive_individuals,0],total_chg[alive_individuals,1])
        corr_data = scipy.stats.pearsonr(total_chg[alive_individuals,0],total_chg[alive_individuals,1])
        print('r={:.3f}, p = {:.3f}'.format(corr_data[0],corr_data[1]))
        my_ax.set_title('r^2={:.3f}\np = {:.3f}'.format(
            corr_data[0]**2,corr_data[1]))
    if len(out_dir) > 0:
        [plotting_tools.clean_plot(my_ax,make_labels=make_labels,suppress_ticklabels=not make_labels) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'intervalresidual_correlations_overtime.png')
        
        
    # Plot MSD
    MSD_health = np.nanmean((indiv_health-indiv_health[:,0][:,np.newaxis])**2,axis=0)
    fig_h,ax_h = plt.subplots(1,1)
    ax_h.plot(indiv_ages,MSD_health)
    ax_h.set_xlabel('Time (d)')
    ax_h.set_ylabel('MSD Prognosis')
    #ax_h.plot(indiv_ages,indiv_ages,'k')
    if len(out_dir) > 0:
        plotting_tools.clean_plot(ax_h,make_labels=make_labels,suppress_ticklabels=not make_labels)
        fig_h.savefig(out_dir+'MSD.png')
        
    MSD_health = np.nanmean((indiv_health-indiv_health[:,0][:,np.newaxis])**2,axis=0)-np.nanmean((indiv_health-indiv_health[:,0][:,np.newaxis]),axis=0)**2
    fig_h,ax_h = plt.subplots(1,1)
    ax_h.plot(indiv_ages,MSD_health)
    ax_h.set_xlabel('Time (d)')
    ax_h.set_ylabel('MSD Prognosis')
    #ax_h.plot(indiv_ages,indiv_ages,'k')
    if len(out_dir) > 0:
        plotting_tools.clean_plot(ax_h,make_labels=make_labels,suppress_ticklabels=not make_labels)
        fig_h.savefig(out_dir+'MSD_detrend.png')
    
    
    # Do everything in relative time?!?!?! (rel. to lifespan?)

if do_demofigs:
    # Illustrations
    norm_sample1 = np.random.normal(size=(1,1000))
    norm_sample2 = np.random.normal(size=(1,1000))
    norm_sample3_cumsum = np.random.normal(size=(2,1000)).cumsum(axis=1)
    
    # Plot realizations of 1D Brownian motion
    fig_h,ax_h = plt.subplots(1,1)
    [ax_h.plot(np.random.normal(size=(30,1)),color=[0.5,0.5,0.5]) for i in range(10)]
    plotting_tools.clean_plot(ax_h,make_labels=make_labels,suppress_ticklabels=not make_labels)
    if len(out_dir)>0:
        fig_h.savefig(out_dir+'/demo_figs/BrownianMotion_1D.png')
    
    # Plot realization of 2D brownian motion
    #fig_h, ax_h = plt.subplots(1,1)
    #for start_pt, end_pt in zip(norm_sample3_cumsum[:10,:-1].T,norm_sample3_cumsum[:10,1:].T):
        #ax_h.plot([start_pt[0], end_pt[0]],[start_pt[1],end_pt[1]],'k')
    #ax_h.scatter(0,0,color='g',s=80)
    #ax_h.scatter(norm_sample3_cumsum[0,-1],norm_sample3_cumsum[1,-1],color='r',s=80)
    #plotting_tools.clean_plot(ax_h,make_labels=False,suppress_ticklabels=True)
    #if len(out_dir)>0: fig_h.savefig(out_dir+'/demo_figs/BrownianMotion_2D.png')
    
    # Plot normal distribution
    pts = np.linspace(-2,2,num=201)
    norm_curve = 1/(2*np.pi)*np.exp(-(pts**2)/2)
    fig_h, ax_h = plt.subplots(1,1)
    ax_h.plot(pts,norm_curve)
    plotting_tools.clean_plot(ax_h,make_labels=False,suppress_ticklabels=True)
    if len(out_dir)>0: fig_h.savefig(out_dir+'/demo_figs/normal_dist.png')
    
    # Plot example scatter of two independent sample paths
    fig_h, ax_h = plt.subplots(1,1)
    ax_h.scatter(norm_sample1,norm_sample2)
    plotting_tools.clean_plot(ax_h,make_labels=False,suppress_ticklabels=True)
    if len(out_dir)>0: fig_h.savefig(out_dir+'/demo_figs/ind_scatter.png')
    
    # Plot actual data as an illustration
    # First, plot health for a set of individuals along with mean over time
    num_individuals = 30
    plot_anexample = True
    group_health = strain_df.display_variables(strain_df.mloc(worms=strain_df.worms[1:num_individuals],measures=['health']), 'health')[0][:,0,:]
    mean_traj = np.nanmean(group_health,axis=0)
    mc_traj = group_health-mean_traj
    
    fig_h,ax_h=plt.subplots(2,1, gridspec_kw = {'height_ratios':[3,1]})
    ax_h[0].plot(strain_df.ages,
        group_health.T,color=[0.5,0.5,0.5])
    if plot_anexample:
        ax_h[0].plot(strain_df.ages, group_health[10,:],'g',linewidth=2)
    ax_h[1].plot(strain_df.ages,np.nanmean(group_health,axis=0))
    if len(out_dir)>0:
        [plotting_tools.clean_plot(my_ax,make_labels=False, suppress_ticklabels=True) for my_ax in ax_h.flatten()]
        fig_h.savefig(out_dir+'/demo_figs/health_trajectories.png')
    
    # Next, get mean-centered trajectories
    
    fig_h,ax_h = plt.subplots(1,1)
    ax_h.plot(strain_df.ages, mc_traj.T,color=[0.5,0.5,0.5])
    if plot_anexample:
        ax_h.plot(strain_df.ages, mc_traj[10,:],'g',linewidth=2)
    if len(out_dir)>0:
        plotting_tools.clean_plot(ax_h,make_labels=False, suppress_ticklabels=True)
        fig_h.savefig(out_dir+'/demo_figs/mc_trajectories.png')
        
    # Plot residuals with horizontal distribution on right
    #time_health = indiv_health[:,24][~np.isnan(indiv_health[:,24])]
    #fig_h,ax_h = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3,1]})
    #ax_h[0].plot(strain_df.ages,
        #group_health.T,color=[0.5,0.5,0.5])
    #plot_residuals([[time_health]],plot_gaussian=False,display_info=False,x_min=0,x_max=1,ax_h=np.array([[ax_h[1]]]),orientation='horizontal')
    #ax_h[1].text(0.95*ax_h[1].get_xlim()[1],np.mean(time_health),'*')
    
