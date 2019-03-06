import numpy as np

from zplib.scalar_stats import kde

import survival_plotting

def get_featurespan(worms, feature, cutoff, dwell_time=2*24, return_crossings=False):
    '''Get amount of time feature is remains above/below a prescribed cutoff value;
        assumes that all ages are in *hours*.

        Parameters
            worms - Worms object of worm data
            feature - str denoting feature of interest in the worms object
            cutoff - float value for cutoff (in units of the feature)
            dwell_time - (optional) minimum time an individual must stay below the
                cutoff to continue in "featurespan" (in hours)
            return_crossings - optional bool denoting whether to return the feature value
                at the time which an individual passed into poor health (good for
                debugging and visualization)

        Returns
            numpy array of featurespans (in hours)
            (optional, if return_crossings is True) crossing values

    '''

    # Ensure that feature decreases with time
    # TODO: Make the features to negate a little less hardcoded?
    if any([label in feature for label in ['af', 'intensity', 'autofluorescence']]):
        negate = True
    else: negate = False

    featurespans, crossing_vals = np.array([]), np.array([])
    for worm in worms:
        ages, feature_data = worm.get_time_range(feature, min_age=0, age_feature='adult_age')
        if negate:
            feature_data *= -1

        adjusted_data = feature_data - cutoff
        crossings = np.where((adjusted_data[:-1] > 0) & (adjusted_data[1:] < 0))[0]

        if len(crossings) == 0:
            if adjusted_data[0] > 0:
                featurespan = worm.adultspan
                crossing_val = feature_data[-1]
            elif adjusted_data[0] <= 0:
                featurespan = 0
                crossing_val = feature_data[0]
        else:
            found_crossing = False
            for crossing_index in crossings:
                if ages[crossing_index] + dwell_time > ages[-1]: # Case where crossing is near death
                    dwell_end = len(ages)
                else:
                    dwell_end = np.where(ages > (ages[crossing_index] + dwell_time))[0][0]

                if (adjusted_data[crossing_index+1:dwell_end] < 0).all():
                    featurespan = ages[crossing_index]
                    crossing_val = feature_data[crossing_index]
                    found_crossing = True
                    break
            if not found_crossing:
                featurespan = ages[crossing_index] # Default to the last crossing in a bad situation
                crossing_val = feature_data[crossing_index]

        featurespans = np.append(featurespans, featurespan)
        crossing_vals = np.append(crossing_vals, crossing_val)

    if return_crossings:
        return featurespans, crossing_vals
    else:
        return featurespans

def plot_expt_hs(worms, feature, cutoff, dwell_time=2*24, ax_h=None,mode='cdf',**plot_kws):
    """Plot survival curves of falling into poor health for one or more separate experiments

        Parameters
            worms
            feature
            cutoff
            dwell_time
            expt_dirs - one or more experiment root directories with collated
                lifespan data (i.e. using the elegant pipeline to produce
                timecourse files)
            ax_h - optional matplotlib axis objects to plot curves on
            import_mode - bool flag that toggles how to import lifespan data;
                if True, searches for a timecourse file to use; otherwise, read
                directly from annotation files
            calc_adultspan - bool flag that toggles whether to calculate lifespan as adultspan;
                if True, uses the 'adult' timepoint as the starting timepoint of interest;
                otherwise, uses the 'larva' timepoint
            plot_kws - optional kw parameters to pass to plt.plot

        Returns
            metadata - list of dicts where each entry corresponds to
                derived metadata about the plotted experiment; this
                includes:
                    n - number of animals
                    mean - mean lifespan of animals
                    median - median lifespan
                    good_ls - numpy array of good lifespans
            (fig_h) - matplotlib figure object if supplied ax_h is None
            (ax_h) - matplotlib axis object if supplied ax_h is None
    """

    if ax_h is None:
        fig_h, ax_h = plt.subplots()
        ax_provided = False
    else: ax_provided = True

    healthspans = get_featurespan(worms,feature, cutoff,dwell_time=dwell_time)/24
    metadata = {'n':len(healthspans),
        'mean': np.nanmean(healthspans),
        'median': np.nanmedian(healthspans),
        'hs': healthspans,}

    if mode == 'cdf':
        survival_plotting.plot_spanseries(healthspans, ax_h=ax_h, **plot_kws)
        ax_h.set_ylabel('Proportion Remaining in Good Health')
    elif mode in ['density']:
        support, density, estimator = kde.kd_distribution(healthspans)
        ax_h.plot(support, density)
        ax_h.set_ylabel('Density')
    ax_h.set_xlabel('Days Adulthood')

    if not ax_provided:
        return (fig_h, ax_h, metadata)
    else:
        return metadata

def get_featurespan_worm(feature, cutoff, dwell_time=2*24, return_crossings=False):
    '''Get amount of time feature is remains above/below a prescribed cutoff value;
        assumes that all ages are in *hours*.

        Parameters
            feature - str denoting feature of interest in the worms object
            cutoff - float value for cutoff (in units of the feature)
            dwell_time - (optional) minimum time an individual must stay below the
                cutoff to continue in "featurespan" (in hours)
            return_crossings - optional bool denoting whether to return the feature value
                at the time which an individual passed into poor health (good for
                debugging and visualization)

        Returns
            featurespans (in hours)
            (optional, if return_crossings is True) crossing value

    '''

    # TODO: Eliminate origin get_featurespan in favor of this function.
    def feature_func(worm):
        # Ensure that feature decreases with time
        # TODO: Make the features to negate a little less hardcoded?
        if any([label in feature for label in ['af', 'intensity', 'autofluorescence']]):
            negate = True
        else: negate = False

        try:
            ages, feature_data = worm.get_time_range(feature, min_age=0, age_feature='adult_age')
        except AttributeError: # 'adult_age' not in origin Willie-style data
            ages, feature_data = worm.get_time_range(feature, min_age=0, age_feature='egg_age')

        if negate:
            feature_data *= -1

        adjusted_data = feature_data - cutoff
        crossings = np.where((adjusted_data[:-1] > 0) & (adjusted_data[1:] < 0))[0]

        if len(crossings) == 0:
            if adjusted_data[0] > 0:
                try:
                    featurespan = worm.adultspan
                except AttributeError: # TEmporary hack for older data
                    featurespan = worm.lifespan
                crossing_val = feature_data[-1]
            elif adjusted_data[0] <= 0:
                featurespan = 0
                crossing_val = feature_data[0]
        else:
            found_crossing = False
            for crossing_index in crossings:
                if ages[crossing_index] + dwell_time > ages[-1]: # Case where crossing is near death
                    dwell_end = len(ages)
                else:
                    dwell_end = np.where(ages > (ages[crossing_index] + dwell_time))[0][0]

                if (adjusted_data[crossing_index+1:dwell_end] < 0).all():
                    featurespan = ages[crossing_index]
                    crossing_val = feature_data[crossing_index]
                    found_crossing = True
                    break
            if not found_crossing:
                featurespan = ages[crossing_index] # Default to the last crossing in a bad situation
                crossing_val = feature_data[crossing_index]

        if return_crossings:
            return featurespan, crossing_val
        else:
            return featurespan
    return feature_func

def plot_HSvsLS(worms, measure,cutoff_value,
                ax_h=None,return_corrs=False,draw_scatter=False,draw_trend=True,draw_quart=False,
                dwell_time=2*24,**plot_kws):
    '''
        cutoff_value - absolute value for measurement (to pass to get_healthspans_unconfusing)
    '''
    ax_provided = ax_h is not None
    if not ax_provided: fig_h, ax_h = plt.subplots()

    lifespans = worms.get_feature('adultspan')
    sl_cohort = (strain_ls <= np.percentile(strain_ls,25)) #& (strain_ls > np.percentile(strain_ls,1))
    ll_cohort = (strain_ls >= np.percentile(strain_ls,75)) #& (strain_ls < np.percentile(strain_ls,97))

    if measure == 'cumulative_eggs':
        raise Warning('span-comparison not defined for cumulative_eggs')

    healthspans = worms.get_feature(get_featurespan_worm(measure, cutoff_value,dwell_time=dwell_time))/24

    corr_data_wp = scipy.stats.linregress(lifespans, healthspans)
    corr_data_sl = scipy.stats.linregress(lifespans[sl_cohort],healthspans[sl_cohort])
    corr_data_ll = scipy.stats.linregress(lifespans[ll_cohort],healthspans[ll_cohort])

    fraction_LSlimited = ((lifespans-healthspans) <= 1).sum()/len(lifespans)
    fraction_HSnull = (healthspans<=1).sum()/len(healthspans)
    print(("{} - {} (cutoff = {:.3f})\n"
                    "WP: Beta={:.3f},R^2={:.3f} (p={:.3f})\n"
                    "SL: Beta={:.3f},R^2={:.3f} (p={:.3f}) (mean={:.3f})\n"
                    "LL: Beta={:.3f},R^2={:.3f} (p={:.3f}) (mean={:.3f})\n"
                    "LS-limited = {:.3f}, HS-null = {:.3f}\n").format(
        worms[0].name, measure, cutoff_value, #Cutoff in absolute units!
        corr_data_wp[0],corr_data_wp[2]**2, corr_data-wp[3],
        corr_data_sl[0],corr_data_sl[2]**2, corr_data_sl[3], healthspans[sl_cohort].mean(),
        corr_data_ll[0],corr_data_ll[2]**2, corr_data_ll[3], healthspans[ll_cohort].mean(),
        fraction_LSlimited, fraction_HSnull))

    if draw_trend:
        x_vals, y_vals = trend(lifespans,healthspans)
        ax_h.plot(x_vals,y_vals,**plot_kws)

    if draw_quart:
        sl_reg = regress.regress(lifespans[sl_cohort],healthspans[sl_cohort])
        ll_reg = regress.regress(lifespans[ll_cohort],healthspans[ll_cohort])
        ax_h.plot(sl_reg[-1],sl_reg[0],**plot_kws)
        ax_h.plot(ll_reg[-1],ll_reg[0],**plot_kws)

    if draw_scatter and (draw_trend or draw_quart): # Fade the scatter if we're superimposing a trendline
        if 'color' in plot_kws:
            s_pkws = plot_kws.copy()
            my_color = np.array(s_pkws.pop('color'))
            if (my_color==0).all():
                s_pkws['color'] = (1-my_color)*0.6
            else:
                s_pkws['color'] = my_color*0.6
            ax_h.scatter(lifespans, healthspans,**s_pkws)
    elif draw_scatter:
        ax_h.scatter(lifespans, healthspans,**plot_kws)

    ax_h.set_xlabel('Lifespan (Days Adulthood)')
    ax_h.set_ylabel('Healthspan (Days Adulthood)')
    ax_h.set_title(("{} (cutoff = {:.3f})\n").format(
        measure, cutoff_value #Cutoff in absolute units!
        ))
#     ax_h.plot(ax_h.get_xlim(),ax_h.get_xlim(),'--k')

    to_return = []
    if ax_provided:
        to_return.append(ax_h)
    else:
        to_return.extend([fig_h,ax_h])
    if return_corrs: to_return.append([corr_data_wp,corr_data_sl,corr_data_ll])
    return to_return

def peak_size(size_measurement,percentile=99):
    '''
    Calc peak size as the nth percentile (default 99th)
    '''
    def feature(worm):
        size = getattr(worm.td, size_measurement)
        try:
            adult_age = getattr(worm.td, 'adult_age')
        except AttributeError:
            adult_age = getattr(worm.td, 'egg_age') # Hack for older/Willie style data.
        if np.isnan(adult_age).all(): # Not dead yet.
            return np.nan

        nan_mask = np.isnan(size) # Need for measurement pipeline including dead timepoints (come back later to fix)
        return np.percentile(size[(adult_age>0) & (~nan_mask)],percentile)
    return feature
