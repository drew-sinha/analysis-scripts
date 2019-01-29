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
