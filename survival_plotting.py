import pandas as pd
import plotting_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import pathlib

from elegant import worm_data, load_data

from utilities import utilities

def plot_spanseries(spans,ax_h = None,**plot_kws):
    """Helper function that plots a survival curve specified by a list
            of times to an event (e.g. lifespans)

        Parameters
            spans - list/numpy array of time to endpoints
            ax_h - optional matplotlib axis to plot survival curve on
            kws - optional kws to pass to matplotlib.pyplot.plt
    """
    if ax_h is None:
        fig_h, ax_h = plt.subplots(1,1)
        ax_provided = False
    else:
        ax_provided = True
    sorted_s = np.sort(spans[~np.isnan(spans)]) #handling the nan eliminates potential leak-through of alive animals when used with annotation code
    prop_alive = (1 - (np.arange(start=0,stop=np.size(sorted_s)))/(np.size(spans)-1))

    # Generate all the points needed for plotting
    # Plot the first point (0,1), then a straight line across and down for each event/death
    x_vals = np.sort(np.concatenate(([0],sorted_s,sorted_s)))
    y_vals = np.sort(np.concatenate(([1,1],prop_alive[:-1],prop_alive[:-1],[prop_alive[-1]])))[::-1] # Reverse this.

    data = ax_h.plot(x_vals,y_vals,**plot_kws)
    ax_h.set_ylim([0,1])

    if not ax_provided:
        return [fig_h,ax_h,data]
    else:
        return data

def plot_manual_individual_ls(*annotation_fns,ax_h=None):
    """Plotting for plotting survival curves from manual lifespan
        experiments where single individuals were tracked
        (i.e. not automated scope experiments);

        Parameters:
            annotation_fns - one or more annotation tsv files containing
                lifespans for a given experiment; files should contain
                at least the following columns: (Worm, Death, Notes)
                where "Worm" is a unique identifier, "Death" is the timepoint
                (int) when an individual died, and "Notes" contains
                additional information on each individual tracked; if no
                lifespan is provided for a given animal, it is assumed to
                have been censored
            ax_h - optional matplotlib axis to plot lifespans on

        Returns:
            data_series - list of matplotlib.plotseries objects plotted on ax_h
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
    compiled_ann_data = [pd.read_csv(ann_file,sep='\t') for ann_file in annotation_fns]
    data_series = []
    metadata = []

    if ax_h is None:
        fig_h, ax_h = plt.subplots(1,1)
        ax_provided = False
    else: ax_provided = True

    for ann_data, ann_color in zip(compiled_ann_data,plotting_tools.qual_colors):
        good_worms = ~np.isnan(ann_data['Death']) # Excludes alive/non-annotated/bad worms!
        good_ls = ann_data['Death'][good_worms]
        data_series.append(plot_spanseries(good_ls, ax_h=ax_h,color=ann_color))
        metadata.append({'n':len(good_ls),
            'mean': np.nanmean(good_ls),
            'median': np.nanmedian(good_ls),
            'ls': good_ls,})
    ax_h.set_ylim([0,1])
    ax_h.set_xlabel('Days Post-Transfer')
    ax_h.set_ylabel('Proportion Surviving')

    if not ax_provided:
        return (fig_h, ax_h, data_series, metadata)
    else:
        return (data_series, metadata)

def plot_expt_ls(*expt_dirs, ax_h=None,import_from_timecourse=True,calc_adultspan=False,**plot_kws):
    """Plot survival curves for one or more separate experiments

        Parameters
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

    if type(expt_dirs[0]) is str:
        expt_dirs = [pathlib.Path(expt_dir) for expt_dir in expt_dirs]

    legend_entries = []
    metadata = []
    for expt_dir in expt_dirs:
        expt_name = expt_dir.name

        if import_from_timecourse:
            timecourse_file = expt_dir / 'derived_data' / f'{expt_name} timecourse.tsv'
            expt_worms = worm_data.read_worms(timecourse_file)
            if calc_adultspan:
                lifespans = expt_worms.get_feature('adultspan')/24
            else:
                lifespans = expt_worms.get_feature('lifespan')/24 #Convert to days
        else:
            experiment_annotations = load_data.read_annotations(expt_dir)
            experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)

            lifespans = np.array([])
            for position, position_annotations in experiment_annotations.items():
                general_annotations, timepoint_annotations = position_annotations
                timepoints = list(timepoint_annotations.keys())
                life_stages = [timepoint_info.get('stage') for timepoint_info in timepoint_annotations.values()]

                if 'dead' in life_stages:
                    if calc_adultspan:
                        birth_timepoint = timepoints[life_stages.index('adult')]
                    else:
                        birth_timepoint = timepoints[life_stages.index('larva')]

                    death_timepoint = timepoints[life_stages.index('dead')]
                    lifespans = np.append(lifespans, (utilities.extract_datetime_fromstr(death_timepoint) - utilities.extract_datetime_fromstr(birth_timepoint)).total_seconds()/(3600*24))
                else: # Catch animals that are still alive....
                    lifespans = np.append(lifespans, np.nan)
            if any(np.isnan(lifespans)):
                print(f'Warning! Some worms in {expt_name} still alive.')
        data_series = plot_spanseries(lifespans,ax_h=ax_h,**plot_kws)
        expt_metadata = {'n':len(lifespans),
            'mean': np.nanmean(lifespans),
            'median': np.nanmedian(lifespans),
            'ls': lifespans,}
        metadata.append(expt_metadata)
        legend_entries.append(f'{expt_name} (n={len(lifespans)})')
    ax_h.legend(legend_entries)
    ax_h.set_ylabel('Proportion Surviving')
    if calc_adultspan:
        ax_h.set_xlabel('Days Adulthood')
    else:
        ax_h.set_xlabel('Days Post-Hatch')

    if not ax_provided:
        return (fig_h, ax_h, metadata)
    else:
        return metadata


def plot_annotation_ls(*annotation_dirs, ax_h=None,calc_adultspan=False,**plot_kws):
    """Plot survival curves for one or more separate experiments

        Parameters
            annotation_dirs - one or more annotation root directories
            ax_h - optional matplotlib axis objects to plot curves on
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

    if type(annotation_dirs[0]) is str:
        annotation_dirs = [pathlib.Path(anno_dir) for anno_dir in annotation_dirs]

    legend_entries = []
    metadata = []
    for anno_dir in annotation_dirs:
        expt_name = anno_dir.name

        experiment_annotations = load_data.read_annotations(anno_dir.parent,annotation_dir=anno_dir.name)
        experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)

        lifespans = np.array([])
        for position, position_annotations in experiment_annotations.items():
            general_annotations, timepoint_annotations = position_annotations
            timepoints = list(timepoint_annotations.keys())
            life_stages = [timepoint_info.get('stage') for timepoint_info in timepoint_annotations.values()]

            if 'dead' in life_stages:
                if calc_adultspan:
                    birth_timepoint = timepoints[life_stages.index('adult')]
                else:
                    birth_timepoint = timepoints[life_stages.index('larva')]

                death_timepoint = timepoints[life_stages.index('dead')]
                lifespans = np.append(lifespans, (utilities.extract_datetime_fromstr(death_timepoint) - utilities.extract_datetime_fromstr(birth_timepoint)).total_seconds()/(3600*24))
            else: # Catch animals that are still alive....
                lifespans = np.append(lifespans, np.nan)
        if any(np.isnan(lifespans)):
            print(f'Warning! Some worms in {expt_name} still alive.')
        data_series = plot_spanseries(lifespans,ax_h=ax_h,**plot_kws)
        expt_metadata = {'n':len(lifespans),
            'mean': np.nanmean(lifespans),
            'median': np.nanmedian(lifespans),
            'ls': lifespans,}
        metadata.append(expt_metadata)
        legend_entries.append(f'{expt_name} (n={len(lifespans)})')
    ax_h.legend(legend_entries)
    ax_h.set_ylabel('Proportion Surviving')
    if calc_adultspan:
        ax_h.set_xlabel('Days Adulthood')
    else:
        ax_h.set_xlabel('Days Post-Hatch')

    if not ax_provided:
        return (fig_h, ax_h, metadata)
    else:
        return metadata


def plot_expt_ls_old(*expt_dirs, ax_h=None, calc_adultspan=False,bad_worm_kws=[], **plot_kws):
    '''Plot lifespans corresponding to one or more single longitudinal scope
        experiments using the old setup; i.e. a .tsv file of annotations
        in an experiment root containing "Hatch", "First Egg Laid", and "Death"
    '''
    import corral_annotations.annotation_file as annotation_file

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
