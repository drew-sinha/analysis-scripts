import pathlib
import collections

import numpy
import matplotlib.pyplot as plt; plt.ion(); plt.show()
import matplotlib.cm

from elegant import worm_data, load_data
from zplib.image import colorize

import plotting_tools

#=================================
# Preprocessing and reading data
#==================================

def read_annotations(experiment_root, annotation_subdir='annotations'):
    """Read annotation data from an experiment directory.
    Parameters:
        experiment_root: the path to an experimental directory.
        annotation_subdir: subdirectory containing annotations of interest
    Returns: an ordered dictionary mapping position names to annotations,
        where each annotation is a (position_annotations, timepoint_annotations)
        pair. In this, position_annotations is a dict of "global" per-position
        annotation information, while timepoint_annotations is an ordered dict
        mapping timepoint names to annotation dictionaries (which themselves map
        strings to annotation data).
    Example:
        positions = read_annotations('my_experiment')
        position_annotations, timepoint_annotations = positions['009']
        life_stage = timepoint_annotations['2017-04-23t0122']['stage']
    """
    experiment_root = pathlib.Path(experiment_root)
    annotation_root = experiment_root / annotation_subdir
    positions = collections.OrderedDict()
    for annotation_file in sorted(annotation_root.glob('*.pickle')):
        worm_name = annotation_file.stem
        positions[worm_name] = load_data.read_annotation_file(annotation_file)
    return positions

def propagate_stages(experiment_root,verbose=False):
    '''
        Modifies experiment annotations by propagating stage information forward
            in time across annotated timepoints.
    '''
    annotations = load_data.read_annotations(experiment_root)
    for position_name, (position_annotations, timepoint_annotations) in annotations.items():
        running_stage = None
        changed = []
        encountered_stages = []

        for timepoint,timepoint_info in timepoint_annotations.items():
            already_encountered = timepoint_info.get('stage') in encountered_stages
            stage_set = timepoint_info.get('stage') is not None

            if running_stage is None: # Either first timepoint or all the annotations up to now are null
                running_stage = timepoint_info.get('stage')
            elif timepoint_info.get('stage') != running_stage and stage_set and not already_encountered:
                running_stage = timepoint_info['stage']

            if stage_set and not already_encountered:
                encountered_stages.append(timepoint_info['stage'])

            if not stage_set and running_stage is not None: # Also handles the case that we are working with an excluded position
                timepoint_info['stage'] = running_stage
                changed.append(timepoint)
            elif stage_set and timepoint_info['stage'] != running_stage and already_encountered:
                timepoint_info['stage'] = running_stage
                changed.append(timepoint)

        if verbose and changed: print(f'{position_name}: {changed}')
    annotations = load_data.write_annotations(experiment_root, annotations)

#====================
# Worm stuff
#==================================

def calculate_ages_and_spans(w):
        """Calculate ages and spans (in hours) from annotated stage data.

        Requires 'stage' and 'timestamp' timecourse data fields. Calculates the
        following timecourse data:
            age
            adult_age
            ghost_age

        Calculates the following summary data:
            lifespan
            [xxx]_span for each non-egg, non-dead stage
        """
        hours = (w.td.timestamp - w.td.timestamp[0]) / 3600
        stages, indices = numpy.unique(w.td.stage, return_index=True)
        order = indices.argsort()
        stages = stages[order]
        indices = indices[order]
        if stages[0] == 'egg':
            # had an egg-visible timepoint: assume hatch was halfway between the
            # last egg-seen time and first larva-seen time.
            hatched_i = indices[1]
            hatch_time = hours[hatched_i-1:hatched_i+1].mean()
            stages = stages[1:]
            indices = indices[1:]
        else:
            # no egg timepoint. Best guess for hatch time is the first larva-seen time.
            hatch_time = hours[0]
        w.td.age = hours - hatch_time
        transition_times = [hatch_time] + [hours[i-1:i+1].mean() for i in indices[1:]]
        transition_times = numpy.array(transition_times)
        spans = transition_times[1:] - transition_times[:-1]
        for stage, span in zip(stages[:-1], spans):
            setattr(w, f'{stage}span', span)
        try:
            adult_i = list(stages).index('adult')
            adult_time = transition_times[adult_i]
            w.td.adult_age = hours - adult_time
        except ValueError:
            # raise ValueError('No timepoint with "adult" label is present; cannot calculate adult_age.')
            pass

        if stages[-1] == 'dead':
            # raise ValueError('No timepoint with "dead" label is present; cannot calculate lifespan.')
            death_time = transition_times[-1]
            w.td.ghost_age = hours - death_time
            w.lifespan = death_time - hatch_time

def collate_data(experiment_root, position_features=('stage_x', 'stage_y', 'starting_stage_z')):
    """Gather all .tsv files produced by measurement runs into a single file.

    This function will concatenate all individual-worm .tsv files for all of the
    different measure_worms runs (which each output their .tsv files into a
    different subdirectory of '{experiment_root}/derived_data/measurements')
    into a single master-file of timecourse data:
        {experiment_root}/derived_data/measurements/{experiment_root.name} timecourse.tsv
    If possible, lifespans and other spans will be calculated for the worms,
    with the results stored in a master-file of summary data:
        {experiment_root}/derived_data/measurements/{experiment_root.name} summary.tsv
    Any features named in the position_features parameter will be transfered
    from the annotations for that position to the worm summary data as well.

    The worms in these files will be renamed as:
        '{experiment_root.name} {position_name}'
    """
    experiment_root = pathlib.Path(experiment_root)
    positions = load_data.read_annotations(experiment_root)
    experiment_name = experiment_root.name
    derived_root = experiment_root / DERIVED_ROOT
    measurement_root = derived_root / 'measurements'
    measurements = []
    name_prefix = experiment_name + ' '
    for measurement_dir in measurement_root.iterdir():
        files = list(measurement_dir.glob('*.tsv'))
        if len(files) > 0:
            measurements.append(worm_data.read_worms(*files, name_prefix=name_prefix, calculate_lifespan=False))
    worms = measurements[0]
    for other_measurement in measurements[1:]:
        worms.merge_in(other_measurement)
    for w in worms:
        try:
            calculate_ages_and_spans(w)
        except (NameError, ValueError):
            print(f'could not calculate lifespan for worm {w.name}')
        position_annotations, timepoint_annotations = positions.get(w.name[len(name_prefix):], ({}, {}))
        for feature in position_features:
            if feature in position_annotations:
                setattr(w, feature, position_annotations[feature])

    worms.write_timecourse_data(derived_root / f'{experiment_name} timecourse.tsv', multi_worm_file=True, error_on_missing=False)
    worms.write_summary_data(derived_root / f'{experiment_name} summary.tsv', error_on_missing=False)


#==================================
# Plotting timecourses
#==================================

def _timecourse_plot_data(worms, feature, min_age=-numpy.inf, max_age=numpy.inf, age_feature='age', color_by='lifespan',
    color_map=None):

    if color_map is None:
        color_map = lambda array: colorize.color_map(array, uint8=False)

    time_ranges = worms.get_time_range(feature, min_age, max_age, age_feature)
    if color_by is None:
        color_vals = [0]*len(time_ranges)
    else:
        color_vals = colorize.scale(worms.get_feature(color_by), output_max=1)
    colors = color_map(color_vals)
    out = []
    for time_range, color in zip(time_ranges, colors):
        x, y = time_range.T
        out.append((x, y, color))
    return out

def plot_timecourse(worms, feature, min_age=-numpy.inf, max_age=numpy.inf,
    age_feature='age', time_units='hours', color_by='lifespan',
    color_map=None, **plot_kws):
    """Plot values of a given feature for each worm, colored by a given
    worm feature (defaults to lifespan).

    Parameters:
        feature: the name of a timecourse feature available in worm.td to
            retrieve (such as "gfp_95th", say), or a function that will be
            called as feature(worm) that will calculate a timecourse feature
            (see examples in Worm.get_time_range).
        min_age, max_age: the beginning and end of the window in which to
            retrieve features. If not specified, features from the very
            beginning and/or to the very end of the timecourse will be
            retrieved.
        age_feature: the name of the feature to use to determine the age
            window and the plot axes. Defaults to 'age', but other
            monotonic features, such as a time-left-alive 'ghost_age' could
            also be used. If this is a function, it will be called as
            age_feature(worm) to generate the ages to examine (see
            examples in Worm.get_time_range).
        time_units: one of "days", "hours", "minutes", or "seconds",
            representing the units in which the time scale should be plotted.
        color_by: worm feature to use for color scale of each timecourse.
    """
    if time_units not in worm_data.TIME_UNITS:
        raise ValueError(f"'time_units' must be one of: {list(TIME_UNITS)}")

    time_scale = worm_data.TIME_UNITS[time_units]
    fig_h, ax_h = plt.subplots()
    ax_h.set_xlabel(f'{age_feature} ({time_units})')
    ax_h.set_ylabel(f'{feature}')
    for x, y, c in _timecourse_plot_data(worms, feature, min_age, max_age, age_feature, color_by,color_map):
        ax_h.plot(x/time_scale, y, color=c, **plot_kws)
    return fig_h, ax_h

'''
    # Code for monochromatic colormap....

    base_color = [] # FILL in.
    cmap = cm.LinearlySegmentedColormap('Monochrome_lab',
        plotting_tools.build_gradient_palette(base_color,256))  # Assume palette is fine enough to linearly segment
    plot_timecourse(worms, feature, min_age=0, age_feature='adult_age', time_units='days', color_map=cmap)
'''
