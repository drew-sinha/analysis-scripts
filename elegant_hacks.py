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


#==================================
# Plotting timecourses
#==================================

def _timecourse_plot_data(worms, feature, min_age=-numpy.inf, max_age=numpy.inf, age_feature='age', color_by='lifespan',
    color_map=None):

    if color_map is None:
        color_map = lambda array: colorize.color_map(array, uint8=False)

    time_ranges = worms.get_time_range(feature, min_age, max_age, age_feature)
    if color_by is None:
        color_vals = [0]*len(times_range)
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
