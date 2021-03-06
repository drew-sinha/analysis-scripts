import pathlib
import collections

import numpy
import matplotlib.cm

from elegant import worm_data, load_data
from zplib.image import colorize

import plotting_tools, elegant_filters

#=================================
# Preprocessing and reading data
#==================================

def check_stage_annotations(annotations, stages):
    """Check for incomplete annotations in an experiment

        Parameters
            annotations - an OrderedDict of experiment_annotations as returned
                by load_data.read_annotations
            stages - A iterable containing the stages that should be annotated
                (e.g. could be ('larva','adult','dead') for a complete experiment, 
                but only ('larva', 'adult') for an ongoing experiment)
        Returns
            annotations OrderedDict for animals with incomplete annotations
    """

    # Create a suitable function to use with filter_positions using a closure
    def select_by_stage_annotation(position_name,position_annotations, timepoint_annotations):
        stage_annotations = [timepoint_annotation.get('stage','')
            for timepoint_annotation in timepoint_annotations.values()]
        return all([stage in stage_annotations for stage in stages])

    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    complete_annotations = load_data.filter_annotations(
        good_annotations,
        select_by_stage_annotation) # Get positions whose stages are not all annotated

    return complete_annotations

def check_experiment_stages(experiment_dir, stages, verbose=False):
    """Check an experiment for incomplete annotations

        Parameters
            experiment_dir - str/pathlib.Path to an experiment directory
            stages - an iterable containing stages that should be annotated

        Returns
            if verbose is False, return a list of positions with incomplete annotations,
                otherwise, return the positions and their annotations
    """

    annotations = load_data.read_annotations(experiment_dir)
    complete_annotations = check_stage_annotations(annotations, stages)
    incomplete_annotations = {position: annotations[position] for position in set(annotations) - set(complete_annotations)}

    if verbose:
        return incomplete_annotations
    else:
        return incomplete_annotations.keys()

def check_for_alive(expt_dir):
    annotations = load_data.read_annotations(expt_dir)
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    dead_annotations = check_stage_annotations(good_annotations, ['dead'])

    print(f'{len(good_annotations)-len(dead_annotations)}/{len(good_annotations)} still alive')

    return set(good_annotations.keys()).difference(set(dead_annotations.keys()))

def check_for_kw(expt_dir, kw, filter_good=True,verbose=True):
    annotations = load_data.read_annotations(expt_dir)
    if filter_good:
        annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    kw_annotations = load_data.filter_annotations(annotations, elegant_filters.filter_by_kw(kw))
    if verbose:
        print(f'{len(kw_annotations)}/{len(annotations)} of animals in experiment has kw {kw} {"(minus excluded)" if filter_good else ""}')
    return set(kw_annotations.keys())

def show_position_notes(experiment_dir):
    assert pathlib.Path(experiment_dir).exists()
    experiment_annotations = load_data.read_annotations(experiment_dir)
    for position, (position_annotations, timepoint_annotations) in experiment_annotations.items():
        print(f'{position} (excluded = {position_annotations["exclude"]}): {position_annotations["notes"]}')

def check_for_null_poses(experiment_root, annotation_dir='annotations'):
    assert pathlib.Path(experiment_root).exists()
    experiment_annotations = load_data.read_annotations(experiment_root, annotation_dir=annotation_dir)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)

    poses = ['pose']
    for position, position_annotations in experiment_annotations.items():
        for timepoint, timepoint_annotations in position_annotations[1].items():
            if 'bf_1 pose' in timepoint_annotations and 'bf_1 pose' not in poses:
                for i in range(7):
                    poses.append(f'bf_{i+1} pose')

            for pose_tag in poses:
                if timepoint_annotations.get(pose_tag, (None, None))[0] is None and timepoint_annotations['stage'] == 'adult':
                    print(f"Position {position}, timepoint {timepoint} doesn't have a vaild {pose_tag} pose")
    print(f'Checked for poses {poses}')

def replace_annotation(experiment_root, annotation_type, old_annotation_values, new_annotation_value, annotation_dir='annotations'):
    if not isinstance(old_annotation_values, collections.Iterable):
        old_annotation_values = list(old_annotation_values)
    if isinstance(old_annotation_values, str):
        old_annotation_values = [old_annotation_values]

    experiment_annotations = load_data.read_annotations(experiment_root, annotation_dir=annotation_dir)
    for position, position_annotations in experiment_annotations.items():
        for timepoint, timepoint_annotations in position_annotations[1].items():
            if annotation_type in timepoint_annotations and timepoint_annotations[annotation_type] in old_annotation_values:
                timepoint_annotations[annotation_type] = new_annotation_value
    load_data.write_annotations(experiment_root, experiment_annotations, annotation_dir=annotation_dir)

def remove_poses(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    for position, position_annotations in experiment_annotations.items():
        timepoint_annotations = position_annotations[1]
        for timepoint, timepoint_annotation in timepoint_annotations.items():
            timepoint_annotation['pose'] = (None, None)
    load_data.write_annotations(experiment_root, experiment_annotations)

def propagate_stages(experiment_root,verbose=False):
    '''
        Modifies experiment annotations by propagating stage information forward
            in time across annotated timepoints.
        Somewhat deprecated by process_data.propagate_worm_stages/update_annotations;
            however, useful for the case in which one wants to propagate stages
            and can't assume that stages are monotonically increasing with time.

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
# Image loading
#==================================
def load_derived_images(experiment_root, derived_dir, *additional_filters):
    experiment_root = pathlib.Path(experiment_root)
    experiment_annotations = load_data.read_annotations(experiment_root)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)

    for filter in additional_filters:
        experiment_annotations = load_data.filter_annotations(experiment_annotations, filter)
    image_filter = elegant_filters.filter_from_elegant_dict(experiment_annotations)

    experiment_images = load_data.scan_experiment_dir(experiment_root, timepoint_filter=image_filter)
    for position, position_images in experiment_images.items():
        for timepoint, timepoint_images in position_images.items():
            timepoint_images.append(experiment_root / 'derived_data' / derived_dir / position / f'{timepoint} bf.png')

    return experiment_images

def get_image_channels(experiment_root):
    experiment_root = pathlib.Path(experiment_root)
    annotations = load_data.read_annotations(experiment_root)
    annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    positions = list(annotations.keys())

    image_channels = {image_file.stem.split()[1]
        for image_file in (experiment_root / positions[0]).iterdir()
        if image_file.suffix[1:] in ['png', 'tif']}
    return image_channels

#==================================
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

    import matplotlib.pyplot as plt; plt.ion(); plt.show()

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



def scatter_features(worms, x_feature, y_feature, color_by='lifespan'):
    """Plot values of a given feature for each worm, colored by a given
    worm feature (defaults to lifespan).

    Parameters:
        x_feature, y_feature: name/callables for two features to compare against each other
        color_by: worm feature to use for color scale of each timecourse.
    """
    def _feature_plot_data(worms, x_feature, y_feature, color_by='lifespan'):
        x_feature_vals = worms.get_feature(x_feature)
        y_feature_vals = worms.get_feature(y_feature)
        color_vals = colorize.scale(worms.get_feature(color_by), output_max=1)
        colors = colorize.color_map(color_vals, uint8=False)
        out = []
        for x, y, color in zip(x_feature_vals, y_feature_vals, colors):
            out.append((x, y, color))
        return out

    import matplotlib.pyplot as plt; plt.ion(); plt.show()
    for x, y, c in _feature_plot_data(worms, x_feature, y_feature, color_by=color_by):
        plt.scatter(x, y, color=c)
