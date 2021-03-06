import numpy
from elegant import load_data

#==============================
# annotation filters
#=============================

def filter_by_kw(kw):
    def filter(position_name, position_annotations, timepoint_annotations):
        return kw in position_annotations.get('notes')
    return filter

def filter_live_animals(position_name, position_annotations, timepoint_annotations):
    return not any([tp.get('stage') == 'dead' for tp in timepoint_annotations.values()])

def filter_dead_animals(position_name, position_annotations, timepoint_annotations):
    return any([tp.get('stage') == 'dead' for tp in timepoint_annotations.values()])

def filter_by_age(min_age, max_age,adult_age=False):
    '''
        min_age, max_age - age in hours
    '''
    def filter(position_name, position_annotations, timepoint_annotations):
        ages = [info['age'] for info in timepoint_annotations.values()]
        stages = [info['stage'] for info in timepoint_annotations.values()]
        if adult_age:
            starting_timestamp = ages[stages.index('adult')]
            starting_age = starting_timestamp + min_age
            stopping_age = starting_timestamp + max_age
        else:
            starting_age, stopping_age = min_age, max_age

        return [timepoint_annotation.get('age') >= starting_age and timepoint_annotation.get('age') <= stopping_age
            for timepoint_annotation in timepoint_annotations.values()]
    return filter

def filter_by_stage(stages):
    if type(stages) is str:
        stages = [stages]

    def stage_filter(position_name, position_annotations, timepoint_annotations):
        return [tp.get('stage') in stages for tp in timepoint_annotations.values()]
    return stage_filter

def filter_living_timepoints(position_name, position_annotations, timepoint_annotations):
    return (
        not load_data.filter_excluded(position_name, position_annotations, timepoint_annotations) or
        [tp.get('stage') not in ['egg', 'dead'] for tp in timepoint_annotations.values()]
    )

def filter_before_timepoint(timepoint):
    def annotation_filter(position_name, position_annotations, timepoint_annotations):
        return [ann_timepoint <= timepoint for ann_timepoint in timepoint_annotations]
    return annotation_filter

def filter_after_timepoint(timepoint):
    def annotation_filter(position_name, position_annotations, timepoint_annotations):
        return [ann_timepoint >= timepoint for ann_timepoint in timepoint_annotations]
    return annotation_filter

def filter_subsample_timepoints(experiment_dir, interval=3):
    '''
        interval - subsampling interval in hours
    '''
    e_metadata = load_data.read_metadata(experiment_dir)
    timestamps = numpy.array(e_metadata['timestamps'])
    timepoints = e_metadata['timepoints']
    expt_start = timestamps[0]
    expt_end = timestamps[-1]

    step = interval*3600 #sec
    i = 1
    timepoints_to_load = []
    while expt_start + step*i < expt_end:
        timepoints_to_load.append(
            timepoints[abs(timestamps-(expt_start + step*i)).argmin()]
        )
        i += 1

    def filter(position_name, position_annotations, timepoint_annotations):
        if position_annotations['exclude']:
            return False
        return [timepoint in timepoints_to_load for timepoint in timepoint_annotations]
    return filter

def filter_range_before_stage(experiment_dir, time_radius,stage='adult'):
    '''
        time_radius - radius in hours
    '''
    experiment_annotations = load_data.read_annotations(experiment_dir)
    timepoints_to_load = {}
    for position, (position_annotations, timepoint_annotations) in experiment_annotations.items():
        if position_annotations['exclude']:
            continue
        for timepoint, time_annotations in timepoint_annotations.items():
            if 'stage' in time_annotations and time_annotations['stage'] == 'adult':
                first_adult_timestamp = time_annotations['timestamp']
                break

        timepoints_to_load[position] = [timepoint
            for timepoint, time_annotations in timepoint_annotations.items()
            if abs(first_adult_timestamp - time_annotations['timestamp']) < time_radius*3600
                and first_adult_timestamp >= time_annotations['timestamp']]

    def filter(position_name, position_annotations, timepoint_annotations):
        if position_name not in timepoints_to_load:
            return False
        return [timepoint in timepoints_to_load[position_name] for timepoint in timepoint_annotations]
    return filter

def compose_timepoint_filters(*filters):
    """Creates a overall timepoint filter that is the composition of AND'ing multiple individual timepoint filters.

    Parameters:
        filters: Iterable of timepoint filters (with standard load_data.filter_annotations signature)
    """

    def composed_filter(position_name, position_annotations, timepoint_annotations):
        return_val = True
        for filter in filters:
            return_val &= numpy.array(filter(position_name, position_annotations, timepoint_annotations)) # Takes care of single boolean and boolean array return values
        return return_val
    return composed_filter

def select_worms(worm_positions):
    def annotation_filter(position_name, position_annotations, timepoint_annotations):
        return position_name in worm_positions
    return annotation_filter

'''
Examples for custom misc. filtering
def select_worms(experiment_dir):
    def annotation_filter(position_name, position_annotations, timepoint_annotations):
        worm_selection = {'20180810_age-1_spe-9_Run_3': ['016'],
            '20180816_age-1_spe-9_Run_4': ['017','037','080','009','087', '057']}
        return position_name in worm_selection[pathlib.Path(experiment_dir).name]
    return annotation_filter




def filter_from_elegant_worms(worms):
    def annotation_filter(position_name, position_annotations, timepoint_annotations):
        return position_name in [worm.name.split()[1] for worm in worms]
    return annotation_filter
'''

#===================================
# scan_experiment_dir filters
#==================================

def filter_latest_images(experiment_root):
    annotations = load_data.read_annotations(experiment_root)
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    def latelife_filter(position_name, timepoint_name):
        return position_name in good_annotations and timepoint_name > good_annotations[position_name][0]['__last_timepoint_annotated__']
    return load_data.scan_experiment_dir(expt_dir, timepoint_filter=latelife_filter)

def filter_adult_images(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    def scan_filter(position_name, timepoint_name):
        return experiment_annotations[position_name][1][timepoint_name].get('stage') == 'adult'
    return scan_filter

def filter_excluded_images(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    def scan_filter(position_name, timepoint_name):
        return not experiment_annotations[position_name][0]['exclude']
    return scan_filter

def filter_from_elegant_dict(annotation_dict):
    '''Use an annotation dictionary to select specific images for further analysis

        annotation_dict - elegant style (ordered)dict which maps positions to
            tuples of position and timepoint annotations
    '''

    def scan_filter(position_name, timepoint_name):
        return position_name in annotation_dict and timepoint_name in annotation_dict[position_name][1]
    return scan_filter

def compose_scan_filters(*filters):
    def composed_filter(position_name, timepoint_name):
        return_val = True
        for filter in filters:
            return_val &= numpy.array(filter(position_name, timepoint_name)) # Takes care of single boolean and boolean array return values
        return return_val
    return composed_filter
