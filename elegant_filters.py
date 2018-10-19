import numpy 
from elegant import load_data

#==============================
# annotation filters
#=============================
def filter_adult_timepoints(position_name, position_annotations, timepoint_annotations):
    # if not load_data.filter_excluded(position_name, position_annotations, timepoint_annotations):
    #     return False
    return [tp.get('stage') == 'adult' for tp in timepoint_annotations.values()]

def filter_by_kw(kw):
    def filter(position_name, position_annotations, timepoint_annotations):
        return kw in position_annotations.get('notes')
    return filter

def filter_by_age(min_age, max_age):
    def filter(position_name, position_annotations, timepoint_annotations):
        return [timepoint_annotation.get('age')/24 >= min_age and timepoint_annotation.get('age')/24 <= max_age
            for timepoint_annotation in timepoint_annotations.values()]
    return filter

def filter_adult_dead_timepoints(position_name, position_annotations, timepoint_annotations):
    # if not load_data.filter_excluded(position_name, position_annotations, timepoint_annotations):
    #     return False
    return [(tp.get('stage') == 'adult') or (tp.get('stage') == 'dead') for tp in timepoint_annotations.values()]

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



'''
Examples for custom misc. filtering
def select_worms(experiment_dir):
    def annotation_filter(position_name, position_annotations, timepoint_annotations):
        worm_selection = {'20180810_age-1_spe-9_Run_3': ['016'],
            '20180816_age-1_spe-9_Run_4': ['017','037','080','009','087', '057']}
        return position_name in worm_selection[pathlib.Path(experiment_dir).name]
    return annotation_filter
'''

#===================================
# scan_experiment_dir filters
#==================================

def filter_adult_images(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    def scan_filter(position_name, timepoint_name):
        return experiment_annotations[position_name][1][timepoint_name].get('stage') == 'adult'
    return scan_filter
