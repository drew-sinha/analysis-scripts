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

#===================================
# scan_experiment_dir filters
#==================================

def filter_adult_images(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    def scan_filter(position_name, timepoint_name):
        return experiment_annotations[position_name][1][timepoint_name].get('stage') == 'adult'
    return scan_filter
