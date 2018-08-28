from elegant import load_data

#==============================
# annotation filters
#=============================
def filter_adult_timepoints(position_name, position_annotations, timepoint_annotations):
    if not load_data.filter_excluded(position_name, position_annotations, timepoint_annotations):
        return False
    return [tp.get('stage') == 'adult' for tp in timepoint_annotations.values()]

def filter_by_kw(kw):
    def filter(position_name, position_annotations, timepoint_annotations):
        return kw in position_annotations.get('notes')
    return filter



#===================================
# scan_experiment_dir filters
#==================================
