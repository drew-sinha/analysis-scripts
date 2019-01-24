import json
import pathlib
from collections import OrderedDict

import numpy as np

from elegant.elegant import load_data


def fill_in_stage_annotations(annotations):
    """Fill in annotations where patches of the timepoints have been annotated"""
    
    for position_annotations in annotations:
        stage = None
        for timepoint, timepoint_annotations in position_annotations.items():
            if timepoint_annotations['stage'] == '': 
                timepoint_annotations['stage'] = stage
            if stage is None: stage = timepoint_annotations['stage'] # Handle first timepoint
            if stage != timepoint_annotations['stage']: # Update running stage
                stage = timepoint_annotations['stage']
    return annotations

def prune_empty_entries(od):
    od_copy = od.copy()
    [od_copy.pop(key) for key,val in od.items() if len(val) == 0]
    return od_copy

def filter_fixed_timepoints(experiment_root):
    experiment_root = pathlib.Path(experiment_root)
    
    annotations = load_data.read_annotations(experiment_root)
    good_annotations = load_data.filter_positions(annotations, load_data.filter_excluded)
    def timepoint_filter(position_name, timepoint_name):
        return position_name in good_annotations
    compiled_images = load_data.scan_experiment_dir(experiment_root, 
        timepoint_filter=timepoint_filter)
    compiled_images = prune_empty_entries(compiled_images) # scan_experiment_dir still returns those entries
        
    with (experiment_root / 'experiment_metadata.json').open('r') as em_fp:
        experiment_metadata = json.load(em_fp)
    timepoints = experiment_metadata['timepoints']
    timestamps = np.array(experiment_metadata['timestamps'])/(3600*24) # timestamps in days
    
    selected_images = OrderedDict()
    for (position, position_annotations), position_images in zip(good_annotations.items(),
        compiled_images.values()):
            
            # Grab timepoints flanking adulthood
            start_timepoint = end_timepoint = None
            for timepoint, timepoint_annotations in position_annotations[1].items():
                if timepoint_annotations['stage'] == 'adult' and start_timepoint is None:
                    start_timepoint = timepoint
                elif timepoint_annotations['stage'] == 'dead':
                    end_timepoint = timepoint
                    break
            start_timestamp = timestamps[timepoints.index(start_timepoint)]
            if end_timepoint is not None:
                end_timestamp = timestamps[timepoints.index(end_timepoint)]
            else: # If animal isn't dead yet....
                end_timestamp = timestamps[-1]
            
            days_adulthood = 0
            selected_images[position] = OrderedDict()
            while True:
                timepoint_idx = np.argmin(np.abs(timestamps - (start_timestamp + days_adulthood)))
                if ((timepoint_idx + 1) == len(timepoints)) or timestamps[timepoint_idx] > end_timestamp:
                    break
                selected_images[position][timepoints[timepoint_idx]] = position_images[timepoints[timepoint_idx]]
                days_adulthood += 1
        
    return selected_images
