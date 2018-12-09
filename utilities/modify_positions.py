import json
import time

import pathlib
import numpy

from elegant import load_data
from zplib import datafile

import elegant_filters

'''
TODO:
1. Better/figure out how to include z-updating
2. Do fine adjustment by doing imaging + lawn_finding to get robust landmarks

'''

def poll_positions(scope, experiment_metadata, positions):
    new_positions = {}
    for position_name in positions:
        scope.stage.x = experiment_metadata[position_name][0]
        scope.stage.y = experiment_metadata[position_name][1]
        
        try:
            input(f'Press enter for position to replace {position_name} ; ctrl-c to abort')
            new_positions[position_name] = scope.stage.position
        except KeyboardInterrupt:
            raise
    return new_positions

def reset_positions(scope, experiment_dir, *annotation_filters, positions=None):
    '''
    Call with annotation filters like so:
    reset_position.reset_positions(scope, experiment_dir, elegant_filters.filter_excluded, elegant_filters.filter_live_animals)
    '''
    
    experiment_dir = pathlib.Path(experiment_dir)
    print(f'Traversing {experiment_dir.name}')
    metadata = load_data.read_metadata(experiment_dir)
    if annotation_filters:
        experiment_annotations = load_data.read_annotations(experiment_dir)
        for filter in annotation_filters:
            experiment_annotations = load_data.filter_annotations(experiment_annotations, filter)
        positions = experiment_annotations.keys()

    if new_positions:
        try:
            new_positions = poll_positions(scope, metadata, positions)
            input(f'\nPress any key to save positions; ctrl-c to abort')
            time_label = time.strftime('%Y%m%d-%H%M-%S')

            with (experiment_dir/f'experiment_metadata_beforechangingpositions_{time_label}.json').open('w') as mdata_file:
                datafile.json_encode_legible_to_file(metadata,mdata_file)

            metadata['positions'].update(new_positions)
            load_data.write_metadata(metadata, experiment_dir)
        except KeyboardInterrupt:
            pass
    else:
        print('No positions found to reset')

def automatic_reset_positions(scope, experiment_dir):
    experiment_dir = pathlib.Path(experiment_dir)
    print(f'Automatically resetting positions for {experiment_dir.name}')
    
    metadata = load_data.read_metadata(experiment_dir)
    old_positions = metadata['positions'].copy()
    old_positions = {position: numpy.array(position_coords) for position, position_coords in old_positions}

    experiment_annotations = load_data.read_annotations(experiment_dir)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_live_animals)
    base_position = list(experiment_annotation.keys())[0]

    farthest_position, max_distance = base_position, 0
    for position in experiment_annotations:
        distance = ((old_positions[base_position][:-1] - old_positions[base_position][:-1])**2).sum()**0.5
        if distance > max_distance:
            max_distance = distance
            farthest_position = position
    positions = [base_position, farthest_position]
    
    new_positions = {position: numpy.array(position_coords) 
        for position, position_coords in poll_positions(scope, metadata, positions).items()}
    
    old_orientation = old_positions[farthest_position] - old_positions[base_position]
    new_orientation = new_positions[farthest_position] - new_positions[base_position]
    old_orientation, new_orientation = old_orientation[:-1], new_orientation[:-1]
    angle_offset = numpy.arctan(*new_orientation) - numpy.arctan(*old_orientation)
    rotation_matrix = numpy.array([[numpy.cos(angle_offset), -numpy.sin(angle_offset)],[numpy.sin(angle_offset), numpy.cos(angle_offset)]])

    adjusted_positions = {}
    for position, position_coords in old_positions:
        adjusted_positions[position] = list(
            rotation_matrix @ (position_coords[:-1] - old_positions[base_position][:-1]) + new_positions[base_position][:-1])
        adjusted_positions[position].append(old_positions[position][-1])
    
    input(f'\nPress any key to save positions; ctrl-c to abort')
    time_label = time.strftime('%Y%m%d-%H%M-%S')

    with (experiment_dir/f'experiment_metadata_beforechangingpositions_{time_label}.json').open('w') as mdata_file:
        datafile.json_encode_legible_to_file(metadata,mdata_file)

    metadata['positions'].update(adjusted_positions)
    load_data.write_metadata(metadata, experiment_dir)
