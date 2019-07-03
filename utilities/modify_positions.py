import pathlib
from elegant import load_data, segment_images
from zplib import datafile

import json
import time

def poll_positions(scope, experiment_metadata, positions):
    new_positions = {}
    for position_name in positions:
        scope.stage.x = experiment_metadata['positions'][position_name][0]
        scope.stage.y = experiment_metadata['positions'][position_name][1]
        
        try:
            input(f'Press enter for position to replace {position_name} ; ctrl-c to finish')
            new_positions[position_name] = scope.stage.position
        except KeyboardInterrupt:
            break
    return new_positions

def reset_positions_manual(scope, experiment_dir, *annotation_filters, positions=None):
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
    else:
        positions = metadata['positions'].keys()

    new_positions = poll_positions(scope, metadata, positions)

    if new_positions:
        try:
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
        
def reset_positions_with_offset(experiment_dir, offset):
    ''' Modify position coordinates based on a fixed x-y offset
    
    Parameters:
        experiment_dir - str/pathlib.Path to experiment root
        offset - list giving the x-y offset
    '''
    
    experiment_dir = pathlib.Path(experiment_dir)
    print(f'Modifying positions for {experiment_dir.name}')
    metadata = load_data.read_metadata(experiment_dir)
    new_metadata = metadata.copy()
    
    if len(offset) == 2 :
        offset.extend([0])

    try:
        time_label = time.strftime('%Y%m%d-%H%M-%S')
        
        for position in metadata['positions']:
            position_coords = metadata['positions'][position]
            new_metadata['positions'][position] = [position_coords[0] + offset[0], position_coords[1] + offset[1], position_coords[2] + offset[2]]

        with (experiment_dir/f'experiment_metadata_beforechangingpositions_{time_label}.json').open('w') as mdata_file:
            datafile.json_encode_legible_to_file(metadata, mdata_file)

        load_data.write_metadata(new_metadata, experiment_dir)
    except KeyboardInterrupt:
        pass

#def automatic_reset_positions(scope, experiment_dir):
    #experiment_dir = pathlib.Path(experiment_dir)
    #print(f'Automatically reseting positions for {experiment_dir.name}')
    
    #metadata = load_data.read_metadata(experiment_dir)
    #experiment_annotations = load_data.read_annotations(experiment_dir)
    #experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)
    #experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_live_animals)
    #base_position = list(experiment_annotation.keys())[0]

    #farthest_position, max_distance = base_position, 0
    #for position in experiment_annotations:
        #distance = ((numpy.array(metadata['positions'][base_position][:-1]) - numpy.array(metadata['positions'][position][:-1]))**2).sum()**0.5
        #if distance > max_distance:
            #max_distance = distance
            #farthest_position = position
    #positions = [base_position, farthest_position]
    
    #old_coordinates = {position: metadata['positions'][position] for position in positions}
    #new_coordinates = poll_positions(scope, metadata, positions)
    
    #offset = numpy.array(new_coordinates[base_position][:-1]) - numpy.array(old_coordinates[base_position][:-1])
    
    #adjusted_coordinates = {position: list(numpy.array(old_coordinates[position]) + offset) for position in old_coordinates}

    #adjusted_coordinates = experiment_metadata['positions'].copy()
    #adjusted_coordinates = {position_coordinates[0] + offset[0]

