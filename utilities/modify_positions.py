import pathlib
from elegant import load_data
from zplib import datafile

import json
import time

def reset_positions(scope, experiment_dir, *annotation_filters):
    '''
    Call with annotation filters like so:
    reset_position.reset_positions(scope, experiment_dir, elegant_filters.filter_excluded, elegant_filters.filter_live_animals)
    '''
    
    experiment_dir = pathlib.Path(experiment_dir)
    print(f'Traversing {experiment_dir.name}')
    metadata = load_data.read_metadata(experiment_dir)
    new_positions = {}
    if annotation_filters:
        experiment_annotations = load_data.read_annotations(experiment_dir)
        for filter in annotation_filters:
            experiment_annotations = load_data.filter_annotations(experiment_annotations, filter)
    
    for position_name, position in metadata['positions'].items():
        if annotation_filters and position_name not in experiment_annotations: continue
        
        scope.stage.x = position[0]
        scope.stage.y = position[1]
        try:
            input(f'Press enter for position to replace {position_name} ; ctrl-c to finish')
            new_positions[position_name] = scope.stage.position
        except KeyboardInterrupt:
            break

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
