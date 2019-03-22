import pathlib
import csv
import json
import shutil

from elegant import process_data, load_data

import elegant_hacks

previous_annotations = ['larva', 'adult', 'dead']

def compile_annotations_from_tsv(experiment_root):
    experiment_root = pathlib.Path(experiment_root)
    #raise Exception()
    process_data.update_annotations(experiment_root)

    with (experiment_root / 'experiment_metadata.json').open('r') as mdata_file:
        experiment_metadata = json.load(mdata_file)

    annotation_data = {}
    with list(experiment_root.glob('*Annotations.tsv'))[0].open('r') as annotation_file:
        reader = csv.reader(annotation_file, delimiter='\t')
        _ = reader.__next__() # Header

        for line in reader:
            position = line[0][1:] # Starts with '\'
            notes = line[-1]

            annotation_data[position] = {}
            if 'DEAD' in notes:
                for field, frame_num in zip(previous_annotations, line[1:]):
                    annotation_data[position][field] = experiment_metadata['timepoints'][int(frame_num)]
            annotation_data[position]['notes'] = line[-1]


    annotations = load_data.read_annotations(experiment_root)
    for position, (position_annotations, timepoint_annotations) in annotations.items():
        if 'DEAD' in annotation_data[position]['notes']:
            first_timepoint = list(timepoint_annotations.keys())[0]
            timepoint_annotations[first_timepoint]['stage'] = 'egg'

            for field in previous_annotations:
                transition_timepoint = annotation_data[position][field]
                timepoint_annotations[transition_timepoint]['stage'] = field
            position_annotations['exclude'] = False
        else:
            position_annotations['exclude'] = True
        position_annotations['notes'] = annotation_data[position]['notes']
    load_data.write_annotations(experiment_root, annotations)

    elegant_hacks.propagate_stages(experiment_root) # Need this stage propagation since the stages aren't monotonically sequential

def fill_empty_stages(experiment_root):
    '''
       Need this since the first pass of compiling annotations didn't do its job correctly.
    '''
    annotations = load_data.read_annotations(experiment_root)
    for position, (position_annotations, timepoint_annotations) in annotations.items():
        if not position_annotations['exclude']:
            for timepoint, timepoint_annotation in timepoint_annotations.items():
                if 'stage' not in timepoint_annotation:
                    timepoint_annotation['stage'] = 'egg'
    load_data.write_annotations(experiment_root, annotations)

def move_great_lawn(experiment_root, remove_lawn=False):
    experiment_root = pathlib.Path(experiment_root)
    (experiment_root / 'derived_data' / 'great_lawns').mkdir(parents=True,exist_ok=True)
    for position_root in [p.parent for p in experiment_root.glob('*/position_metadata.json')]:
        #print(position_root)
        if (position_root / 'great_lawn.png').exists():
            shutil.copyfile(str(position_root / 'great_lawn.png'), str(experiment_root / 'derived_data' / 'great_lawns' / f'{position_root.name}.png'))
        if remove_lawn:
            (position_root / 'great_lawn.png').unlink()

def prep_experiment_for_processing(experiment_root):
    compile_annotations_from_tsv(experiment_root)
    move_great_lawn(experiment_root, remove_lawn=False)
