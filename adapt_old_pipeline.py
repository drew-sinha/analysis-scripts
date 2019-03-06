import pathlib
import csv
import json

from elegant import process_data

import elegant_hacks

previous_annotations = ['hatch', 'adult', 'dead']

def compile_annotations_from_tsv(experiment_root):
    experiment_root = pathlib.Path(experiment_root)
    process_data.update_annotations(experiment_root)

    with (experiment_root / 'experiment_metadata.json').open('r') as mdata_file:
        experiment_metadata = json.load(mdata_file)

    annotation_data = {}
    with experiment_root.glob('*Annotations.tsv')[0].open('r') as annotaton_file:
        reader = csv.reader(annotation_file, delimiter='\t')
        _ = reader.__next__() # Header

        for line in reader:
            position = line[0][1:] # Starts with '\'
            for field, frame_num in zip(previous_annotations, line[1:])
                annotation_data[position][field] = experiment_metadata['timepoints'][frame_num]

    annotations = load_data.read_annotations(experiment_root)
    for position, (position_annotations, timepoint_annotations) in annotations.items():
        for field in previous_annotations:
            transition_timepoint = annotation_data[position][field]
            timepoint_annotations[transition_timepoint]['stage'] = field
    load_data.write_annotations(experiment_root, annotations)

    elegant_hacks.propagate_stages(experiment_root) # Need this stage propagation since the stages aren't monotonically sequential
