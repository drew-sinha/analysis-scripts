import pathlib
import csv
import json
import shutil
import sys

from elegant import process_data, load_data, process_experiment

previous_annotations = ['larva', 'adult', 'dead']

def compile_annotations_from_tsv(experiment_root):
    if type(experiment_root) is str:
        experiment_root = pathlib.Path(experiment_root.replace('\\ ', ' '))
    _check_metadata_for_timepoints(experiment_root)
    process_data.update_annotations(experiment_root)

    with (experiment_root / 'experiment_metadata.json').open('r') as mdata_file:
        experiment_metadata = json.load(mdata_file)

    annotation_data = {}
    with list(experiment_root.glob('*.tsv'))[0].open('r') as annotation_file:
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

    # Note: process_data.propagate_worm_stages assumes that stages are monotonically increasing.
    # To make use of propagate_worm_stages, prepopulate only the timepoints where there's a transition.
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
    process_data.annotate(experiment_root,position_annotators=[process_data.propagate_worm_stage])

def _check_metadata_for_timepoints(experiment_root):
    pm_files = list(experiment_root.glob('*/position_metadata.json'))
    with pm_files[0].open('r') as pm_fp:
        position_metadata = json.load(pm_fp)
    if 'timepoint' not in position_metadata[0]:
        experiment_metadata = load_data.read_metadata(experiment_root)
        for pm_file in pm_files:
            position_root = pm_file.parent
            if not (position_root / 'position_metadata_original.json').exists():
                shutil.copyfile(
                    str(pm_file), 
                    str(position_root / 'position_metadata_original.json'))
            with pm_file.open('r') as pm_fp:
                position_metadata = json.load(pm_fp)
            for metadata_entry, timepoint in zip(position_metadata, experiment_metadata['timepoint']):
                # Was there a bug with the purging code that breaks this?
                metadata_entry['timepoint'] = timepoint
            with pm_file.open('w') as pm_fp:
                json.dump(position_metadata, pm_fp)

def move_great_lawn(experiment_root, remove_lawn=False):
    if type(experiment_root) is str:
        experiment_root = pathlib.Path(experiment_root.replace('\\ ', ' '))
    (experiment_root / 'derived_data' / 'great_lawns').mkdir(parents=True,exist_ok=True)
    for position_root in [p.parent for p in experiment_root.glob('*/position_metadata.json')]:
        #print(position_root)
        if (position_root / 'great_lawn.png').exists():
            shutil.copyfile(str(position_root / 'great_lawn.png'), str(experiment_root / 'derived_data' / 'great_lawns' / f'{position_root.name}.png'))
            if remove_lawn:
                (position_root / 'great_lawn.png').unlink()

def prep_experiment_for_processing(experiment_root, nominal_temperature=None):
    move_great_lawn(experiment_root, remove_lawn=True)
    compile_annotations_from_tsv(experiment_root)
    if nominal_temperature:
        process_experiment.auto_update_metadata_file(experiment_root, nominal_temperature)

if __name__ == "__main__":
    # signature: python adapt_old_pipeline EXPERIMENT_ROOT [NOMINAL TEMPERATURE]
    experiment_root = pathlib.Path(sys.argv[1].replace('\\ ', ' ')) # Handle spaces in older experiment names
    print(experiment_root)
    if len(sys.argv) > 2:
        assert sys.argv[2] in ['20','25']
        nominal_temperature = sys.argv[2]
    else:
        nominal_temperature = None
    prep_experiment_for_processing(experiment_root, nominal_temperature=nominal_temperature)
