import pkg_resources
import sys
import pickle
import pathlib
import time

from elegant import process_experiment, load_data, segment_images, process_data

def filter_adult_images(experiment_root):
    '''Filter for only adult timepoints from non-excluded animals'''
    experiment_annotations = load_data.read_annotations(experiment_root)
    def scan_filter(position_name, timepoint_name):
        return not experiment_annotations[position_name][0]['exclude'] and experiment_annotations[position_name][1][timepoint_name].get('stage') == 'adult'
    return scan_filter

def propagate_stages(experiment_root,verbose=False):
    '''
        Modifies experiment annotations by propagating stage information forward
            in time across annotated timepoints.
    '''
    annotations = load_data.read_annotations(experiment_root)
    # annotations = load_data.filter_annotations(annotations,load_data.filter_excluded)
    for position_name, (position_annotations, timepoint_annotations) in annotations.items():
        running_stage = None
        changed = []
        for timepoint,timepoint_info in timepoint_annotations.items():
            if running_stage is None: # Either first timepoint or all the annotations up to now are null
                running_stage = timepoint_info.get('stage')
            elif timepoint_info.get('stage') != running_stage:
                running_stage = timepoint_info.get('stage')

            if timepoint_info.get('stage') is None: # Also handles the case that we are working with an excluded position
                timepoint_info['stage'] = running_stage
                changed.append(timepoint)

        if verbose and changed: print(f'{position_name}: {changed}')
    annotations = load_data.write_annotations(experiment_root, annotations)


def process_experiment_with_filter(experiment_root, model, image_filter):
    propagate_stages(experiment_root)

    start_t = time.time()
    positions = load_data.scan_experiment_dir(experiment_root,
        timepoint_filter=image_filter)
    scan_t = time.time()
    print(f'scanning done after {(scan_t-start_t)} s') #3 s once, 80s another, taking a while to load up the segmenter....
    segment_images.segment_positions(positions, model, use_gpu=True,
        overwrite_existing=False)
    segment_t = time.time()
    print(f'scanning done after {(segment_t-scan_t)} s')


    mask_root = pathlib.Path(experiment_root) / 'derived_data' / 'mask'
    with (mask_root / 'notes.txt').open('w') as notes_file:
        notes_file.write(f'These masks were segmented with model {model}\n')

    process_data.annotate(experiment_root, [process_data.annotate_poses]) # ~5-6 hr for single bf
    annotation_t = time.time()
    print(f'annotation done after {(annotation_t - segment_t)} s') # ~3.5 hr

if __name__ == "__main__":
    experiment_root = sys.argv[1]
    if len(sys.argv) >= 2:
        model = sys.argv[2]
    else:
        model = 'default_CF.mat'

    #process_experiment.segment_images(experiment_root,segmenter_path,overwrite_existing=True) # If you stop the job early, the largest component does get pulled out
    process_experiment.segment_experiment(experiment_root,model,overwrite_existing=False) # New elegant no longer finds largest components during segmentation; instead this now happens when annotations are updated.
