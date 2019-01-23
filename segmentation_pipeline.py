import pkg_resources
import sys
import pickle
import pathlib
import time
import datetime

from elegant import process_experiment, load_data, segment_images, process_data, worm_widths

import elegant_hacks

def filter_adult_images(experiment_root):
    '''Filter for only adult timepoints from non-excluded animals'''
    experiment_annotations = load_data.read_annotations(experiment_root)
    def scan_filter(position_name, timepoint_name):
        return not experiment_annotations[position_name][0]['exclude'] and experiment_annotations[position_name][1][timepoint_name].get('stage') == 'adult'
    return scan_filter

def process_experiment_with_filter(experiment_root, model, image_filter, mask_root=None, overwrite_existing=False, channels='bf'):
    '''
         image_filter - filter for scan_experiment_dir
    '''

    if mask_root is None:
        mask_root = pathlib.Path(experiment_root) / 'derived_data' / 'mask'

    elegant_hacks.propagate_stages(experiment_root)

    start_t = time.time()
    positions = load_data.scan_experiment_dir(experiment_root,
        timepoint_filter=image_filter, channels=channels)
    scan_t = time.time()
    print(f'scanning done after {(scan_t-start_t)} s') #3 s once, 80s another, taking a while to load up the segmenter....

    process = segment_images.segment_positions(positions, model, mask_root, use_gpu=True,
        overwrite_existing=False)
    if process.stderr:
        print(f'Errors during segmentation: {process.stderr}') #raise Exception)
        #raise Exception()
    segment_t = time.time()
    print(f'segmenting done after {(segment_t-scan_t)} s')
    with (mask_root / 'notes.txt').open('a+') as notes_file:
        notes_file.write(f'{datetime.datetime.today().strftime("%Y-%m-%dt%H%M")} These masks were segmented with model {model}\n')

    annotations = load_data.read_annotations(experiment_root)
    metadata = load_data.read_metadata(experiment_root)
    age_factor = metadata.get('age_factor', 1)
    width_estimator = worm_widths.WidthEstimator.from_experiment_metadata(metadata, age_factor)
    segment_images.annotate_poses_from_masks(positions, mask_root, annotations,
        overwrite_existing, width_estimator)
    load_data.write_annotations(experiment_root, annotations)

    annotation_t = time.time()
    print(f'annotation done after {(annotation_t - segment_t)} s') # ~3.5 hr

if __name__ == "__main__":
    '''Call signature %run segmentation_pipeline.py EXPERIMENT_ROOT MODEL_PATH'''

    experiment_root = sys.argv[1]
    if len(sys.argv) >= 2:
        model = sys.argv[2]
    else:
        model = 'default_CF.mat'

    process_experiment.segment_experiment(experiment_root,model,overwrite_existing=False) # New elegant no longer finds largest components during segmentation; instead this now happens when annotations are updated.
