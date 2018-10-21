import pathlib
import sys

from elegant import load_data, process_data,worm_data, segment_images
import elegant_filters

def filter_living_timepoints(position_name, position_annotations, timepoint_annotations):
    """Filter-function for filter_annotations() to return only non-excluded worms
    which have been completely annotated with life stages; of these, all timepoints
    annotated as "egg" or "dead", except the last "egg" and first "dead" will be
    excluded. (The non-excluded "egg" and "dead" allow us to define the hatch and
    death times more carefully.) Any timepoints which have been annotated "exclude"
    will also be excluded."""
    if not load_data.filter_excluded(position_name, position_annotations, timepoint_annotations):
        return False
    stages = [tp.get('stage') for tp in timepoint_annotations.values()]
    # if not all(stages) or stages[-1] != 'dead':
    #     return False
    good_stages = []
    n = len(stages)
    for i, stage in enumerate(stages):
        if stage == 'egg':
            keep = i < n-1 and stages[i+1] != 'egg'
        elif stage == 'dead':
            keep = i > 0 and stages[i-1] != 'dead'
        else:
            keep = True
        good_stages.append(keep)
    excludes = [tp.get('exclude', False) for tp in timepoint_annotations.values()]
    all_good = [good_stage and not exclude for good_stage, exclude in zip(good_stages, excludes)]
    return all_good

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
            elif timepoint_info.get('stage') != running_stage and timepoint_info.get('stage') is not None:
                running_stage = timepoint_info.get('stage')

            if timepoint_info.get('stage') is None and running_stage is not None: # Also handles the case that we are working with an excluded position
                timepoint_info['stage'] = running_stage
                changed.append(timepoint)

        if verbose and changed: print(f'{position_name}: {changed}')
    annotations = load_data.write_annotations(experiment_root, annotations)

def make_basic_measurements(experiment_root):
    measures = [process_data.BasicMeasurements()] #, process_data.PoseMeasurements(microns_per_pixel=5)]
    measurement_name = 'core_measures'

    #process_data.update_annotations(experiment_root)
    propagate_stages(experiment_root,verbose=True)
    positions = load_data.read_annotations(experiment_root)
    # to_measure = load_data.filter_annotations(positions, load_data.filter_living_timepoints)
    to_measure = load_data.filter_annotations(positions, filter_living_timepoints)
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_movement_measurements(experiment_root, update_poses=True, adult_only=True):
    measures = [process_data.PoseMeasurements(microns_per_pixel=1.3)]
    measurement_name = 'pose_measures'

    annotations = load_data.read_annotations(experiment_root)
    if adult_only:
        to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_adult_timepoints)
    else:
        to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_excluded)

    if update_poses:
        images = load_data.scan_experiment_dir(experiment_root, 
            timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
    	segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_af_measurements(experiment_root):
    measures = [process_data.FluorMeasurements('green_yellow_excitation_autofluorescence')]
    measurement_name = 'autofluorescence_measures'

    #process_data.update_annotations(experiment_root)
    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_adult_timepoints)
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def remove_poses(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    for position, position_annotations in experiment_annotations.items():
        timepoint_annotations = position_annotations[1]
        for timepoint, timepoint_annotation in timepoint_annotations.items():
            timepoint_annotation['pose'] = (None, None)
    load_data.write_annotations(experiment_root, experiment_annotations)

if __name__ == "__main__":
    # Call make_measurements EXPT_DIR
    expt_dir = pathlib.Path(sys.argv[1])
    process_data.update_annotations(expt_dir)
    make_basic_measurements(expt_dir)
    make_movement_measurements(expt_dir)
    #make_af_measurements(expt_dir)
    process_data.collate_data(expt_dir)
