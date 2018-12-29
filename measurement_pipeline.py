import pathlib
import sys

from elegant import load_data, process_data,worm_data, segment_images
import elegant_filters

def make_basic_measurements(experiment_root):
    measures = [process_data.BasicMeasurements()] #, process_data.PoseMeasurements(microns_per_pixel=5)]
    measurement_name = 'core_measures'

    elegant_hacks.propagate_stages(experiment_root,verbose=True)
    positions = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(positions, elegant_filters.filter_living_timepoints)
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_movement_measurements(experiment_root, update_poses=True, adult_only=True):
    measures = [process_data.PoseMeasurements(microns_per_pixel=1.3)]
    measurement_name = 'pose_measures'

    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_excluded)

    if adult_only:
        to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_by_stage('adult'))

    if update_poses:
        images = load_data.scan_experiment_dir(experiment_root,
            timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
        segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_af_measurements(experiment_root, fl_measurement_name='green_yellow_excitation_autofluorescence'):
    measures = [process_data.FluorMeasurements(fl_measurement_name)]
    measurement_name = 'autofluorescence_measures'

    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_by_stage('adult'))
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

if __name__ == "__main__":
    # Call make_measurements EXPT_DIR
    expt_dir = pathlib.Path(sys.argv[1])
    process_data.update_annotations(expt_dir)
    make_basic_measurements(expt_dir)
    make_movement_measurements(expt_dir)
    make_af_measurements(expt_dir)
    process_data.collate_data(expt_dir)
