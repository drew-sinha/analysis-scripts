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

class MultipassMovementMeasurements:
    feature_names = ['summed_multipass_centroid_dist']
    POSE_ANNOTATIONS = ['pose'] + [f'bf_{i+1} pose' for i in range(7)]

    def __init__(self, microns_per_pixel):
        self.microns_per_pixel = microns_per_pixel

    def measure(self, position_root, timepoint, annotations, before, after):
        measures = {}
        centroid_distances = []
        for initial_pose_annotation, next_pose_annotation in zip(self.POSE_ANNOTATIONS[:-1],self.POSE_ANNOTATIONS[1:]):
            initial_center_tck, initial_width_tck = annotations.get(initial_pose_annotation, (None, None))
            next_center_tck, next_width_tck = annotations.get(next_initial_pose_annotation, (None, None))
            if initial_center_tck is None or next_center_tck is None:
                centroid_distances.append(numpy.nan)
                break
            else:
                centroid_distances.append(spline_geometry.centroid_distance(initial_center_tck, next_center_tck, num_points=300))
        measures['summed_multipass_centroid_dist'] = numpy.sum(centroid_distances) * self.microns_per_pixel
        return [measures.get(feature, numpy.nan) for feature in feature_names]

def make_multipass_movement_measurements(experiment_root, update_poses=True, adult_only=True):
    measures = [MultipassMovementMeasurements(microns_per_pixel=1.3)]
    measurement_name = 'multipass_measures'

    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_excluded)

    if adult_only:
        to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_by_stage('adult'))

    if update_poses:
        images = load_data.scan_experiment_dir(experiment_root,
            timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
        segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)


if __name__ == "__main__":
    # Call make_measurements EXPT_DIR
    expt_dir = pathlib.Path(sys.argv[1])
    process_data.update_annotations(expt_dir)
    make_basic_measurements(expt_dir)
    make_movement_measurements(expt_dir)
    make_af_measurements(expt_dir)
    make_multipass_movement_measurements(expt_dir)
    process_data.collate_data(expt_dir)
