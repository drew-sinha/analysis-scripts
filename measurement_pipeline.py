import pathlib
import sys

import numpy
import skimage.measure as ski_measure

import freeimage
import zplib.image.mask as zpl_mask
from zplib.curve import spline_geometry
from elegant import load_data, process_data, worm_data, segment_images
import elegant_filters, elegant_hacks

def make_basic_measurements(experiment_root):
    measures = [process_data.BasicMeasurements()]
    measurement_name = 'core_measures'

    elegant_hacks.propagate_stages(experiment_root,verbose=True)
    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_excluded)
    #to_measure = load_data.filter_annotations(to_measure, load_data.filter_living_timepoints)
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_pose_measurements(experiment_root, update_poses=False, adult_only=True):
    measures = [process_data.PoseMeasurements(microns_per_pixel=1.3)]
    measurement_name = 'pose_measures'

    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_excluded)

    if adult_only:
        to_measure = load_data.filter_annotations(to_measure, elegant_filters.filter_by_stage('adult'))

    if update_poses:
        images = load_data.scan_experiment_dir(experiment_root,
            timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
        segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_af_measurements(experiment_root, fl_measurement_name='green_yellow_excitation_autofluorescence'):
    measures = [process_data.FluorMeasurements(fl_measurement_name)]
    measurement_name = 'autofluorescence_measures'

    annotations = load_data.read_annotations(experiment_root)
    annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_by_stage('adult'))
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_gfp_measurements(experiment_root, fl_measurement_name='gfp'):
    measures = [process_data.FluorMeasurements(fl_measurement_name)]
    measurement_name = 'fluorescence_measures'

    annotations = load_data.read_annotations(experiment_root)
    annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    to_measure = load_data.filter_annotations(annotations, elegant_filters.filter_by_stage('adult'))
    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

class MultipassPoseMeasurements:
    feature_names = ['summed_multipass_centroid_dist', 'multipass_length', 'multipass_max_width', 'multipass_area', 'multipass_volume']
    POSE_ANNOTATIONS = ['pose'] + [f'bf_{i+1} pose' for i in range(7)]

    def __init__(self, microns_per_pixel):
        self.microns_per_pixel = microns_per_pixel

    def measure(self, position_root, timepoint, annotations, before, after):
        measures = {}
        centroid_distances = []
        lengths, widths, areas, volumes = [], [], [], []

        first_center_tck, first_width_tck = annotations.get(self.POSE_ANNOTATIONS[0], (None, None))
        if first_center_tck is None or first_width_tck is None:
            lengths.append(numpy.nan)
            widths.append(numpy.nan)
            areas.append(numpy.nan)
            volumes.append(numpy.nan)
            centroid_distances.append(numpy.nan)
        else:
            volume, surface_area = spline_geometry.volume_and_surface_area(first_center_tck, first_width_tck)
            area = spline_geometry.area(first_center_tck, first_width_tck)
            length, max_width = spline_geometry.length_and_max_width(first_center_tck, first_width_tck)

            lengths.append(length * self.microns_per_pixel)
            widths.append(max_width * 2 * self.microns_per_pixel)
            areas.append(area * self.microns_per_pixel**2)
            volumes.append(volume * self.microns_per_pixel**3)

            for initial_pose_annotation, next_pose_annotation in zip(self.POSE_ANNOTATIONS[:-1],self.POSE_ANNOTATIONS[1:]):
                initial_center_tck, initial_width_tck = annotations.get(initial_pose_annotation, (None, None))
                next_center_tck, next_width_tck = annotations.get(next_pose_annotation, (None, None))

                if initial_center_tck is None or next_center_tck is None:
                    lengths.append(numpy.nan)
                    widths.append(numpy.nan)
                    areas.append(numpy.nan)
                    volumes.append(numpy.nan)
                    centroid_distances.append(numpy.nan)
                    break
                else:
                    volume, surface_area = spline_geometry.volume_and_surface_area(next_center_tck, next_width_tck)
                    length, max_width = spline_geometry.length_and_max_width(next_center_tck, next_width_tck)
                    area = spline_geometry.area(next_center_tck, next_width_tck)

                    lengths.append(length * self.microns_per_pixel)
                    widths.append(max_width * 2 * self.microns_per_pixel)
                    areas.append(area * self.microns_per_pixel**2)
                    volumes.append(volume * self.microns_per_pixel**3)
                    centroid_distances.append(spline_geometry.centroid_distance(initial_center_tck, next_center_tck, num_points=300))
        measures['summed_multipass_centroid_dist'] = numpy.sum(centroid_distances) * self.microns_per_pixel
        measures['multipass_length'] = numpy.median(lengths)
        measures['multipass_max_width'] = numpy.median(widths)
        measures['multipass_area'] = numpy.median(areas)
        measures['multipass_volume'] = numpy.median(volumes)

        return [measures.get(feature, numpy.nan) for feature in self.feature_names]

def make_multipass_measurements(experiment_root, update_poses=False, adult_only=True):
    measures = [MultipassPoseMeasurements(microns_per_pixel=1.3)]
    measurement_name = 'multipass_measures'

    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_excluded)

    if update_poses:
        images = load_data.scan_experiment_dir(experiment_root,
            timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
        segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

class MaskPoseMeasurements:
    """Provide data columns based on annotated worm pose information.

    Given the pose data, each worm's length, volume, surface_area, and maximum
    width, plus the area of the 2D projection of the worm into the image plane
    (i.e. the area of the worm region in the image).

    If no pose annotation is present, Nones are returned.

    Note: the correct microns_per_pixel conversion factor MUST passed to the
    constructor of this class.
    """
    feature_names = ['mask_area', 'mask_centroid_dist']

    def __init__(self, microns_per_pixel, pose_annotation='pose'):
        self.microns_per_pixel = microns_per_pixel
        self.pose_annotation = pose_annotation
        self.mask_name = 'bf'

    def get_mask(self, position_root, derived_root, timepoint, annotations):
        mask_file = derived_root / 'mask' / position_root.name / f'{timepoint} {self.mask_name}.png'
        if not mask_file.exists():
            #raise Exception()
            print(f'No mask file found for {position_root.name} at {timepoint}.')
            return None
        else:
            mask = freeimage.read(mask_file)
            if mask.sum() == 0:
                print(f'No worm region defined for {position_root.name} at {timepoint}')
                return None
            else:
                mask = zpl_mask.get_largest_object(mask, structure=numpy.ones((3,3)))
                return mask

    def measure(self, position_root, timepoint, annotations, before, after):
        print(f'Measuring position {position_root.name} - {timepoint}')
        measures = {}
        derived_root = position_root.parent / 'derived_data'

        mask = self.get_mask(position_root, derived_root, timepoint, annotations)
        if mask is None:
            return [numpy.nan] * len(self.feature_names)

        measures['mask_area'] = mask.sum() * self.microns_per_pixel**2

        moments = ski_measure.moments(mask, order=1)
        centroid = numpy.array([moments[1,0] / moments[0,0], moments[0,1] / moments [0,0]])

        centroid_distances = []
        for adjacent in (before, after):
            if adjacent is not None:
                adj_mask = self.get_mask(position_root, derived_root, adjacent['timepoint'], annotations)
                if adj_mask is None:
                    centroid_distances.append(numpy.nan)
                    break
                adj_moments = ski_measure.moments(adj_mask, order=1)
                adj_centroid = numpy.array([adj_moments[1,0] / adj_moments[0,0], adj_moments[0,1] / adj_moments [0,0]])
                adj_dist = ((centroid - adj_centroid)**2).sum()**0.5
                centroid_distances.append(adj_dist)
        measures['mask_centroid_dist'] = numpy.sum(centroid_distances) * self.microns_per_pixel
        return [measures[feature] for feature in self.feature_names]

def annotate_timepoints(experiment_root, position, timepoint, metadata, annotations):
    annotations['timepoint'] = metadata['timepoint']

def make_mask_measurements(experiment_root, update_poses=False):
    process_data.annotate(experiment_root, annotators=[annotate_timepoints])

    measures = [MaskPoseMeasurements(microns_per_pixel=1.3)]
    measurement_name = 'mask_measures'

    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_excluded)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_living_timepoints)
    to_measure = load_data.filter_annotations(to_measure, elegant_filters.filter_by_stage('adult'))

    if update_poses:
        images = load_data.scan_experiment_dir(experiment_root,
            timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
        segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def make_lawn_measurements(experiment_root, remake_lawns=False):
    if not (pathlib.Path(experiment_root) / 'derived_data' / 'lawn_masks').exists() or remake_lawns:
        process_data.annotate(experiment_root, position_annotators=[process_data.annotate_lawn])

    measures = [process_data.LawnMeasurements()]
    measurement_name = 'lawn_measures'
    
    annotations = load_data.read_annotations(experiment_root)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_excluded)
    to_measure = load_data.filter_annotations(annotations, load_data.filter_living_timepoints)

    process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

def run_canonical_measurements(experiment_dir):
    '''Run standard measurements on the specified experiment directory'''
    experiment_dir = pathlib.Path(experiment_dir)

    process_data.update_annotations(experiment_dir)

    position_features = ['stage_x','stage_y','starting_stage_z','notes']
    annotations = load_data.read_annotations(experiment_dir)
    annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    if any(['lawn_area' in position_annotations for (position_annotations, timepoint_annotations) in annotations.items()]):
        position_features.append('lawn_area')

    make_basic_measurements(experiment_dir)
    make_pose_measurements(experiment_dir)

    process_data.collate_data(experiment_dir, position_features=position_features)

    make_mask_measurements(experiment_dir)

    image_channels = elegant_hacks.get_image_channels(experiment_dir)
    print(f'Image channels: {image_channels}')

    if 'bf_1' in image_channels:
        print('Found multipass movement channel bf_1; making measurements')
        make_multipass_measurements(experiment_dir, update_poses=False)

    process_data.collate_data(experiment_dir,position_features=position_features) # For convenience since autofluorescence can take a little while....

    if 'green_yellow_excitation_autofluorescence' in image_channels or 'autofluorescence' in image_channels:
        fl_measurement_name = 'autofluorescence' if 'autofluorescence' in image_channels else 'green_yellow_excitation_autofluorescence'
        print(f'Found autofluorescence channel {fl_measurement_name}; making measurements')

        make_af_measurements(experiment_dir, fl_measurement_name=fl_measurement_name)

    process_data.collate_data(experiment_dir, position_features=position_features)


if __name__ == "__main__":
    # Call make_measurements EXPT_DIR
    expt_dir = pathlib.Path(sys.argv[1])
    run_canonical_measurements(expt_dir)
