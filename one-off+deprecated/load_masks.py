import pathlib
import sys

from ris_widget import ris_widget
from elegant import load_data
from elegant.gui import stage_field, experiment_annotator, pose_annotation

import elegant_filters

def load_masks(experiment_root):
    experiment_root = pathlib.Path(experiment_root)
    experiment_annotations = load_data.read_annotations(experiment_root)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_subsample_timepoints(experiment_root))
    experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_adult_timepoints)
    image_filter = elegant_filters.filter_from_elegant_dict(experiment_annotations)

    experiment_images = load_data.scan_experiment_dir(experiment_root, timepoint_filter=image_filter)
    # experiment_images_masks = load_data.scan_experiment_dir(experiment_root / 'derived_data' / 'mask', timepoint_filter=image_filter)

    for position, position_images in experiment_images.items():
        for timepoint, timepoint_images in position_images.items():
            timepoint_images.append(experiment_root / 'derived_data' / 'mask' / position / f'{timepoint} bf.png')

    # experiment_images = experiment_images_bf.copy()
    # for position, position_images in experiment_images.items():
    #     for timepoint, timepoint_images in position_images.items():
    #         timepoint_images.append(experiment_images_masks[position][timepoint])

    return experiment_images

if __name__ == "__main__":
    expt_dir = pathlib.Path(sys.argv[1])

    show_poses = True

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()

    if hasattr(rw, 'annotator'):
        rw.annotator.close()
        del(rw.annotator)

    # measurement_pipeline.propagate_stages(expt_dir)

    experiment_images = load_masks(expt_dir)

    annotation_fields = []
    annotation_fields.append(stage_field.StageField())

    if show_poses:
        metadata = load_data.read_metadata(expt_dir)
        pa = pose_annotation.PoseAnnotation.from_experiment_metadata(metadata, rw)
        annotation_fields.append(pa)

    ea = experiment_annotator.ExperimentAnnotator(rw, expt_dir.parts[-1],
            experiment_images, annotation_fields)
