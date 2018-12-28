import pathlib
import sys

from ris_widget import ris_widget
from elegant import load_data
from elegant.gui import experiment_annotator, stage_field, pose_annotation

import elegant_hacks, elegant_filters

if __name__ == "__main__":
    '''Call as python run_annotator.py EXPT_DIR'''
    expt_dir = pathlib.Path(sys.argv[1])

    show_poses = False
    annotation_dir = 'annotations'

    timepoint_filters = [load_data.filter_excluded,elegant_filters.filter_live_animals, elegant_filters.filter_adult_timepoints] #, , elegant_filters.filter_after_timepoint('2018-10-10t1529')] #, , elegant_filters.filter_after_timepoint('2018-11-12t1200')  #[elegant_filters.filter_subsample_timepoints(expt_dir)]#elegant_filters.filter_range_before_stage(expt_dir, 3)] #load_data.filter_excluded] #[select_worms(expt_dir)] # [elegant_filters.filter_adult_dead_timepoints]#load_data.filter_excluded]
    channels = ['bf'] #, 'gfp', 'autofluorescence'] #, 'green_yellow_excitation_autofluorescence'] # First one is the one used to load poses when specified.

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()

    if hasattr(rw, 'annotator'):
        rw.annotator.close()
        del(rw.annotator)

    elegant_hacks.propagate_stages(expt_dir)

    experiment_annotations = load_data.read_annotations(expt_dir, annotation_dir=annotation_dir)

    if timepoint_filters:
        experiment_annotations = load_data.filter_annotations(
            experiment_annotations,
            elegant_filters.compose_timepoint_filters(*timepoint_filters))

    expt_pos = load_data.scan_experiment_dir(expt_dir,channels=channels,timepoint_filter = lambda position_name, timepoint_name: position_name in experiment_annotations and timepoint_name in experiment_annotations[position_name][1])

    annotation_fields = []
    annotation_fields.append(stage_field.StageField())
    if show_poses:
        metadata = load_data.read_metadata(expt_dir)
        if channels[0] == 'bf':
            pose_name = 'pose'
        else:
            pose_name = f'{channels[0]} pose'
        pa = pose_annotation.PoseAnnotation.from_experiment_metadata(metadata, rw, name=pose_name)
        annotation_fields.append(pa)

    ea = experiment_annotator.ExperimentAnnotator(rw, expt_dir.parts[-1],
            expt_pos, annotation_fields,annotation_dir=annotation_dir)
