import pathlib
import sys

from collections import OrderedDict

import numpy

from ris_widget import ris_widget
from elegant import load_data, worm_widths
from elegant.gui import experiment_annotator, stage_field, pose_annotation

import elegant_hacks

def check_stage_annotations(annotations, stages):
    """Check that a set of annotations are complete

        Parameters
            annotations - An OrderedDict mapping position names to corresponding
                annotations (returned by load_data.read_annotations)
            stages - A iterable containing the stages that should be annotated
                for this experiment (e.g. could be ('larva','adult','dead')
                for a complete experiment, but only ('larva', 'adult') for
                an ongoing experiment)
        Returns
            bad_positions - a list of positions with incomplete annotations
    """

    # Create a suitable function to use with filter_positions using a closure
    def select_by_stage_annotation(position_name,position_annotations, timepoint_annotations):
        stage_annotations = [timepoint_annotation.get('stage','')
            for timepoint_annotation in timepoint_annotations.values()]
        return all([stage in stage_annotations for stage in stages])

    return load_data.filter_annotations(
        annotations,
        select_by_stage_annotation) # Get positions whose stages are not all annotated

def check_for_alive(expt_dir):
    annotations = load_data.read_annotations(expt_dir)
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    dead_annotations = check_stage_annotations(good_annotations, ['dead'])

    print(f'{len(good_annotations)-len(dead_annotations)}/{len(good_annotations)} still alive')

    return set(good_annotations.keys()).difference(set(dead_annotations.keys()))

if __name__ == "__main__":
    '''Call as python run_annotator.py EXPT_DIR [MODE]'''
    import elegant_filters
    expt_dir = pathlib.Path(sys.argv[1])

    show_poses = False
    adult_only = False
    annotation_dir = 'annotations'


    # additional_filters = [elegant_filters.filter_by_age(9,10)]
    additional_filters = [load_data.filter_excluded,elegant_filters.filter_live_animals, elegant_filters.filter_adult_timepoints] #, , elegant_filters.filter_after_timepoint('2018-10-10t1529')] #, , elegant_filters.filter_after_timepoint('2018-11-12t1200')  #[elegant_filters.filter_subsample_timepoints(expt_dir)]#elegant_filters.filter_range_before_stage(expt_dir, 3)] #load_data.filter_excluded] #[select_worms(expt_dir)] # [elegant_filters.filter_adult_dead_timepoints]#load_data.filter_excluded]
    channels = ['bf'] #, 'gfp', 'autofluorescence'] #, 'green_yellow_excitation_autofluorescence']

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()

    if hasattr(rw, 'annotator'):
        rw.annotator.close()
        del(rw.annotator)

    elegant_hacks.propagate_stages(expt_dir)

    experiment_annotations = load_data.read_annotations(expt_dir, annotation_dir=annotation_dir)
    # experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)

    if adult_only:
        experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_adult_timepoints) #time_a.get('stage') == 'adult')

    if additional_filters:
        for position_filter in additional_filters:
            experiment_annotations = load_data.filter_annotations(experiment_annotations, position_filter)

    expt_pos = load_data.scan_experiment_dir(expt_dir,channels=channels,timepoint_filter = lambda position_name, timepoint_name: position_name in experiment_annotations and timepoint_name in experiment_annotations[position_name][1])

    # if (len(sys.argv) == 1) or (sys.argv[2] == 'adult'):
    #     sf = stage_field.StageField()
    # elif sys.argv[2] == 'larval':
    #     stages = ['egg','L1','L2','L3','L4','young_adult','adult','dead']
    #     transitions = ['hatch','L1m','L2m','L3m','L4m','egg-laying','death']
    #     shortcuts = ['h','1','2','3','4','v','d']
    #     sf = stage_field.StageField(stages=stages,transitions=transitions,shortcuts=shortcuts)

    annotation_fields = []
    annotation_fields.append(stage_field.StageField())
    if show_poses:
        metadata = load_data.read_metadata(expt_dir)
        pa = pose_annotation.PoseAnnotation.from_experiment_metadata(metadata, rw)
        annotation_fields.append(pa)

    ea = experiment_annotator.ExperimentAnnotator(rw, expt_dir.parts[-1],
            expt_pos, annotation_fields,annotation_dir=annotation_dir)
