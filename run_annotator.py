import pathlib
import sys

from collections import OrderedDict

from ris_widget import ris_widget
from elegant import load_data
from elegant.gui import experiment_annotator, stage_field, pose_annotation

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

def scan_late_life(expt_dir, datestr):
    annotations = load_data.read_annotations(expt_dir)
    good_annotations = load_data.filter_positions(annotations, load_data.filter_excluded)
    def latelife_filter(position_name, timepoint_name):
        return timepoint_name > datestr and position_name in good_annotations
    return load_data.scan_experiment_dir(expt_dir, timepoint_filter=latelife_filter)

def scan_latest(expt_dir):
    annotations = load_data.read_annotations(expt_dir)
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    def latelife_filter(position_name, timepoint_name):
        return position_name in good_annotations and timepoint_name > good_annotations[position_name][0]['__last_timepoint_annotated__']
    return load_data.scan_experiment_dir(expt_dir, timepoint_filter=latelife_filter)

def check_for_alive(expt_dir):
    annotations = load_data.read_annotations(expt_dir)
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    dead_annotations = check_stage_annotations(good_annotations, ['dead'])

    print(f'{len(good_annotations)-len(dead_annotations)}/{len(good_annotations)} still alive')

    return set(good_annotations.keys()).difference(set(dead_annotations.keys()))

if __name__ == "__main__":
    '''Call as python run_annotator.py EXPT_DIR [MODE]'''

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()
        rw.show()

    if hasattr(rw, 'annotator'):
        rw.annotator.close()
        del(rw.annotator)

    sf = stage_field.StageField()

    expt_dir = pathlib.Path(sys.argv[1])
    expt_pos = load_data.scan_experiment_dir(expt_dir)

    # if (len(sys.argv) == 1) or (sys.argv[2] == 'adult'):
    #     sf = stage_field.StageField()
    # elif sys.argv[2] == 'larval':
    #     stages = ['egg','L1','L2','L3','L4','young_adult','adult','dead']
    #     transitions = ['hatch','L1m','L2m','L3m','L4m','egg-laying','death']
    #     shortcuts = ['h','1','2','3','4','v','d']
    #     sf = stage_field.StageField(stages=stages,transitions=transitions,shortcuts=shortcuts)

    # width_estimator, width_pca_basis = pose_annotation.default_width_data(pixels_per_micron=1/1.3, experiment_temperature=20)
    # pa = pose_annotation.PoseAnnotation(rw, mean_widths=width_estimator, width_pca_basis=width_pca_basis)

    ea = experiment_annotator.ExperimentAnnotator(rw, expt_dir.parts[-1],
        expt_pos, [sf])
