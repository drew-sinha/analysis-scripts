import pathlib
import sys

from collections import OrderedDict

from ris_widget import ris_widget
from elegant import load_data
from elegant.gui import experiment_annotator, stage_field

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
    print(len(annotations))
    good_annotations = load_data.filter_positions(annotations, filter_good)
    print(len(good_annotations))
    def latelife_filter(position_name, timepoint_name):
        return timepoint_name > datestr and position_name in good_annotations
    return load_data.scan_experiment_dir(expt_dir, timepoint_filter=latelife_filter)

def scan_latest(expt_dir):
    annotations = load_data.read_annotations(expt_dir)
    print(len(annotations))
    good_annotations = load_data.filter_positions(annotations, filter_good)
    print(len(good_annotations))
    def latelife_filter(position_name, timepoint_name):
        return position_name in good_annotations and timepoint_name > good_annotations[position_name][0]['__last_timepoint_annotated__']
    return load_data.scan_experiment_dir(expt_dir, timepoint_filter=latelife_filter)

def check_for_alive(expt_dir):
    annotations = load_data.read_annotations(expt_dir)
    print(len(annotations))
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
    print(len(good_annotations))
    dead_annotations = check_stage_annotations(good_annotations, ['dead'])

    print(f'{len(good_annotations)-len(dead_annotations)}/{len(good_annotations)} still alive')
    
    return set(good_annotations.keys()).difference(set(dead_annotations.keys()))

if __name__ == "__main__":
    '''Call as python run_annotator.py EXPT_DIR'''
    
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
    
    ea = experiment_annotator.ExperimentAnnotator(rw, expt_dir.parts[-1],
        expt_pos, [sf])
