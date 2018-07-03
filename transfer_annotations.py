import pathlib
import collections

from elegant import load_data

def read_annotations(experiment_root, annotation_subdir='annotations'):
    """Read annotation data from an experiment directory.

    Parameters:
        experiment_root: the path to an experimental directory.

    Returns: an ordered dictionary mapping position names to annotations,
        where each annotation is a (position_annotations, timepoint_annotations)
        pair. In this, position_annotations is a dict of "global" per-position
        annotation information, while timepoint_annotations is an ordered dict
        mapping timepoint names to annotation dictionaries (which themselves map
        strings to annotation data).

    Example:
        positions = read_annotations('my_experiment')
        position_annotations, timepoint_annotations = positions['009']
        life_stage = timepoint_annotations['2017-04-23t0122']['stage']
    """
    experiment_root = pathlib.Path(experiment_root)
    positions = collections.OrderedDict()
    for annotation_file in sorted(experiment_root.glob(annotation_subdir+'/*.pickle')):
        worm_name = annotation_file.stem
        positions[worm_name] = load_data.read_annotation_file(annotation_file)
    return positions

def transfer_timepoint_annotations(target_annotations, source_annotations, kw):
    for position in target_annotations:
        assert target_annotations[position][1].keys() == source_annotations[position][1].keys()
        for timepoint in target_annotations[position][1]:
            target_annotations[position][1][timepoint][kw] = source_annotations[position][1][timepoint][kw]
    return target_annotations

expt_dir = '/mnt/9karray/Sinha_Drew/20180518_spe-9_Run_3'
DS_annotations = read_annotations(expt_dir)
LT_annotations = read_annotations(expt_dir, 'annotations_larval')
DS_annotations = transfer_timepoint_annotations(DS_annotations, LT_annotations, 'stage')
load_data.write_annotations(expt_dir,DS_annotations)

