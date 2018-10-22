import multiprocessing
import pathlib

import numpy

from elegant import load_data, segment_images

import elegant_filters


def annotate_poses(experiment_root, to_measure):
    images = load_data.scan_experiment_dir(experiment_root,
        timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
    segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure,)
    return to_measure

def update_poses(experiment_root):
    experiment_annotations = load_data.read_annotations(experiment_root)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)
    experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_adult_timepoints)

    n_jobs = 4
    job_position_names = numpy.array_split(list(experiment_annotations.keys()), n_jobs)
    job_positions = [{name: experiment_annotations[name] for name in names} for names in job_position_names]
    job_args = [(experiment_root, job_pos) for job_pos in job_positions]  # Make separate working copy for in-place modification
    with multiprocessing.Pool(processes=n_jobs) as pool:
        annotation_dicts = pool.starmap(annotate_poses, job_args)

    # for annotations in annotation_dicts:  # This was the old. Should work fine either way....
    for annotations in job_positions:
        load_data.merge_annotations(experiment_annotations, annotations)
    load_data.write_annotations(experiment_root, experiment_annotations)
    raise Exception() # Keep this state around in case something doesn't work as advertised
