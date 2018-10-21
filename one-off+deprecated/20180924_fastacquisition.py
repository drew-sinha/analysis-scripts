import multiprocessing
import pathlib

import numpy

from elegant import load_data, segment_images

import elegant_filters


def annotate_poses(experiment_root, to_measure):
    images = load_data.scan_experiment_dir(experiment_root, 
        timepoint_filter=lambda position_n, timepoint_n: position_n in to_measure and timepoint_n in to_measure[position_n][1])
    segment_images.annotate_poses_from_masks(images, pathlib.Path(experiment_root) / 'derived_data' / 'mask', to_measure,) 


def update_poses(experiment_root): 
    annotations = load_data.read_annotations(experiment_root) 
    to_measure = load_data.filter_annotations(annotations, load_data.filter_excluded) 
    to_measure = load_data.filter_annotations(to_measure, elegant_filters.filter_adult_timepoints) 

    n_jobs = 4
    job_position_names = numpy.array_split(list(to_measure.keys()), n_jobs)
    job_positions = [{name: to_measure[name] for name in names} for names in job_position_names]
    job_args = [(experiment_root, to_measure) for job_pos in job_positions]
    with multiprocessing.Pool(processes=n_jobs) as pool:
        pool.starmap(annotate_poses, job_args)
