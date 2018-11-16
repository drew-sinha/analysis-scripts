import multiprocessing
import pathlib
import collections

import numpy

from elegant import load_data, segment_images

import elegant_filters, elegant_hacks


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



#===============================================
# Faster acquisition code (201811 - experiment)

def segment_faster_acquisition(experiment_root):
    elegant_hacks.propagate_stages(experiment_root)

    experiment_annotations = load_data.read_annotations(experiment_root)
    filters = [load_data.filter_excluded, elegant_filters.filter_adult_timepoints]
    for filter in filters:
        experiment_annotations = load_data.filter_annotations(experiment_annotations, filter)
    timepoint_filter = elegant_filters.filter_from_elegant_dict(experiment_annotations)

    mask_root = pathlib.Path(experiment_root) / 'derived_data' / 'mask'

    positions = load_data.scan_all_images(experiment_root)
    filtered_positions = collections.OrderedDict()
    for position_name, timepoints in positions.items():
        filtered_timepoints = collections.OrderedDict()
        for timepoint_name, timepoint_images in timepoints.items():
            if timepoint_filter is None or timepoint_filter(position_name, timepoint_name):
                channel_images = [image_path for channel, image_path in sorted(timepoint_images.items())]
                if len(channel_images) > 0:
                    filtered_timepoints[timepoint_name] = channel_images
        if len(filtered_timepoints) > 0:
            filtered_positions[position_name] = filtered_timepoints
    positions = filtered_positions
    print('done scanning')

    process = segment_images.segment_positions(positions, model, mask_root, use_gpu=True,
        overwrite_existing=False)
    print('done segmenting')
    print(f'segmenting errors: {process.stderr}')


    mask_root = pathlib.Path(experiment_root) / 'derived_data' / 'mask'
    with (mask_root / 'notes.txt').open('w') as notes_file:
        notes_file.write(f'These masks were segmented with model {model}\n')
