import pathlib

import freeimage
from zplib.image import mask

from elegant import clean_timepoint_data

def clean_experiment(experiment_root,dry_run):
    clean_timepoint_data.remove_excluded_positions(experiment_root, dry_run=dry_run)
    clean_timepoint_data.remove_dead_timepoints(experiment_root, postmortem_timepoints=16, dry_run=dry_run)

def save_mask_lcc(experiment_root):
    experiment_root = pathlib.Path(experiment_root)
    mask_root = experiment_root / 'derived_data' / 'mask'

    for position_mask_root in sorted(mask_root.iterdir()):
        for mask_file in sorted(positin_mask_root.iterdir()):
            mask_image = freeimage.read(str(mask_file)) > 0
            new_mask = mask.get_largest_object(mask_image).astype(numpy.uint8)
            freeimage.write(new_mask*255, str(mask_file))
