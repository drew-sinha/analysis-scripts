import pathlib
import collections

import freeimage

from elegant.elegant import load_data, worm_spline

def benchmark_masks_DS(expt_dir):
    expt_dir = pathlib.Path(expt_dir)
    annotations = load_data.read_annotations(expt_dir)
    annotations = load_data.filter_excluded(annotations, load_data.filter_excluded)
    
    ious = collections.OrderedDict()
    for position, position_annotations in annotations.items():
        for timepoint, timepoint_annotations in position_annotations.items()
            center_tck, width_tck = timepoint_annotations.get('pose', (None, None))
            if center_tck is not None and width_tck is not None:
                image_key = position + '_' + timepoint
                mask_image = freeimage.read(str(expt_dir / position / (timepoint + ' bf_mask.png'))) > 0
                manual_mask = worm_spline.lab_frame_mask(center_tck, width_tck, mask_image.shape) > 0
                
                ious[image_key] = (mask_image & manual_mask).sum() / (mask_image | manual_mask).sum()
    return ious
