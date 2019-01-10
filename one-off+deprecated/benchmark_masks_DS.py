import pathlib
import collections

import freeimage

from elegant import load_data, worm_spline

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
                mask_image = freeimage.read(str(expt_dir / position / (timepoint + ' bf_mask.png'))) > 0 # Saved mask made by segmenter
                manual_mask = worm_spline.lab_frame_mask(center_tck, width_tck, mask_image.shape) > 0 # Mask regenerated from manual annotations
                
                ious[image_key] = (mask_image & manual_mask).sum() / (mask_image | manual_mask).sum()
    return ious

def overlay_masks(rw, position_directory):
    position_directory = pathlib.Path(position_directory)
    expt_dir = position_directory.parent
    position_annotations = load_data.read_annotations(expt_dir)[position_directory.name]
    
    files_to_load = []
    page_names = []
    global_positions, timepoint_annotations = position_annotations
    for timepoint, timepoint_data in timepoint_annotations.items():
        image_key = position_directory.name + '_' + timepoint
        image = freeimage.read(str(position_directory / (timepoint + ' bf.png')))
        mask_file = expt_dir / 'derived_data' / 'mask' / position_directory.name / (timepoint + ' bf.png')
        if mask_file.exists():
            mask_image = freeimage.read(str(expt_dir / 'derived_data' / 'mask' / position_directory.name / (timepoint + ' bf.png'))) > 0
            files_to_load.append([image, mask_image])
        else:
            files_to_load.append([image])
    rw.flipbook_pages = files_to_load


'''
# Example code for plotting a distribution of ious
from zplib.scalar_stats import kde
import matplotlib.pyplot as plt

fig_h, ax_h = plt.subplots()
ax_h.set_title('Accuracy of LRR pipeline')
iou_vals = np.array([v for k,v in ious.items()])
ax_h.plot(*kde.kd_distribution(iou_vals)[:-1], color = 'k')
ax_h.set_xlim([0,1])
ax_h.set_xlabel('IoU')
ax_h.set_ylabel('Density')
'''