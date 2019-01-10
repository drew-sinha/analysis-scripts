import pathlib
import numpy as np

from zplib.image import write_movie
import freeimage
from elegant import load_data, process_images

import elegant_filters

def yield_rgb(image_generator):
    for image in image_generator:
        yield np.repeat(image, 3, axis=2)

def shrink(image_generator, factor=2, fast=False):
    """Shrink images produced by a generator by the specified factor.

    Parameters:
        factor: amount to shrink the image by (fold-change)
        fast: if True and if factor is an integer, perform no smoothing.
            If False, smooth the image before downsampling to avoid aliasing.
    """
    if fast:
        fast = factor = int(factor) # Need this....?
    go_fast = fast and int_factor == factor
    for image in image_generator:
        if go_fast:
            yield image[::int(factor), ::int(factor)]
        else:
            yield pyramid.pyr_down(image, factor).astype(numpy.uint8)

# image_dir = pathlib.Path('/mnt/9karray/Sinha_Drew/20180810_spe-9_Control/012/')
# image_filepaths = sorted(image_dir.glob('*bf.png'))
# image_generator = write_movie.generate_images_from_files(image_filepaths)
# image_generator = write_movie.shrink(yield_rgb(image_generator), factor=4) #pyramid.pyr_down returning [x,y,1] shaped image leads to an assertion error.... Should the assertion be around?
# write_movie.write_movie(image_generator, 'test.mp4',framerate=5)


# Helper to combine images and add both to each; also add timepoint+relevant age measurement
def double_image_layout(image_generator1, image_generator2, *info):
    for image1, image2 in zip(image_generator1, image_generator2):
        assert (image1.shape == image2.shape).all()
        combined_image = numpy.zeros((image1.shape[0], 2*image1.shape[1]))
        combined_image[:image1.shape[0]] = image1
        combined_image[image1.shape[0]:] = image2
        yield combined_image

# Actual code to pull relevant images
experiment_root = pathlib.Path('/mnt/9karray/Sinha_Drew/20180924_spe-9_fastmovement')
out_dir = pathlib.Path('/home/drew')
worm = '00'
min_age = 12    # hr
max_age = 24

experiment_annotations = load_data.read_annotations(experiment_root)
experiment_annotations = load_data.filter_annotations(experiment_annotations, load_data.filter_excluded)
experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_adult_timepoints)
experiment_annotations = load_data.filter_annotations(experiment_annotations, elegant_filters.filter_by_age(min_age,max_age,adult_age=True))
image_filter = elegant_filters.filter_from_elegant_dict(experiment_annotations)

# scan_experiment_dir for bf images
print('Scanning...')
bf_images = load_data.scan_experiment_dir(experiment_root, timepoint_filter=image_filter)
bf_generator = write_movie.generate_images_from_files(bf_images['00'])
rep_image = freeimage.read(bf_images['00'][list(bf_images['00'].keys())[0]])

print('Making lab frame mask...')
# scan_experiment_dir for mask images or helper to grab pose masks
position_annotations, timepoint_annotations = experiment_root[worm]
mask_generator = (worm_spline.lab_frame_mask(timepoint_info['bf pose'][0], timepoint_info['bf pose'][1], rep_image.shape) for timepoint_info in timepoint_annotations)

image_generator = write_movie.shrink(yield_rgb(double_image_layout(bf_generator, mask_generator)), factor=4)
write_movie.write_movie(image_generator, out_dir / (f'{worm} {min_age/24:.1f}-{max_age/24:.1f}.replace('.','_')' + '.mp4'),framerate=4)
