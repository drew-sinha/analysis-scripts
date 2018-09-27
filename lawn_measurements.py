import numpy
import freeimage
from skimage import feature, filters
from scipy.ndimage import morphology

import zplib.image.mask as zpl_mask
from elegant import load_data, process_images, worm_spline

class LawnMeasurements:
    feature_names = ['summed_lawn_intensity', 'median_lawn_intensity', 'background_intensity']

    def measure(self, position_root, timepoint, annotations, before, after):
        measures = {}
        timepoint_imagepath = position_root / (timepoint + ' bf.png')
        timepoint_image = freeimage.read(str(timepoint_imagepath))
        rescaled_image = process_images.pin_image_mode(timepoint_image, optocoupler=metadata['optocoupler'])

        experiment_root, position_name = position_root.parent, position_root.name
        lawn_mask = freeimage.read(experiment_root / 'derived_data' / 'lawn_masks' / f'{position_name}.png')

        # Remove the animal from the lawn if possible.
        center_tck, width_tck = annotations.get('pose', (None, None))
        if center_tck is None:
            animal_mask = numpy.zeros(timepoint_image.shape).astype('bool')   # Or should this just return Nones in the measures???
        else:
            animal_mask = worm_spline.lab_frame_mask(center_tck, width_tck, timepoint_image.shape)
        lawn_mask = lawn_mask & ~animal_mask

        vignette_mask = process_images.vignette_mask(optocoupler, timepoint_image.shape)

        measures['summed_lawn_intensity'] = numpy.sum(timepoint_image[lawn_mask])
        measures['median_lawn_intensity'] = numpy.median(timepoint_image[lawn_mask])
        measures['background_intensity'] = numpy.median(timepoint_image[~lawn_mask & vignette_mask])

##################
# Annotate the lawn information with process_data.annotate

def annotate_lawn(experiment_root, position, metadata, annotations):
    '''Position annotator used to find the lawn and associated metadata about it'''

    position_root = experiment_root / position
    lawn_mask_root = experiment_root / 'derived_data' / 'lawn_masks'

    microns_per_pixel = process_images.microns_per_pixel(metadata['objective'],metadata['optocoupler'])
    num_images_for_lawn = 10

    position_images = load_data.scan_experiment_dir(position_root.parent)[position_name]
    first_imagepaths = []
    for timepoint, timepoint_images in list(position_images.items()):
        if timepoint + ' bf.png' in timepoint_images:
            first_imagepaths.append(str(position_root / (timepoint + ' bf.png')))
        if len(first_imagepaths) == num_images_for_lawn:
            break

    first_images = list(map(freeimage.read, first_imagepaths))
    first_images = [process_images.pin_image_mode(image, optocoupler=metadata['optocoupler'])
        for image in first_images]
    median_first_images = numpy.median(first_images, axis=0)

    individual_lawns = [edge_lawn_maker(image, metadata['optocoupler']) for image in first_images]
    lawn_mask = numpy.max(individual_lawns, axis=0)
    vignette_mask = process_images.vignette_mask(metadata['optocoupler'], lawn_mask.shape)

    freeimage.write(lawn_mask, str(lawn_mask_root / f'{position}.png')) # Some better way to store this mask in the annotations?
    annotations['lawn_area'] = lawn_mask.sum() * microns_per_pixel**2


def edge_lawn_maker(image, optocoupler):
    '''Find a lawn in an image using Canny edge detection

    Parameters:
        image - numpy ndarray of the image to find the lawn from
        optocoupler - optocoupler magnification (as a float) used for the specified image

    Returns:
        lawn mask as a bool ndarray
    '''

    filtered_image = filters.median(image)  # Preliminary median filtering
    image_edges = feature.canny(image, sigma=0.02)
    vignette_mask = process_images.vignette_mask(optocoupler, image.shape)
    image_edges[~vignette_mask] = False

    image_edges = morphology.binary_dilation(image_edges, iterations = 10)
    lawn_mask = morphology.binary_fill_holes(image_edges)
    try:
        lawn_mask = zpl_mask.get_largest_object(lawn_mask)
        lawn_mask = morphology.binary_erosion(lawn_mask, iterations = 10)
    except:
        lawn_mask = numpy.zeros(lawn_mask.shape).astype('bool')

    return lawn_mask

'''
Scratch

        # initial_lawn_intensity = #FILL IN from annotations....
        # background_intensity = process_images._image_mode_numpy(rescaled_image) # Or get from original lawn during annotations
        #measures['background_diff'] = background_intensity - lawn_intensity
        #measures['relative_intensity'] = measures['background_diff']/(background_intensity - initial_lawn_intensity)

import numpy
from elegant import worm_spline, process_images
import celiagg

def circle_mask(cx, cy, r, shape):
    cx, cy, r = int(cx * shape[0]), int(cy * shape[1]), int(r * shape[0])
    path = celiagg.Path()
    path.ellipse(cx, cy, r, r)
    return worm_spline._celiagg_draw_mask(shape, path, antialias=False)
#image = rw.flipbook_pages[0][0].data
image = freeimage.read('/mnt/9karray/Sinha_Drew/20180816_spe-9_Control/087/2018-08-19t0703 bf00.png')
image = process_images.pin_image_mode(image)
mask = circle_mask(0.5,0.5,0.4,image.shape)
mask = mask > 0
#rw.flipbook_pages[0].append(mask)

image2 = image.copy()
image2[mask] = 2**16-1
#rw.flipbook_pages[0].append(image2)

nonlawn_hist = numpy.bincount(image2.flatten())
nonlawn_hist[-1] = 0
lawn_hist = numpy.bincount(image[mask].flatten())

    # annotations['initial_summed_lawn_intensity'] = numpy.sum(median_first_images[lawn_mask])    # Need to multiply by microns_per_pixel**2?
    # annotations['initial_median_lawn_intensity'] = numpy.median(median_first_images[lawn_mask])
    # annotations['initial_background_intensity'] = numpy.median(median_first_images[~lawn_mask & vignette_mask])
'''
