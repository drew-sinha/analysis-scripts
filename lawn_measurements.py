import json
import pickle

import numpy
import freeimage
from skimage import feature, filters
from scipy.ndimage import morphology
import scipy.ndimage.filters as ndi_filters
from sklearn import mixture

import zplib.image.mask as zpl_mask
from elegant import load_data, process_images, worm_spline

class LawnMeasurements:
    feature_names = ['summed_lawn_intensity', 'median_lawn_intensity', 'background_intensity']

    def measure(self, position_root, timepoint, annotations, before, after):
        measures = {}
        print(f'Working on position {position_root.name} - {timepoint}')

        # Load metadata.... TODO ask Zach if this can be a feature of the measurement signature
        experiment_root, position_name = position_root.parent, position_root.name
        with (experiment_root / 'experiment_metadata.json').open('r') as md_file:
            metadata = json.load(md_file)

        timepoint_imagepath = position_root / (timepoint + ' bf.png')
        timepoint_image = freeimage.read(timepoint_imagepath)
        rescaled_image = process_images.pin_image_mode(timepoint_image, optocoupler=metadata['optocoupler'])

        lawn_mask = freeimage.read(experiment_root / 'derived_data' / 'lawn_masks' / f'{position_name}.png').astype('bool')
        vignette_mask = process_images.vignette_mask(metadata['optocoupler'], timepoint_image.shape)

        # Remove the animal from the lawn if possible.
        center_tck, width_tck = annotations.get('pose', (None, None))
        if center_tck is None:
            animal_mask = numpy.zeros(timepoint_image.shape).astype('bool')   # Or should this just return Nones in the measures???
        else:
            animal_mask = worm_spline.lab_frame_mask(center_tck, width_tck, timepoint_image.shape).astype('bool')
        lawn_mask = lawn_mask & ~animal_mask

        measures['summed_lawn_intensity'] = numpy.sum(rescaled_image[lawn_mask])
        measures['median_lawn_intensity'] = numpy.median(rescaled_image[lawn_mask])
        measures['background_intensity'] = numpy.median(rescaled_image[~lawn_mask & vignette_mask])

        return [measures[feature_name] for feature_name in self.feature_names]

##################
# Annotate the lawn information with process_data.annotate

def annotate_lawn(experiment_root, position, metadata, annotations):
    '''Position annotator used to find the lawn and associated metadata about it'''
    num_images_for_lawn = 3

    print(f'Working on position {position}')
    position_root = experiment_root / position
    lawn_mask_root = experiment_root / 'derived_data' / 'lawn_masks'
    lawn_mask_root.mkdir(parents=True, exist_ok=True)

    lawn_model_root = experiment_root / 'derived_data' / 'lawn_models'
    lawn_model_root.mkdir(parents=True, exist_ok=True)

    microns_per_pixel = process_images.microns_per_pixel(metadata['objective'],metadata['optocoupler'])

    position_images = load_data.scan_experiment_dir(experiment_root)[position]
    first_imagepaths = []
    for timepoint, timepoint_images in position_images.items():
        if position_root / (timepoint + ' bf.png') in timepoint_images:
            first_imagepaths.append(str(position_root / (timepoint + ' bf.png')))
        if len(first_imagepaths) == num_images_for_lawn:
            break

    first_images = list(map(freeimage.read, first_imagepaths))
    first_images = [process_images.pin_image_mode(image, optocoupler=metadata['optocoupler'])
        for image in first_images]
    median_first_images = numpy.median(first_images, axis=0)

    individual_lawns = [gmm_lawn_maker(image, metadata['optocoupler']) for image in first_images]
    lawn_mask = numpy.max(individual_lawns, axis=0)

    median_lm, gmm_model = gmm_lawn_maker(median_first_images.astype('int'), metadata['optocoupler'], return_model=True)

    vignette_mask = process_images.vignette_mask(metadata['optocoupler'], lawn_mask.shape)

    with (lawn_model_root / f'{position}.pickle').open('wb') as lm_file:
        gmm_means = gmm_model.means_.flatten()
        gmm_stds = numpy.sqrt(gmm_model.covariances_.flatten())
        sorting_idxs = gmm_means.argsort()
        gmm_means, gmm_stds = gmm_means[sorting_idxs], gmm_stds[sorting_idxs] # First one is lawn (since less intensity is darker)

        model_features = {'lawn_mean':gmm_means[0], 'background_mean':gmm_means[1], 
            'lawn_std': gmm_stds[0], 'background_std': gmm_stds[1]}
        pickle.dump(model_features, lm_file)
    freeimage.write(lawn_mask.astype('uint8')*255, str(lawn_mask_root / f'{position}.png')) # Some better way to store this mask in the annotations?
    annotations['lawn_area'] = lawn_mask.sum() * microns_per_pixel**2

def gmm_lawn_maker(image, optocoupler, return_model=False):
    '''Find a lawn in an image use Gaussian mixture modeling (GMM)

    This lawn maker models an image (i.e. its pixel intensities) as as mixture
        of two Gaussian densities. Each corresponds to either the background & lawn.

    Parameters:
        image - numpy ndarray of the image to find the lawn from
        optocoupler - optocoupler magnification (as a float) used for the specified image

    Returns:
        lawn mask as a bool ndarray
        fitted GMM model
    '''
    scaled_image = ndi_filters.median_filter(image, size=(3,3), mode='constant')
    vignette_mask = process_images.vignette_mask(optocoupler, image.shape)

    img_data = scaled_image[vignette_mask]
    img_hist = numpy.bincount(img_data.flatten())
    img_hist = img_hist/img_hist.sum()

    gmm = mixture.GaussianMixture(n_components=2)
    gmm.fit(numpy.expand_dims(img_data,1))

    # Calculate boundary point for label classification as intensity threshold
    gmm_support = numpy.linspace(0,2**16-1,2**16)

    labels = gmm.predict(numpy.reshape(gmm_support, (-1,1)))
    thr = numpy.argmax(numpy.abs(numpy.diff(labels)))
    lawn_mask = (scaled_image < thr) & vignette_mask

    lawn_mask = morphology.binary_erosion(lawn_mask, iterations=10)
    lawn_mask = zpl_mask.get_largest_object(lawn_mask)
    lawn_mask = morphology.binary_fill_holes(lawn_mask)
    lawn_mask = morphology.binary_dilation(lawn_mask, iterations=10)

    if return_model:
        return lawn_mask, gmm
    else:
        return lawn_mask

if __name__ == "__main__":
    import sys
    from elegant import process_data
    import measurement_pipeline

    expt_dir = sys.argv[1]

    def make_lawn_measurements(experiment_root):
        measures = [LawnMeasurements()] #, process_data.PoseMeasurements(microns_per_pixel=5)]
        measurement_name = 'lawn_measures'

        #process_data.update_annotations(experiment_root)
        measurement_pipeline.propagate_stages(experiment_root,verbose=True)
        positions = load_data.read_annotations(experiment_root)
        # to_measure = load_data.filter_annotations(positions, load_data.filter_living_timepoints)
        to_measure = load_data.filter_annotations(positions, measurement_pipeline.filter_living_timepoints)
        process_data.measure_worms(experiment_root, to_measure, measures, measurement_name)

    process_data.annotate(expt_dir, position_annotators=[annotate_lawn])
    make_lawn_measurements(expt_dir)
    process_data.collate_data(expt_dir)
