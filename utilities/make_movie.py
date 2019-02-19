import pathlib
import numpy as np

from zplib.image import write_movie
import freeimage

def yield_rgb(image_generator):
    for image in image_generator:
        if image.ndim == 2:
            image = image[:,:,np.newaxis]
        yield np.repeat(image, 3, axis=2)

def shrink(image_generator, factor=2, fast=False):
    from zplib.image import pyramid
    """Shrink images produced by a generator by the specified factor.
    Parameters:
        factor: amount to shrink the image by (fold-change)
        fast: if True and if factor is an integer, perform no smoothing.
            If False, smooth the image before downsampling to avoid aliasing.
    """
    # can only do fast downsampling if fast was requested AND the factor is integral
    fast = fast and int(factor) == factor
    for image in image_generator:
        if fast:
            if image.ndims == 3:
                yield image[::int(factor), ::int(factor),:]
            else:
                yield image[::int(factor), ::int(factor)]
        else:
            # TODO handle 3D images/ change color tint to properly handle 3d images
            yield pyramid.pyr_down(image, factor).astype(np.uint8)

'''
    # TODO: Why are these above methods needed again? (I know there was a reason!)
'''

def double_image_layout(image_generator1, image_generator2):
    for image1, image2 in zip(image_generator1, image_generator2):
        assert (image1.shape == image2.shape).all()
        combined_image = np.zeros((image1.shape[0], 2*image1.shape[1]))
        combined_image[:image1.shape[0]] = image1
        combined_image[image1.shape[0]:] = image2
        yield combined_image


def make_experiment_movies(experiment_root, output_dir,
    shrink_factor=4, framerate=5,glob_str='*bf.png',num_frames=-1,positions=None, **scale_params):
    experiment_root = pathlib.Path(experiment_root)
    if positions is None:
        position_roots = sorted(p.parent for p in experiment_root.glob('*/position_metadata.json'))
    else:
        position_roots = experiment_root / position for position in positions
        for position_root in position_roots:
            assert position_root.exists()

    for position_root in position_roots:
        image_paths = sorted(position_root.glob(glob_str))
        if num_frames != -1: image_paths = image_paths[:num_frames]

        image_generator = write_movie.generate_images_from_files(image_paths,**scale_params)
        image_generator = write_movie.shrink(image_generator, factor=shrink_factor, fast=True)
        write_movie.write_movie(image_generator, output_file, framerate=framerate)


# image_dir = pathlib.Path('/mnt/9karray/Sinha_Drew/20180810_spe-9_Control/012/')
# image_filepaths = sorted(image_dir.glob('*bf.png'))
# image_generator = write_movie.generate_images_from_files(image_filepaths)
# image_generator = write_movie.shrink(yield_rgb(image_generator), factor=4) #pyramid.pyr_down returning [x,y,1] shaped image leads to an assertion error.... Should the assertion be around?
# write_movie.write_movie(image_generator, 'test.mp4',framerate=5)
