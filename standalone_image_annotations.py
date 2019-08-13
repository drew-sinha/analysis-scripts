import pathlib
import csv
import shutil

import numpy

from zplib.image import mask as zplib_image_mask
import freeimage
from elegant import worm_spline
from zplib.curve import spline_geometry

def load_wormsizer_results(result_directory):
    results = {}
    for result_filepath in pathlib.Path(result_directory).glob('*results.csv'):
        with result_filepath.open('r') as result_file:
            reader = csv.reader(result_file)
            header = reader.__next__()
            for line in reader:
                experiment_id = line[0]
                results.setdefault(experiment_id, {measurement:[] for measurement in header[4:]})
                is_good = line[3] == 'true'
                if is_good:
                    for measurement, measurement_value in zip(header[4:], line[4:]):
                        results[experiment_id][measurement].append(float(measurement_value))
    return results

def parse_aggregate_centerlines(aggregate_mask):
    centerlines = []
    temp_mask = aggregate_mask.copy()
    for layer in numpy.moveaxis(temp_mask, -1, 0):
        while layer.any():
            new_centerline = zplib_image_mask.get_largest_object(layer>0,structure=numpy.ones((3,3)))
            new_centerline[new_centerline>0] = -1
            layer[new_centerline>0] = 0
            centerlines.append(new_centerline)
    return centerlines

def process_centerline_dir(source_dir, microns_per_pixel):
    source_dir = pathlib.Path(source_dir)
    out_dir = source_dir / 'individual_centerlines'
    out_dir.mkdir(exist_ok=True)

    centerline_data = {}
    centerline_data_entries = ['name', 'length']
    
    mask_data = {}
    [mask_data.setdefault(entry,[]) for entry in centerline_data_entries]
    for centerline_image_path in sorted(source_dir.iterdir()):
        if centerline_image_path.suffix[1:] != 'png':
            continue
        aggregate_centerline_image = freeimage.read(centerline_image_path)
        masks = parse_aggregate_centerlines(aggregate_centerline_image)
        for mask_num, mask in enumerate(masks):
            mask_name = centerline_image_path.stem + f'_{mask_num}'
            print(mask_name)
            freeimage.write(mask.astype('uint8')*255, out_dir / (mask_name+'.png'))
            center_tck, _ = worm_spline.pose_from_mask(mask) # Toss the widths
            try:
                length = spline_geometry.arc_length(center_tck) * microns_per_pixel
                mask_data['name'].append(mask_name)
                mask_data['length'].append(length)
            except TypeError:
                print(f'Warning: couldn\'t find centerline for {mask_name} (px location {list(numpy.where(mask))})')

    with (out_dir / 'measurements.txt').open('w+') as measurement_file:
        measurement_file.write('\t'.join(centerline_data_entries)+'\n')
        for data in zip(*[mask_data[entry] for entry in centerline_data_entries]):
            measurement_file.write('\t'.join([str(item) for item in data])+'\n')

def copy_experiment_images(experiment_dir, timepoint):
    '''Move brightfield images for a particular experiment + timepoint; 
        useful for collecting images for separate annotating (e.g. centerlines)
    '''
    experiment_dir = pathlib.Path(experiment_dir)
    (experiment_dir / 'derived_data' / timepoint).mkdir(exist_ok=True,parents=True)

    for position_dir in [pmd_file.parent for pmd_file in sorted(experiment_dir.glob('*/position_metadata.json'))]:
        try:
            shutil.copyfile(position_dir / f'{timepoint} bf.png', experiment_dir / 'derived_data' / timepoint / f'{position_dir.name} {timepoint} bf.png')
        except FileNotFoundError:
            print(f'No timepoint image found for {position_dir.name} at {timepoint}')

