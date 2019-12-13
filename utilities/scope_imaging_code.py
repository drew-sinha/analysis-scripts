import os
import freeimage
import pathlib

import numpy as np

from zplib import datafile

lamp_dict = {'cyan':'gfp','green_yellow':'RedmChr', 'teal':'yfp'}

def get_image_sequence_simple(scope, position_data, out_dir, lamp=None):
    '''
        scope - ScopeClient object
        position_data - positions acquired using get_objpositions
        out_dir - String path
        lamp - List of the form (exposure time, lamp_name)
    '''

    scope.camera.live_mode=False
    scope.camera.acquisition_sequencer.new_sequence()
    scope.camera.acquisition_sequencer.add_step(2,'TL', tl_intensity=255)
    if lamp is not None: scope.camera.acquisition_sequencer.add_step(lamp[0],lamp[1])
    if not os.path.isdir(out_dir): os.mkdir(out_dir)
    for pos_num, this_position_data in enumerate(position_data):
        scope.nosepiece.position = this_position_data[0]
        scope.stage.position = this_position_data[1:]
        my_images = scope.camera.acquisition_sequencer.run()
        freeimage.write(my_images[0], out_dir+os.path.sep+'_{:03d}_bf.png'.format(pos_num))
        if lamp is not None: freeimage.write(my_images[1], out_dir+os.path.sep+'_{:03d}_'.format(pos_num)+lamp_dict[lamp[1]]+'.png')

def get_image_sequence(scope, position_data, out_dir, lamp=None):
    '''
        scope - ScopeClient object
        position_data - positions acquired using get_objpositions
        out_dir - String path
        lamp - List of lists of the form [[lamp_exposure1,lamp_name1],...]
    '''

    if lamp is None:
        lamp = [[2, 'TL']]
    if not any([arg[1] is 'TL' for arg in lamp]):
        lamp.append([2,'TL'])


    scope.camera.live_mode=False
    scope.camera.acquisition_sequencer.new_sequence()
    for lamp_exposure, lamp_name in lamp:
        if lamp_name is not 'TL': scope.camera.acquisition_sequencer.add_step(lamp_exposure, lamp_name)
        else: scope.camera.acquisition_sequencer.add_step(lamp_exposure, lamp_name, tl_intensity=255)
    if not os.path.isdir(out_dir): os.mkdir(out_dir)
    for pos_num, this_position_data in enumerate(position_data):
        scope.nosepiece.position = this_position_data[0]
        scope.stage.position = this_position_data[1:]
        my_images = scope.camera.acquisition_sequencer.run()
        for (lamp_exposure, lamp_name, this_image) in zip([arg[0] for arg in lamp],[arg[1] for arg in lamp], my_images):
            if lamp_name is 'TL': freeimage.write(this_image, out_dir+os.path.sep+'_{:03d}_bf_{}_ms'.format(pos_num,lamp_exposure)+'.png')
            else: freeimage.write(this_image, out_dir+os.path.sep+'_{:03d}_'.format(pos_num)+lamp_dict[lamp_name]+'_{}_ms'.format(lamp_exposure)+'.png')

def get_objpositions(scope):
    """Return a list of interactively-obtained scope stage positions."""
    positions = []
    print('Press enter after each position has been found; press control-c to end')
    while True:
        try:
            input()
        except KeyboardInterrupt:
            break
        positions.append(scope.stage.position)
        positions[-1].insert(0,scope.nosepiece.position)
        print('Position {}: {}'.format(len(positions), tuple(positions[-1])), end='')
    return positions

def take_a_picture(scope,position_data,out_dir, lamp=None):
    '''
        Take a picture ignoring info. on objective

        scope - ScopeClient object
        position_data - positions acquired using get_objpositions
        out_dir - String path
        lamp - List of the form (exposure time, lamp_name)
    '''

    scope.camera.live_mode=False
    scope.camera.acquisition_sequencer.new_sequence()
    for lamp_exposure, lamp_name in lamp:
        if lamp_name is not 'TL': scope.camera.acquisition_sequencer.add_step(lamp_exposure, lamp_name)
        else: scope.camera.acquisition_sequencer.add_step(lamp_exposure, lamp_name, tl_intensity=255)
    if not os.path.isdir(out_dir): os.mkdir(out_dir)
    for pos_num, this_position_data in enumerate(position_data):
        #scope.nosepiece.position = this_position_data[0]
        scope.stage.position = this_position_data[1:]
        my_images = scope.camera.acquisition_sequencer.run()
        for (lamp_exposure, lamp_name, this_image) in zip([arg[0] for arg in lamp],[arg[1] for arg in lamp], my_images):
            if lamp_name is 'TL': freeimage.write(this_image, out_dir+os.path.sep+'_{:03d}_bf_{}_ms'.format(pos_num,lamp_exposure)+'.png')
            else: freeimage.write(this_image, out_dir+os.path.sep+'_{:03d}_'.format(pos_num)+lamp_dict[lamp_name]+'_{}_ms'.format(lamp_exposure)+'.png')

def take_sequential_images(scope, out_dir, tl_intensity):
    '''
        Take sequential images as one specifies positions on a stage
    '''

    out_dir = pathlib.Path(out_dir)
    scope.camera.acquisition_sequencer.new_sequence()
    scope.camera.acquisition_sequencer.add_step(2, 'TL', tl_intensity=tl_intensity)
    out_dir.mkdir(exist_ok=True)
    pos_num = 0
    print('Press enter after each position has been found; press control-c to end')

    while True:
        try:
            input()
        except KeyboardInterrupt:
            break
        bf_image = scope.camera.acquisition_sequencer.run()[0]
        freeimage.write(bf_image, out_dir / f'_{pos_num:03d}.png')
        pos_num += 1

    if pos_num > 0:
        imaging_parameters = {'lamp':'TL', 'exposure':2, 'intensity':tl_intensity}
        with (out_dir / 'imaging_parameters.json').open('w') as param_file:
            datafile.json_encode_legible_to_file(imaging_parameters, param_file)

def well_plate_acquisition(scope, out_dir, tl_intensity, grid_size = [4,6], well_spacing_cc=12.92):
    '''
        grid_size, well_spacing_cc for nominal Falcon 48-well plate (latter in mm)

    '''


    out_dir = pathlib.Path(out_dir)
    scope.camera.acquisition_sequencer.new_sequence()
    scope.camera.acquisition_sequencer.add_step(2, 'TL', tl_intensity=tl_intensity)
    out_dir.mkdir()

    try:
        print('Specify xy-position for top-left well; press ctrl-c to abort.')
        input()
        topleft_position = scope.stage.position

        pos_num = 0
        for row_num in range(grid_size[0]):
            for column_num in range(grid_size[1]):
                scope.stage.position = [
                    topleft_position[0]+column_num*well_spacing_cc,
                    topleft_position[1]+row_num*well_spacing_cc,
                    topleft_position[2],
                ]
                print('Adjust to desired z-position; press ctrl-c to abort.')
                input()

                bf_image = scope.camera.acquisition_sequencer.run()[0]
                freeimage.write(bf_image, out_dir / f'_{pos_num:03d}.png')
                pos_num += 1
        imaging_parameters = {'lamp':'TL', 'exposure':2, 'intensity':tl_intensity}
        with (out_dir / 'imaging_parameters.json').open('w') as param_file:
            datafile.json_encode_legible_to_file(imaging_parameters, param_file)
    except KeyboardInterrupt:
        return

def take_automated_plate_images(scope, out_dir):
    out_dir = pathlib.Path(out_dir)
    scope.camera.acquisition_sequencer.new_sequence()
    scope.camera.acquisition_sequencer.add_step(2, 'TL', tl_intensity=255)
    out_dir.mkdir(exist_ok=True)
    field_spacing = 2160 * 0.0065 * 1 / 2.5 # FILL ME IN WITH APPROPRIATE FIELD SIZE BASED ON 2.5X OBJECTIVE
        
    try:
        input('Specify center of plate')
        center_position = scope.stage.position
        input('Specify outer extent of plate')
        outer_position = scope.stage.position
        
        roi_radius = ((center_position[0] - outer_position[0])**2 + (center_position[1] - outer_position[1])**2)**0.5
        
        # Define function to interpolate z - assume the plate surface is a parabola with radial symmetry about its center w.r.t. both x & y
        scale_param = (outer_position[2] - center_position[2]) / (roi_radius**2)
        interpolate_z = lambda x, y: scale_param * ((x - center_position[0])**2 + (y - center_position[1])**2) + center_position[2]
        
        grid_x = np.arange(center_position[0] - roi_radius/np.sqrt(2), center_position[0] + roi_radius/np.sqrt(2), field_spacing)
        grid_y = np.arange(center_position[1] - roi_radius/np.sqrt(2), center_position[1] + roi_radius/np.sqrt(2), field_spacing)
        xcoor, ycoor = np.meshgrid(grid_x, grid_y)
        xcoor = np.array([pts if num % 2 else pts[::-1] for num, pts in enumerate(xcoor)]) # invert the x-coordinates appropriately so that we make take min time to traverse the slide
        #raise Exception()
        
        pos_num = 0
        for x, y in zip(xcoor.flatten(), ycoor.flatten()):
            scope.stage.position = [x, y, interpolate_z(x,y)]
            
            bf_image = scope.camera.acquisition_sequencer.run()[0]
            freeimage.write(bf_image, out_dir / f'_{pos_num:03d}.png')
            pos_num += 1
            
        imaging_parameters = {'lamp':'TL', 'exposure':2, 'intensity':255}
        with (out_dir / 'imaging_parameters.json').open('w') as param_file:
            datafile.json_encode_legible_to_file(imaging_parameters, param_file)
        scope.stage.position = center_position
    except KeyboardInterrupt:
        return
