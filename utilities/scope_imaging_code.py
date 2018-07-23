import os
import freeimage

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
