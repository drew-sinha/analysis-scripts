import pathlib
import numpy as np
import freeimage
import json
import zplib.util
from ris_widget.ris_widget import RisWidget
from ris_widget.qwidgets.flipbook import ImageList, Image
from ris_widget.point_list_picker import PointListPicker
#from ris_widget.examples import SimplePointPicker
# from simple_point_picker import SimplePointPicker
from ris_widget.overlay import point_set

from skimage.morphology import binary_erosion,label
from skimage.measure import find_contours
import skimage.filters

from zplib.util import json_encode_legible_to_file

import matplotlib.pyplot as plt

# Work around for promiscuous PyQt/IPython
import IPython
import prompt_toolkit

import time

def ip_input(message=''):
    ip = IPython.get_ipython()
    el = prompt_toolkit.shortcuts.create_eventloop(ip.inputhook)
    return prompt_toolkit.prompt(message, eventloop=el)

'''
def add_event_to_metadata(expt_dir,event):
    ''
        event - dictionary with items to append to corresponding entry in metadata
            key - entry title in metadata
            value - to add
    ''
    if type(expt_dir) is not pathlib.Path: expt_dir = pathlib.Path(expt_dir)
    with (expt_dir/'experiment_metadata.json').open('r') as mdata_file:
        metadata = json.load(mdata_file)
    event_entry = metadata.get(event['key'],[]).append(event['value']) # Workaround if a key yet
    metadata[event['key']] = event_entry
    with (expt_dir/'experiment_metadata.json').open('w') as mdata_file:
        json_encode_legible_to_file(metadata, mdata_file)
'''

def take_img(scope,expt_dir,poll=True):
    if poll:
        try:
            ip_input('Press enter when position acquired; press ctrl-c to abort.')
        except KeyboardInterrupt:
            raise
    position = scope.stage.position

    scope.camera.acquisition_sequencer.new_sequence()
    scope.camera.acquisition_sequencer.add_step(2,'TL')
    my_image = scope.camera.acquisition_sequencer.run()

    return [position,my_image[0]]

def make_align_img(scope,expt_dir):
    expt_dir = pathlib.Path(expt_dir)
    try:
        scope_pos, my_image = take_img(scope,expt_dir)
    except KeyboardInterrupt:
        return

    time_label = time.strftime('%Y%m%d-%H%M-%S')

    with (expt_dir/'experiment_metadata.json').open('r') as mdata_file:
        metadata = json.load(mdata_file)
    with (expt_dir/f'experiment_metadata_noalign_{time_label}.json').open('w') as mdata_file:
        json_encode_legible_to_file(metadata,mdata_file)

    metadata['align_position'] = scope_pos
    with (expt_dir/'experiment_metadata.json').open('w') as mdata_file:
        json_encode_legible_to_file(metadata, mdata_file)
    freeimage.write(my_image,
        expt_dir/'calibrations'/'align_image.png')


def perform_alignment(scope,expt_dir,rw):
    px_conversion = 1.304/1000   #mm/px
    expt_dir = pathlib.Path(expt_dir)

    before_img = freeimage.read(str(expt_dir/'calibrations'/'align_image.png'))
    with (expt_dir/'experiment_metadata.json').open('r') as mdata_file:
        metadata = json.load(mdata_file)
    before_scope_pos = np.array(metadata['align_position'])
    if scope is not None:
        scope.stage.position = before_scope_pos   # Set z to same position
        after_scope_pos, after_img = take_img(scope,expt_dir)
    else:   # old debugging not on scope
        after_scope_pos = np.array([163.6972, 3.085, 23.739])
        #~ after_img = freeimage.read('/mnt/scopearray/Sinha_Drew/testing/align_pics/after_img.png')
        #after_img = freeimage.read('/media/Data/Documents/Dropbox/align_pics/after_img.png')
        #after_img = freeimage.read('/home/zplab/Desktop/align_pics/after_img.png')

    if len(rw.flipbook.pages)>0: rw.flipbook.pages.clear()
    rw.flipbook.pages.append(ImageList([Image(before_img),Image(after_img)]))
    #rw.layers[0].tint=[255,232,185,0.5]
    #rw.layers[1].tint = [270,70,255,0.5]
    # my_ptpicker = PointListPicker(rw.main_view,rw.main_scene.layer_stack_item)
    # my_ptpicker = SimplePointPicker(rw.main_view,rw.main_scene.layer_stack_item)
    my_ptpicker = point_set.PointSet(rw)

    ip_input('Click on first landmark in old and new images; press Enter when done.')
    before_first_pos, after_first_pos = [np.array(point) for point in my_ptpicker.points]   # IN PIXELS
    print('Acquired first landmark')
    # my_ptpicker.points = []
    for point in my_ptpicker.points: point.remove()

    ip_input('Click on second landmark in old and new images; press Enter when done.')
    before_next_pos, after_next_pos = [np.array(point) for point in my_ptpicker.points]
    print('Acquired second landmark')
    # my_ptpicker.points = []
    for point in my_ptpicker.points: point.remove()


    xy_offset = after_scope_pos[:2]+after_first_pos*px_conversion - (before_scope_pos[:2]+before_first_pos*px_conversion)
    xy_offset = xy_offset*[-1,1]    #Correct for directionality of stage on ISCOPE    CHECK THIS ON SCOPE!!!!!!!!!!!
    print(xy_offset)

    before_next_adj = before_next_pos-before_first_pos
    after_next_adj = after_next_pos-after_first_pos

    before_incline = np.arctan(before_next_adj[1]/before_next_adj[0])
    after_incline = np.arctan(after_next_adj[1]/after_next_adj[0])

    incline_offset = after_incline-before_incline
    print('Before incline:{:.2f}\nAfter incline:{:.2f}\nOffset:{:.2f}'.format(before_incline,after_incline,incline_offset)) # GOOD UP TO HERE.

    #raise Exception("xy:{} (px), incline:{} (deg)".format(xy_offset,incline_offset*360/2*np.pi))


    old_metadata = metadata.copy()
    metadata['align_position'].append(after_scope_pos)  # Save the old!
    for well in metadata['positions']:
        # Do rotation first
        rotated_point = np.dot(np.array([[np.cos(incline_offset),-np.sin(incline_offset)],[np.sin(incline_offset),np.cos(incline_offset)]]),
            np.array(metadata['positions'][well][:2]-np.array(before_scope_pos[:2]+before_first_pos*px_conversion)))+np.array(before_scope_pos[:2]+before_first_pos*px_conversion)
        metadata['positions'][well][:2]=list(rotated_point+xy_offset)

    wells = [well_dir.parts[-1] for well_dir in expt_dir.iterdir() if str(well_dir.parts[-1]).isnumeric()]
    metadata['realignment_time']=metadata.get('realignment_time',[]).append(time.strftime("%Y-%m-%dt%H%M"))

    # Compare before and after positions (offset only)
    scope.stage.position = list(np.append(before_scope_pos[:2]+xy_offset,before_scope_pos[2]))
    trash,final_landmark_img = take_img(scope,expt_dir,poll=False)
    rw.flipbook.pages.append(ImageList([Image(before_img),Image(final_landmark_img)]))
    ip_input('Press anything if landmarks looks aligned (xy-offset); ctrl-c to abort')

    # Compare first and last wells
    rw.flipbook.pages.clear()
    before_first_img =freeimage.read([bf_file for bf_file in (expt_dir/wells[0]).iterdir() if 'bf.png' in bf_file.parts[-1]][0])
    scope.stage.position = list(np.append(np.array(old_metadata['positions'][wells[0]][:2])+xy_offset,old_metadata['positions'][wells[0]][2]))  # Do one with the offset for debugging
    trash, after_first_offset_img = take_img(scope,expt_dir,poll=False)
    rw.flipbook.pages.append(ImageList([Image(before_first_img),Image(after_first_offset_img)]))    #,name='Before After Well'+str(wells[0]))
    scope.stage.position = metadata['positions'][wells[0]]
    trash, after_first_img = take_img(scope,expt_dir,poll=False)
    rw.flipbook.pages.append(ImageList([Image(before_first_img),Image(after_first_img)]))    #,name='Before After Well'+str(wells[0]))
    ip_input('Press anything if first well looks aligned; ctrl-c to abort')

    before_last_img =freeimage.read([bf_file for bf_file in (expt_dir/wells[-1]).iterdir() if 'bf.png' in bf_file.parts[-1]][0])
    scope.stage.position = metadata['positions'][wells[-1]]
    trash, after_last_img = take_img(scope,expt_dir,poll=False)
    rw.flipbook.pages.clear()
    rw.flipbook.pages.append(ImageList([Image(before_last_img),Image(after_last_img)]))  # ,name='Before After Well'+str(wells[-1]))
    ip_input('Press anything if last well looks aligned; ctrl-c to abort')

    with (expt_dir/'experiment_metadata_oldprealignment.json').open('w') as mdata_file:
        json_encode_legible_to_file(old_metadata, mdata_file)
    with (expt_dir/'experiment_metadata.json').open('w') as mdata_file:
        json_encode_legible_to_file(metadata, mdata_file)
    freeimage.write(before_img,
        str(expt_dir/'calibrations'/'old_align_image.png'))
    freeimage.write(after_img,
        str(expt_dir/'calibrations'/'align_image.png'))
    #~ except KeyboardInterrupt:
        #~ return

# def auto_alignment(expt_dir,rw=None):
#     if type(expt_dir) is not pathlib.Path: expt_dir = pathlib.Path(expt_dir)
#
#     def get_image_border(img):
#         img_mask = img > 1000
#         contours = find_contours(skimage.filters.median(img,np.ones((3,3))),1000)
#         #~ prelim_border = img_mask & (~binary_erosion(img_mask))
#         #~ labeled_border,num_labels = label(prelim_border,return_num=True)
#         #~ label_counts = np.array([np.count_nonzero(labeled_border==num) if num !=0 else 0 for num in range(num_labels)])
#         #~ return labeled_border == np.argmax(label_counts)
#         return contours
#
#     px_distance = lambda x,y: np.sqrt(x**2+y**2)
#
#     before_img = freeimage.read(str(expt_dir/'calibrations'/'align_image.png'))
#     with (expt_dir/'experiment_metadata.json').open('r') as mdata_file:
#         metadata = json.load(mdata_file)
#     before_scope_pos = np.array(metadata['align_position'])
#     #if scope is not None: scope.stage.position[2] = before_scope_pos[2]   # Set z to same position
#     #after_scope_pos, after_img = take_img(scope,expt_dir)
#     after_scope_pos = np.array([163.6972, 3.085, 23.739])
#     #~ after_img = freeimage.read('/mnt/scopearray/Sinha_Drew/testing/align_pics/after_img.png')
#     after_img = freeimage.read('/media/Data/Documents/Dropbox/align_pics/after_img.png')
#     fig_h, ax_h = plt.subplots(1,2,sharex=True,sharey=True)
#     ax_h[0].imshow(before_img)
#     ax_h[1].imshow(after_img)
#
#     # Get image borders
#     before_border = get_image_border(before_img)
#     after_border = get_image_border(after_img)
#     before_border_pos = np.array(np.where(before_border))
#     after_border_pos = np.array(np.where(after_border))
#     fig_h, ax_h = plt.subplots(1,2,sharex=True,sharey=True)
#     ax_h[0].imshow(before_border)
#     ax_h[1].imshow(after_border)
#
#     before_corner_idx = np.argmin(px_distance(*before_border_pos))
#     after_corner_idx = np.argmin(px_distance(*after_border_pos))
#
#     print(before_border_pos)
#     print(after_border_pos)
#
#     xy_offset = after_border_pos[:,after_corner_idx]-before_border_pos[:,before_corner_idx]
#
#
#     ax_h[0].scatter(*reversed(before_border_pos[:,before_corner_idx]))
#     ax_h[1].scatter(*reversed(after_border_pos[:,after_corner_idx]))
#     print(xy_offset)

# def add_offset(expt_dir,offset,scope,rw):
#     '''
#         offset - list of three numbers (e.g. [0,0,0])
#     '''
#     if type(expt_dir) is not pathlib.Path: expt_dir = pathlib.Path(expt_dir)
#     with (expt_dir/'experiment_metadata.json').open('r') as mdata_file:
#         metadata = json.load(mdata_file)
#     old_metadata = metadata.copy()
#     for well in metadata['positions']:
#         metadata['positions'][well]=list(np.array(metadata['positions'][well])+np.array(offset))
#     wells = [well_dir.parts[-1] for well_dir in expt_dir.iterdir() if str(well_dir.parts[-1]).isnumeric()]
#
#
#     if len(rw.flipbook.pages)>0: rw.flipbook.pages.clear()
#     before_first_img =freeimage.read([bf_file for bf_file in (expt_dir/wells[0]).iterdir() if 'bf.png' in bf_file.parts[-1]][0])
#
#     scope.stage.position = metadata['positions'][wells[0]]  #TODO Actually grab first well with padding
#     trash, after_first_offset_img = take_img(scope,expt_dir,poll=False)
#     rw.flipbook.pages.append(ImageList([Image(after_first_offset_img)]))
#     ip_input('Press anything if the well looks good; ctrl-c to abort')
#
#     with (expt_dir/'experiment_metadata_oldpreoffset.json').open('w') as mdata_file:
#         json_encode_legible_to_file(old_metadata, mdata_file)
#     with (expt_dir/'experiment_metadata.json').open('w') as mdata_file:
#         json_encode_legible_to_file(metadata, mdata_file)

# if __name__ is "__main__":
#     rw = RisWidget()
#     expt_dir = '/media/Data/Documents/Dropbox/align_pics/20170307_alignment/'
#     perform_alignment(None, expt_dir,rw)


'''
    rw = RisWidget()
    expt_dir = '/media/Data/Documents/Dropbox/align_pics/20170307_alignment/'
    scope_alignment.perform_alignment(None, expt_dir,rw)
'''
