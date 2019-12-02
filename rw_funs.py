from ris_widget import ris_widget
import pathlib

from elegant import load_data, process_data

def make_global_riswidget():
    rw_defined = 'rw' in globals()
    global rw
    if not rw_defined:
        rw = ris_widget.RisWidget()
    return rw

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

def rw_load_image_stack_fromfilenames(rw_obj, *dir_args, filter_kw='', start_idx=0,stop_idx=None):
    '''
        dir_args - list of directories to pull images from (variable length)
        filter_kw - one or more strings for filtering images based on their names
        start/stop_idx
    '''
    
    img_fns = []
    if(type(filter_kw)==str): filter_kw = [filter_kw]*len(dir_args)
    if type(dir_args[0]) is not pathlib.Path:   # Force Path objects
        dir_args = [pathlib.Path(my_dir) for my_dir in dir_args]
    
    for file_dir,kw in zip(dir_args,filter_kw):
        img_fns.append([img_file for img_file in sorted(file_dir.iterdir()) if img_file.is_file() and kw in str(img_file)])
    rw_obj.flipbook.add_image_files([collected_fns for collected_fns in zip(*img_fns)][start_idx:(stop_idx if stop_idx is not None and stop_idx < len(img_fns[0]) else len(img_fns[0]))])
    
def rw_load_last_images_fromexpt(expt_dir, channels=['bf'], rw=None):
    '''
        expt_dir - directory containing multiple animals/fields of interest
        filter_kw - Used to filter through images based on name
    '''
    
    if type(expt_dir) is not pathlib.Path: expt_dir = pathlib.Path(expt_dir)

    rw_supplied = rw is not None
    if not rw_supplied:
        rw = make_global_riswidget()

    img_fns = []
    for subdir in sorted([position_dir.parent for position_dir in expt_dir.glob('*/position_metadata.json')]):
        last_timepoint = [img_file.stem.split()[0] for img_file in subdir.glob('*.png')][-1]
        img_fns.append([img_file 
            for channel in channels 
            for img_file in subdir.glob('*.png') 
            if (channel in img_file.name and last_timepoint in img_file.name)])
    rw.flipbook.add_image_files(img_fns)
    
    if not rw_supplied:
        return rw

def rw_load_timepoint_fromexpt(rw, expt_dir, timepoint, channel='bf'):
    '''
        Loads the single bf image for a timepoint into the flipbook
        timepoint_str - timepoint in standard str format to use when searching for images
    '''
    expt_dir = pathlib.Path(expt_dir)

    '''
    def timepoint_filter(position_name, timepoint_name):
        return timepoint_name == timepoint
    position_images = load_data.scan_experiment_dir(expt_dir,timepoint_filter=timepoint_filter)
    image_names, image_paths = [], []
    for position_name, images in position_images.items():
        if images[timepoint]:   # Skip positions not having this timepoint
            image_names.append(position_name + images[timepoint][0])
            image_paths.append(images[timepoint][0])
    rw.add_image_files_to_flipbook(image_paths,page_names=image_names)

    '''
    image_files = expt_dir.glob(f'*/{timepoint} {channel}.png')
    rw.add_image_files_to_flipbook(image_files)

def rw_load_position_images(rw, expt_dir, position, channels='bf'):
    def position_filter(position_name, position_annotations, timepoint_annotations):
        return position_name == position

    process_data.update_annotations(expt_dir)
    position_images = load_data.scan_positions(expt_dir, position_filter, channels=channels)
    image_files = position_images[position].values()
    rw.add_image_files_to_flipbook(image_files)

def get_labeled_positions(rw_obj, labels):
    '''
        labels - list of labels to look for in risWidget flipbook pages (e.g. ['c'])
    '''
    
    labeled_positions = {my_label:[] for my_label in labels}
    for my_label in labels:
        for time_idx, flipbook_page in enumerate(rw_obj.flipbook.pages):
            if flipbook_page.name == my_label: labeled_positions[my_label].append(time_idx)
    return labeled_positions
