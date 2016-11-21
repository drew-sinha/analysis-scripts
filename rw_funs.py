import ris_widget
import pathlib


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
    
def rw_load_last_images_fromexpt(rw_obj, expt_dir, filter_kw='bf.png'):
    '''
        expt_dir - directory containing multiple animals/fields of interest
        filter_kw - Used to filter through images based on name
    '''
    
    if type(expt_dir) is not pathlib.Path: expt_dir = pathlib.Path(expt_dir)
    
    img_fns = []
    for subdir in expt_dir.iterdir():
        if subdir.is_dir() and subdir.parts[-1].isnumeric():
            img_fns.append([img_file for img_file in subdir.iterdir() if img_file.is_file() and filter_kw in str(img_file)][-1])
    rw_obj.flipbook.add_image_files(img_fns)

def get_labeled_positions(rw_obj, labels):
    '''
        labels - list of labels to look for in risWidget flipbook pages (e.g. ['c'])
    '''
    
    labeled_positions = {my_label:[] for my_label in labels}
    for my_label in labels:
        for time_idx, flipbook_page in enumerate(rw_obj.flipbook.pages):
            if flipbook_page.name == my_label: labeled_positions[my_label].append(time_idx)
    return labeled_positions
