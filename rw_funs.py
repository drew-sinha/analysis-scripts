import ris_widget
import os
import freeimage


def rw_load_image_stack_fromfilenames(rw_obj, *dir_args, filter_kw='', start_idx=0,stop_idx=None):
    img_fns = []
    if(type(filter_kw)==str): filter_kw = [filter_kw]*len(dir_args)
    print(dir_args)
    for file_dir,kw in zip(dir_args,filter_kw):
        print(file_dir)
        print(kw)
        img_fns.append([file_dir+os.path.sep+file_n for file_n in sorted(os.listdir(file_dir)) if os.path.isfile(file_dir+os.path.sep+file_n) and kw in file_n])
    print(img_fns)
    print(len(img_fns[0]))
    rw_obj.flipbook.add_image_files([collected_fns for collected_fns in zip(*img_fns)][start_idx:(stop_idx if stop_idx is not None and stop_idx < len(img_fns[0]) else len(img_fns[0]))])
    #rw_obj.flipbook.add_image_files([[freeimage.read(img_f) for img_f in collected_fns] for collected_fns in zip(*img_fns)])
    
def rw_load_last_images_fromexpt(rw_obj, expt_dir, filter_kw='bf.png'):
    '''
        expt_dir - make this a Path object
    '''
    
    img_fns = []
    for subdir in expt_dir.iterdir():
        if subdir.is_dir() and subdir.parts[-1].isnumeric():
            img_fns.append([img_file for img_file in subdir.iterdir() if img_file.is_file() and filter_kw in str(img_file)][-1])
    rw_obj.flipbook.add_image_files(img_fns)

def get_labeled_positions(rw_obj, labels = ['c']):
    labeled_positions = {my_label:[] for my_label in label}
    for my_label in labels:
        for time_idx, flipbook_page in enumerate(rw.flipbook.pages):
            if flipbook_page.name == my_label: labeled_positions['my_label'].append(times)
    return labeled_positions
