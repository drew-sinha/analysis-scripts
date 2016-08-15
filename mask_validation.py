'''
    Tools for working with masks
'''
from collections import OrderedDict
import json
import os
import freeimage

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

# Load masks from each directory GOOD
# Compare individual masks to each other and calc summary statistics per image

def find_char(string,ch):
    '''
        Pulls out positions for a particular character in a string
    '''
    return [idx for (idx, letter) in enumerate(string) if letter == ch]

def load_images_fromdir(image_dir,filter_kw = '', master_vig = None):
    '''
        Return an ordered dictionary with key as file_name
    '''
    image_list = []
    #image_fn_list = [image_fn[:find_char(image_fn,'.')[-1]] for image_fn in sorted(os.listdir(image_dir)) \
    #    if os.path.isfile(image_dir+os.path.sep+image_fn) and filter_kw in image_fn]
    image_fn_list = [image_fn for image_fn in sorted(os.listdir(image_dir)) \
        if os.path.isfile(image_dir+os.path.sep+image_fn) and filter_kw in image_fn]
    for image_fn in image_fn_list:
        image = freeimage.read(image_dir+os.path.sep+image_fn)
        if master_vig is not None:
            image[master_vig==0]=0
        image_list.append(image)

    return OrderedDict((image_fn[:find_char(image_fn,'.')[-1]],image) for image_fn,image in zip(image_fn_list,image_list))

def load_images_fromfns(image_fns, filter_kw = '', master_vig = None):
    image_list = []
    for image_fn in image_fns:
        image = freeimage.read(image_fn)
        if master_vig is not None:
            image[master_vig==0]=0
        image_list.append(image)
    return image_list

if __name__ == "__main__":
    expt_superdir_hmasks = '/media/Data/Work/ZPLab/WormImages/spe9Acquisition_DevImages_hmasks/'
    hmask_info = json.loads(open(expt_superdir_hmasks+os.path.sep+'data_info.json').read())
    
    all_worm_masks = {'dt':[], 'h_masks':[],'auto_masks':[]}
    
    for expt_dir in sorted(os.listdir(expt_superdir_hmasks)):
        print(expt_dir)
        if os.path.isdir(expt_superdir_hmasks+os.path.sep+expt_dir):
            for worm_dir in sorted(os.listdir(expt_superdir_hmasks+os.path.sep+expt_dir)):
                if os.path.isdir(expt_superdir_hmasks+os.path.sep+expt_dir+os.path.sep+worm_dir):
                    h_masks = load_images_fromdir(expt_superdir_hmasks+os.path.sep+expt_dir+os.path.sep+worm_dir,filter_kw='hmask')
                    auto_masks = load_images_fromfns([
                        hmask_info['expt_fps'][int(expt_dir)]+os.path.sep+worm_dir+os.path.sep+img_fn[0:15]+' mask.png' \
                        for img_fn in h_masks.keys()])
                    all_worm_masks['dt'].extend([img_fn[0:15] for img_fn in h_masks.keys()])
                    all_worm_masks['h_masks'].extend(h_masks.values())
                    all_worm_masks['auto_masks'].extend(auto_masks)
    hmask_sizes = np.array([np.count_nonzero(mask>0) for mask in all_worm_masks['h_masks']])
    auto_sizes = np.array([np.count_nonzero(mask>0) for mask in all_worm_masks['auto_masks']])
    
    # Do summary over all images
    plt.ion()
    fig_h = plt.figure()
    ax_h = fig_h.gca()
    ax_h.scatter(hmask_sizes, auto_sizes)
    ax_h.set_xlabel('Human mask size (px)')
    ax_h.set_ylabel('auto mask size (px)')
    ax_h.set_title('Spearman r^2: {:.2f}\nPearson r^2: {:.2f}'.format(
        scipy.stats.spearmanr(hmask_sizes, auto_sizes)[0]**2,
        scipy.stats.pearsonr(hmask_sizes, auto_sizes)[0]**2))
    
    # Do summary per larval stage
    fig_h, ax_h = plt.subplots(5,1)
    for stg in range(5):    # first is hatch
        ax_h[stg].scatter(hmask_sizes[stg+1::6],auto_sizes[stg+1::6])
        ax_h[stg].set_xlabel('Human mask size (px)')
        ax_h[stg].set_ylabel('auto mask size (px)')
        ax_h[stg].set_title('{} Spearman r^2: {:.2f}\nPearson r^2: {:.2f}'.format(
            'L'+str(stg) if stg>0 else 'Hatch', \
            scipy.stats.spearmanr(hmask_sizes[stg+1::6], auto_sizes[stg+1::6])[0]**2, \
            scipy.stats.pearsonr(hmask_sizes, auto_sizes)[0]**2))
