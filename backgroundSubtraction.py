# -*- coding: utf-8 -*-
"""
Created on Sun May 31 19:06:30 2015

@author: Willie
"""

import os
import freeimage
import numpy as np
import shutil

from zplib.image import mask as zplib_image_mask
import scipy.ndimage.morphology


def overallBackgroundSubtract(data_dir, match_string, new_string, temporal_radius, save_dir = '', save_dir2 = '', save_dir3 = '', demonstration_mode = False):
    '''
    Do background subtraction to find worms. This uses only past data, masking out the worms to create a background that won't disappear once the worm stops  moving.
    '''
    if save_dir == '':
        save_dir = data_dir 
    
    my_files = sorted(os.listdir(data_dir))
    my_files = [a_file for a_file in my_files if match_string in a_file]

    # Intialize my special background context.
    if save_dir != data_dir:
        temp_folder = save_dir + os.path.sep + 'temp'
        try:
            os.stat(temp_folder)
        except: 
            os.mkdir(temp_folder)
        for i in range(0, temporal_radius):
            shutil.copy(data_dir + os.path.sep + my_files[i], temp_folder + os.path.sep + my_files[i])

    # Run the actual simple subtraction, saving out masked files.
    context_files = [freeimage.read(data_dir + os.path.sep + my_files[j]) for j in range(0, temporal_radius)]
    for i in range(temporal_radius, len(my_files)):
        real_raw_file = freeimage.read(data_dir + os.path.sep + my_files[i])
        raw_file = real_raw_file.copy()     
        context_files.append(raw_file)
        (foreground_file, background_file) = simple_running_median_subtraction(raw_file, context_files)
        
        thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
        final_mask = clean_dust_and_holes(thresholded_mask)

        raw_file[final_mask.astype('bool')] = background_file[final_mask.astype('bool')]        

        if demonstration_mode:
            freeimage.write(real_raw_file, save_dir + os.path.sep + my_files[i])
            freeimage.write(background_file, save_dir2 + os.path.sep + my_files[i])
            freeimage.write(final_mask, save_dir3 + os.path.sep + my_files[i])

        if not demonstration_mode:
            if save_dir != data_dir:
                freeimage.write(real_raw_file, temp_folder + os.path.sep + my_files[i])
            freeimage.write(final_mask, save_dir + os.path.sep + my_files[i].replace(match_string, new_string))
            #freeimage.write(background_file, save_dir + os.path.sep + my_files[i].replace(match_string, 'bgfile'))
            #freeimage.write(np.abs(raw_file-background_file).astype('uint16'), save_dir + os.path.sep + my_files[i].replace(match_string, 'bgdiff'))
            #freeimage.write(np.min([[raw_file.astype('int16') - background_file.astype('int16')],[-1*raw_file.astype('int16') + background_file.astype('int16')]],axis=0)[0].astype('uint16'), save_dir + os.path.sep + my_files[i].replace(match_string, 'bgdiff'))
        context_files = context_files[1:]

    # Run another small chunk of background subtraction backwards to fill out the early range.
    context_files = [freeimage.read(data_dir + os.path.sep + my_files[j]) for j in reversed(range(temporal_radius, temporal_radius*2))]
    for i in reversed(range(0, temporal_radius)):
        real_raw_file = freeimage.read(data_dir + os.path.sep + my_files[i])
        raw_file = real_raw_file.copy()     
        context_files.append(raw_file)
        (foreground_file, background_file) = simple_running_median_subtraction(raw_file, context_files)
        
        thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
        final_mask = clean_dust_and_holes(thresholded_mask)

        raw_file[final_mask.astype('bool')] = background_file[final_mask.astype('bool')]        

        if demonstration_mode:
            freeimage.write(real_raw_file, save_dir + os.path.sep + my_files[i])
            freeimage.write(background_file, save_dir2 + os.path.sep + my_files[i])
            freeimage.write(final_mask, save_dir3 + os.path.sep + my_files[i])

        if not demonstration_mode:
            if save_dir != data_dir:
                freeimage.write(real_raw_file, temp_folder + os.path.sep + my_files[i])
            freeimage.write(final_mask, save_dir + os.path.sep + my_files[i].replace(match_string, new_string))    
            
        context_files = context_files[1:]
    return
    
def overallBackgroundSubtract_lawnremoving(data_dir, match_string, new_string, temporal_radius, lawn_image, save_dir = '', save_dir2 = '', save_dir3 = '', debug_mode = False):
    '''
    Do background subtraction to find worms. This uses only past data, masking out the worms to create a background that won't disappear once the worm stops  moving.
    '''
    if save_dir == '':
        save_dir = data_dir 
    
    my_files = sorted(os.listdir(data_dir))
    my_files = [a_file for a_file in my_files if match_string in a_file]

    # Intialize my special background context.
    if save_dir != data_dir:
        try:
            os.stat(save_dir)
        except:
            os.mkdir(save_dir)
        
        temp_folder = save_dir + os.path.sep + 'temp'
        try:
            os.stat(temp_folder)
        except: 
            os.mkdir(temp_folder)
        for i in range(0, temporal_radius):
            shutil.copy(data_dir + os.path.sep + my_files[i], temp_folder + os.path.sep + my_files[i])

    # Run the actual simple subtraction, saving out masked files.
    context_files = [freeimage.read(data_dir + os.path.sep + my_files[j]) for j in range(0, temporal_radius)]
    for i in range(temporal_radius, len(my_files)):
        real_raw_file = freeimage.read(data_dir + os.path.sep + my_files[i])
        raw_file = real_raw_file.copy()     
        context_files.append(raw_file)
        (foreground_file, background_file) = simple_running_median_subtraction_lawnremoving(raw_file, context_files,lawn_image)
        
        thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
        final_mask = clean_dust_and_holes(thresholded_mask)

        raw_file[final_mask.astype('bool')] = background_file[final_mask.astype('bool')]        

        if debug_mode:
            last_dot = [idx for idx, ch in enumerate(my_files[i]) if ch == '.'][-1]
            
            freeimage.write(real_raw_file, save_dir + os.path.sep + my_files[i][:last_dot]+'_raw'+my_files[i][last_dot:])
            freeimage.write(background_file, save_dir + os.path.sep + my_files[i][:last_dot]+'_back'+my_files[i][last_dot:])
            freeimage.write(foreground_file, save_dir + os.path.sep + my_files[i][:last_dot]+'_fore'+my_files[i][last_dot:])
            freeimage.write(final_mask, save_dir + os.path.sep + my_files[i][:last_dot]+'_final'+my_files[i][last_dot:])
        else:
            if save_dir != data_dir:
                freeimage.write(real_raw_file, temp_folder + os.path.sep + my_files[i])
            freeimage.write(final_mask, save_dir + os.path.sep + my_files[i].replace(match_string, new_string).replace('tiff','png'))
            #freeimage.write(background_file, save_dir + os.path.sep + my_files[i].replace(match_string, 'bgfile'))
            #freeimage.write(np.abs(raw_file-background_file).astype('uint16'), save_dir + os.path.sep + my_files[i].replace(match_string, 'bgdiff'))
            #freeimage.write(np.min([[raw_file.astype('int16') - background_file.astype('int16')],[-1*raw_file.astype('int16') + background_file.astype('int16')]],axis=0)[0].astype('uint16'), save_dir + os.path.sep + my_files[i].replace(match_string, 'bgdiff'))
        context_files = context_files[1:]

    # Run another small chunk of background subtraction backwards to fill out the early range.
    context_files = [freeimage.read(data_dir + os.path.sep + my_files[j]) for j in reversed(range(temporal_radius, temporal_radius*2))]
    for i in reversed(range(0, temporal_radius)):
        real_raw_file = freeimage.read(data_dir + os.path.sep + my_files[i])
        raw_file = real_raw_file.copy()     
        context_files.append(raw_file)
        (foreground_file, background_file) = simple_running_median_subtraction(raw_file, context_files)
        
        thresholded_mask = percentile_floor(foreground_file, threshold_proportion = 0.975)
        final_mask = clean_dust_and_holes(thresholded_mask)

        raw_file[final_mask.astype('bool')] = background_file[final_mask.astype('bool')]        

        if debug_mode:
            last_dot = [idx for idx, ch in enumerate(my_files[i]) if ch == '.'][-1]
            
            freeimage.write(real_raw_file, save_dir + os.path.sep + my_files[i][:last_dot]+'_raw'+my_files[i][last_dot:])
            freeimage.write(background_file, save_dir + os.path.sep + my_files[i][:last_dot]+'_back'+my_files[i][last_dot:])
            freeimage.write(foreground_file, save_dir + os.path.sep + my_files[i][:last_dot]+'_fore'+my_files[i][last_dot:])
            freeimage.write(final_mask, save_dir + os.path.sep + my_files[i][:last_dot]+'_final'+my_files[i][last_dot:])
        else:
            if save_dir != data_dir:
                freeimage.write(real_raw_file, temp_folder + os.path.sep + my_files[i])
            freeimage.write(final_mask, save_dir + os.path.sep + my_files[i].replace(match_string, new_string).replace('tiff','png'))    
            
        context_files = context_files[1:]
    return

def clean_dust_and_holes(dusty_pic):
    '''
    Picks out the largest object in dusty_mask and fills in its holes, returning cleaned_mask.
    '''
    my_dtype = dusty_pic.dtype
    dust_mask = np.invert(zplib_image_mask.get_largest_object(dusty_pic))       
    dusty_pic[dust_mask] = 0
    cleaned_mask = simple_floor(zplib_image_mask.fill_small_area_holes(dusty_pic, 90000).astype(my_dtype), 1)
    return cleaned_mask

def simple_floor(focal_image, threshold_value):
    ''' 
    Takes a grayscale focal image (in the form of a numpy array), and sets to zero (black) all values below the threshold value.
    '''
    max_value = -1
    binary_image = focal_image.copy()
    binary_image[binary_image < threshold_value] = 0 
    binary_image[binary_image >= threshold_value] = max_value 
    return binary_image
    
def percentile_floor(focal_image, threshold_proportion):
    '''
    Takes a grayscale focal image (in the form of a numpy array), and sets to zero (black) all values below the percentile indicated by threshold_proportion.
    '''
    max_value = -1
    binary_image = focal_image.copy()
    threshold_value = int(np.percentile(binary_image, threshold_proportion*100))
    binary_image[binary_image < threshold_value] = 0 
    binary_image[binary_image >= threshold_value] = max_value   
    return binary_image

def simple_running_median_subtraction(focal_image, background_images):
    ''' 
    Takes a focal image and a list of background images (grayscale, in the form of numpy arrays), and returns the focal image with the background subtracted. This simply takes the median value of each pixel to construct a background.   
    '''
    median_image = median_image_from_list(background_images)
    foreground_only = abs(focal_image.astype('int16') - median_image.astype('int16')).astype('uint16')
    #foreground_only = np.min([[focal_image.astype('int16') - median_image.astype('int16')],[-1*focal_image.astype('int16') + median_image.astype('int16')]],axis=0)[0].astype('uint16')
    return (foreground_only, median_image)

def simple_running_median_subtraction_lawnremoving(focal_image, background_images, lawn_image):
    ''' 
    Takes a focal image and a list of background images (grayscale, in the form of numpy arrays), and returns the focal image with the background subtracted. This simply takes the median value of each pixel to construct a background.   
    '''
    median_image = median_image_from_list(background_images)
    foreground_only = abs(focal_image.astype('int16') - median_image.astype('int16')).astype('uint16')
    foreground_only[np.invert(scipy.ndimage.morphology.binary_erosion(lawn_image>0,iterations=60))]=0
    #foreground_only = np.min([[focal_image.astype('int16') - median_image.astype('int16')],[-1*focal_image.astype('int16') + median_image.astype('int16')]],axis=0)[0].astype('uint16')
    return (foreground_only, median_image)

def median_image_from_list(background_images):
    ''' 
    Takes a list of background images (grayscale, in the form of numpy arrays), and returns an image constructed by taking the median value of each pixel.
    '''
    big_array = np.array([background_image for background_image in background_images])
    median_image = np.median(big_array, axis = 0).astype('uint16')
    return median_image


def main():
    return

if __name__ == "__main__":
    main()
