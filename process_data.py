import datetime
import freeimage
import numpy as np
import json
import os
import shutil

import zplib
import zplib.image.mask as zplib_image_mask

import skimage.feature
import skimage.morphology
import skimage.filters.rank
import scipy.ndimage.morphology
import cv2

#import MeasureFluorescence.wormFinding.edgeMorphology
#import MeasureFluorescence.wormFinding.backgroundSubtraction as backgroundSubtraction
import backgroundSubtraction

import time

import concurrent.futures
import multiprocessing

import annotation_file

def import_json_fromfn(fn):
    return json.loads(open(fn).read())

def extract_datetime_fromstr(time_str):
    return datetime.datetime.strptime(time_str,'%Y-%m-%dt%H%M')
    
def list_image_files(dir_path):
    return [dir_file for dir_file in sorted(os.listdir(dir_path)) if os.path.isfile(dir_path+os.path.sep+dir_file) and ('.png' in dir_file or '.tiff' in dir_file)]
'''
    For making measurements:
        X Calibrate images
        X Generate masks (using old bgsubtract code)
        X Measure areas of masks and write out information per timepoint
'''


class Experiment:
    def __init__(self, expt_path):
        self.expt_path = expt_path
        self.expt_mdata = import_json_fromfn(expt_path+os.path.sep+'experiment_metadata.json')
        
        #self.acquired_wells = [well for well in sorted(os.listdir(expt_path)) if os.path.isdir(expt_path+os.path.sep+well+os.path.sep) and well != 'calibrations']
        self.acquired_wells = [well for well in sorted(os.listdir(expt_path)) if os.path.isdir(expt_path+os.path.sep+well+os.path.sep) and not well[0].isalpha()]
        position_mdata = {}
        for well in self.acquired_wells:
            position_mdata[well]=import_json_fromfn(expt_path+os.path.sep+well+os.path.sep+'position_metadata.json')
    
    #def get_expt_path(self):
        #return self.expt_path
    
    def get_expt_times(self):
        return [extract_datetime_fromstr(time_point) for time_point in self.expt_mdata['timepoints']]
    
    #def get_acquired_wells(self):
        #return self.acquired_wells
    
    def get_expt_time_offsets(self):
        expt_times = self.get_expt_times()
        return [(expt_tp - expt_times[0]).total_seconds() for expt_tp in expt_times]
        
    # TODO: Continue working on (have to figure out how to handle multiple images (e.g. bf and fl)
    def get_well_images(self, well, start_idx = 0, stop_idx = None):
        #well_image_files = list_image_files(self.expt_path+os.path.sep+well)
        #if stop_idx is None: stop_idx = len(well_image
        pass
    
    #TODO FINISH AND FIX (looks like issue with use of wells below)
    def measure_growth_forexpt(self):
        '''
        worm_size = ...
        for well in self.acquired_wells:
            well_masks = self.build_masks_forwell_bgsubtract()
            for idx, time_point in enumerate(self.expt_mdata['timepoints']):
                worm_size[well, idx] = np.count_nonzero(well_masks[idx])
        '''
        worm_sizes = np.array([np.count_nonzero(worm_mask) for well in self.acquired_wells for worm_mask in self.build_masks_forwell_bgsubtract()])
        
        expt_times = self.get_expt_times()
        time_diffs = [(right_tp-left_tp).total_seconds()/3600 for left_tp,right_tp in zip(expt_times[0:len(expt_times)-1],expt_times[1:len(expt_times)])]
        worm_growth_rates = np.diff(worm_sizes, axis=1)/np.tile(time_diffs, len(self.acquired_wells))
        return (worm_growth_rates, worm_sizes)
    
    def measure_growthrate_forexpt_batch(self, match_string='mask'):
        worm_sizes = self.measure_size_forexpt_batch()
        expt_times = self.get_expt_times()
        time_diffs = [(right_tp-left_tp).total_seconds()/3600 for left_tp,right_tp in zip(expt_times[0:len(expt_times)-1],expt_times[1:len(expt_times)])]
        #worm_growth_rates = np.diff(worm_sizes, axis=1)/np.tile(time_diffs, len(self.acquired_wells))
        worm_growth_rates = [np.diff(individual_worm_size)/time_diffs[0:len(individual_worm_size)-1] for individual_worm_size in worm_sizes]  # Limit time_diffs in case worm wasn't acquired over entire duration of expt
        return worm_growth_rates  # row/axis 0 - well; col/axis 1 - time
        
    def measure_size_forexpt_batch(self, match_string='mask'):
        return [np.array([np.count_nonzero(worm_mask) for worm_mask in [freeimage.read(self.expt_path+os.path.sep+well+os.path.sep+mask_file) for mask_file in sorted(os.listdir(self.expt_path+os.path.sep+well)) if match_string in mask_file]]) \
            for well in self.acquired_wells]
    
    def calibrate_bf_images_batch(self, data_match_str='', cali_match_str='', save_path = '', save_str = '_corrected', save_images = False, remake_master = False):
        if save_images and save_path == '': save_path = self.expt_path
        if not os.path.exists(save_path): 
            print('save path not found; making save path at:'+save_path)
            os.mkdir(save_path)
        print('saving images in: '+save_path)

        if remake_master:
            master_v_mask = self.make_master_vignette_mask(save_mask=True)
        else:
            master_v_mask = self.read_master_vignette_mask()
        print('loaded master vignette mask')
        
        median_ref_intensity = np.median([self.expt_mdata['brightfield metering'][dataset]['ref_intensity'] for dataset in self.expt_mdata['brightfield metering'].keys()])
        
        cali_path = self.expt_path+os.path.sep+'calibrations'+os.path.sep
        for well in self.acquired_wells:   # Use directory names (i.e. worm numbers) as the labels to iterate over
            print('Calibrating well '+well)
            well_path = self.expt_path+os.path.sep+well
            for well_file in list_image_files(well_path):
                if data_match_str in well_file:
                    metering_data = self.expt_mdata['brightfield metering'][well_file[0:15]]
                    trans_img = freeimage.read(well_path+os.path.sep+well_file)*median_ref_intensity/metering_data['ref_intensity']* \
                        freeimage.read(cali_path+os.path.sep+metering_data['cal_image_prefix']+' bf_flatfield.tiff')*master_v_mask
                    
                    if save_images:
                        if not os.path.exists(save_path+os.path.sep+well): os.mkdir(save_path+os.path.sep+well)
                        save_fn = save_path+os.path.sep+well+os.path.sep+well_file[:well_file.index('.')]+save_str+well_file[well_file.index('.'):]
                        freeimage.write(trans_img.astype('uint16'), save_fn)
    
    def make_master_vignette_mask(self, save_mask=False):
        cali_path = self.expt_path+os.path.sep+'calibrations'
        master_v_mask = np.logical_and.reduce(np.array([freeimage.read(cali_path+os.path.sep+cali_file) > 0 for cali_file in list_image_files(cali_path) if ('bf_flatfield') in cali_file]))
        if save_mask: 
            freeimage.write(master_v_mask.astype('uint8'),cali_path+os.path.sep+'master_v_mask.tiff')
        return master_v_mask
    
    def read_master_vignette_mask(self):
        return freeimage.read(self.expt_path+os.path.sep+'calibrations'+os.path.sep+'master_v_mask.tiff')
    
    def read_well_lawn_mask(self, well):
        return freeimage.read(self.expt_path+os.path.sep+well+os.path.sep+'great_lawn.png')
    
    def build_mask_bgsubtract_batch(self):
        # Assumes 'great_lawn.png' saved in each worm directory
        for well in self.acquired_wells:
            #print('building masks for well {}'.format(well))
            backgroundSubtraction.overallBackgroundSubtract_lawnremoving(
                self.expt_path+os.path.sep + well+os.path.sep,'bf_corrected', 'mask', 10,self.read_well_lawn_mask(well))
    
    def build_mask_bgsubtract_perwell_batch(self, well, save_dir, debug_mode):
        backgroundSubtraction.overallBackgroundSubtract_lawnremoving(
            self.expt_path+os.path.sep + well+os.path.sep,'bf_corrected', 'mask', 10,self.read_well_lawn_mask(well), save_dir = save_dir+os.path.sep+well, debug_mode = debug_mode)

    
    # TODO GET TO A GOOD STATE...
    def make_lawn_forwell(self, well, master_v_mask):
        '''
        #Make a mega lawn mask for use with all images of a worm. This avoids problems with detecting the relatively faint edge of the lawn, which depends on the exact focus plane and the brightness of the lamp.
        '''
        bf_files = [self.expt_path+os.path.sep+well+os.path.sep+bf_file for bf_file in list_image_files(self.expt_path+os.path.sep+well) if 'bf_corrected' in bf_file][0:10]
        file_lawns = np.array([make_lawn_forimage(bf_file, master_v_mask) for bf_file in bf_files])
        well_lawn = np.logical_or.reduce(file_lawns, axis=0) 
        freeimage.write(well_lawn.astype('uint16'), self.expt_path + os.path.sep + 'calibrations' + os.path.sep + 'lawn_'+well+'.png')
        return well_lawn
    
    def label_time_perstage(self, a_file, metadata_list):
        '''
            Return, for all worms in the experiment, labels for each timepoint denoting what development stage the individual was in
            
            Input:
                a_file - annotation file for experiment
                metadata_list: (character delimiter: experiment metadata file paths)
            
            Output:
                time_labels: Matrix (shape = num_worms x num_timepoints) giving labels for dev. stage at the given timepoint
                    Labels: 0 - Egg; 1-4 - larval stgs; 5 - adult
        '''
        time_labels = np.zeros((len(self.acquired_wells), len(self.get_expt_times())))
        events=a_file.data_as_timestamps(metadata_list)
        expt_times = self.get_expt_time_offsets()
        event_times_forbinning = np.array([events[a_tag] for a_tag in ['Hatch', 'L1 ecdysis', 'L2 ecdysis','L3 ecdysis', 'L4 ecdysis']])
        for worm in self.acquired_wells:
            worm_idx = np.where(events['Worm']==('/'+worm))[0]
            time_labels[worm_idx,:] = np.digitize(expt_times, event_times_forbinning[:,worm_idx][:,0])
        
        return time_labels
        
    '''
    #TODO: Update based on batch
    def calibrate_single_image(self, time_point, well, master_v_mask, microscopy_type='bf', save_path = '', save_str = '', image_type = 'tiff'):
        if save_path != '' and not os.path.exists(save_path):
            print('save path not found; making save path at:'+save_path)
            os.mkdir(save_path)
        
        cali_path = self.expt_dir+os.path.sep+'calibrations'+os.path.sep
        median_ref_intensity = np.median([self.expt_mdata['brightfield_metering'][dataset]['ref_intensity'] for dataset in expt_mdata['brightfield metering'].keys()])
        
        # Build a master vignette mask from all of the calibration files
        #master_v_mask = np.logical_and.reduce(np.array([freeimage.read(cali_path+os.path.sep+cali_file) > 0 for cali_file in os.listdir(cali_path) if (microscopy_type+'_flatfield') in cali_file]))
        metering_data = self.expt_mdata['brightfield metering'][well]
        image_file = self.expt_dir+os.path.sep+well+os.path.sep+time_point+' '+microscopy_type+'.'+image_type
        cali_file = cali_path+os.path.sep+(self.expt_mdata['brightfield metering'][well]['cal_image_prefix']+' '+microscopy_type+'_flatfield.tiff')
        out_img = freeimage.read(image_file)*median_ref_intensity/metering_data['ref_intensity']* \
            freeimage.read(cali_file)*master_v_mask
        if save_path != '':
            save_fn = save_path+os.path.sep+d_file[:-4]+save_str+d_file[-4:]
            freeimage.write(out_img.astype('uint16'), save_fn)
        return out_img
    
    #TODO: Update based on batch
    def calibrate_images_forwell(self, well, master_v_mask, microscopy_type, save_path = '', save_str = '', image_type = 'tiff', return_images = False):
        if not return_images:
            [self.calibrate_single_image(time_point, well, microscopy_type, save_path, save_str, image_type) for time_point in self.expt_mdata['timepoints']]
        else:
            return np.array([self.calibrate_single_image(time_point, well, save_path, save_str, image_type) for time_point in self.expt_mdata['timepoints']])
    
    #TODO: Update based on batch
    def calibrate_images_forexpt(self, microscopy_type, save_path = '', save_str = '', image_type = 'tiff', return_images = False):
        cali_path = self.expt_dir+os.path.sep+'calibrations'+os.path.sep
        master_v_mask = np.logical_and.reduce(np.array([freeimage.read(cali_path+os.path.sep+cali_file) > 0 for cali_file in os.listdir(cali_path) if (microscopy_type+'_flatfield') in cali_file]))
        if not return_images:
            [self.calibrate_images_forwell(well, master_v_mask, microscopy_type, save_path, save_str, image_type) for well in self.acquired_wells]
        else:
            return np.array([self.calibrate_images_forwell(well, master_v_mask, microscopy_type, save_path, save_str, image_type, return_images=True) for well in self.acquired_wells])
    
    
    def build_mask_bgsubtract_singleframe(self, data_image, bg_files, save_path = '', image_type = 'tiff', return_images = False):
        #TODO
        background_image = np.median(bg_files,axis=0)
        foreground_image = np.abs(data_image-background_image)
        pass
    
    def build_mask_bgsubtract_timepoint(self, well, time_point, temporal_radius, save_path = '', image_type = 'tiff', return_images=False):
        
        image_file = self.expt_path+os.path.sep+well+os.path.sep+time_point+' bf_corrected.'+image_type
        
        # Find timepoint in list of timepoints as well as timepoints within temporal radius
        bg_images = np.array([])
        time_idx = self.expt_mdata['timepoints'].index(time_point)
        if time_idx < temporal_radius-1: bg_time_idx = [time_idx+1:
            
        else: #Use trailing images
            bg_images = [freeimage.read(
        
        pass
        
    
    def build_masks_forwell_bgsubtract():
        #TODO
        pass
    '''

#for worm_idx in np.arange(stop = np.len(events['Hatch'])):
     ##Since annotations use frame number (0:) to give stage, provide the range to use for binning)
    #if str(worm_idx).zfill(len(self.acquired_wells[0])) in self.acquired_wells:
        #time_labels[worm_idx,:] = np.digitize(np.arange(stop=len(self.get_expt_times())), event_times_forbinning[worm_idx])


#For making a lawn from a single image
#def make_lawn_forimage(image_fn, master_v_mask, debug_mode = False):
    #'''
    #Find the bacterial lawn in one worm image.
    #'''
    ## Prepare a worm image for use in lawn-finding.
    #worm_image = freeimage.read(image_fn)
    ##proc_image = skimage.filters.rank.median(zplib.image.colorize.scale(worm_image, output_max=2**16-1).astype('uint16'), skimage.morphology.disk(1))
    #proc_image = cv2.medianBlur(zplib.image.colorize.scale(worm_image, output_max=2**16-1).astype('uint16'),3)  # Do blur over 8-neighborhood
    
    ##Remove extraneous edges and out-of-lawn junk by finding the lawn and also applying an "ultra-vignette" mask.
    ##my_edges = skimage.feature.canny(proc_image, sigma = 0.05)
    
    ##ultra_vignette = scipy.ndimage.morphology.binary_erosion(master_v_mask, iterations = 10)
    #ultra_vignette = scipy.ndimage.morphology.binary_erosion(master_v_mask, iterations = 200)
    #edge_image[np.invert(ultra_vignette)] = False 
    #edge_image = scipy.ndimage.morphology.binary_dilation(edge_image, iterations = 10)
    
    #lawn_image = scipy.ndimage.morphology.binary_fill_holes(edge_image)
    #lawn_image = zplib.image.mask.get_largest_object(lawn_image).astype('bool')
    #lawn_image = scipy.ndimage.morphology.binary_erosion(lawn_image, iterations = 10)
    #lawn_image = zplib.image.mask.get_largest_object(lawn_image).astype('bool')
    
    #if debug_mode: return (lawn_image, edge_image
    #return my_lawn


def make_lawn_forimage(image_fn, master_v_mask, debug_mode = False):
    worm_image = freeimage.read(image_fn)
    thr_image = worm_image < np.percentile(worm_image,5)
    
    ultra_vignette = scipy.ndimage.morphology.binary_erosion(master_v_mask, iterations = 10)
    thr_image[np.invert(ultra_vignette)] = False
    
    lawn_image = scipy.ndimage.morphology.binary_dilation(thr_image, iterations = 3)
    lawn_image = scipy.ndimage.morphology.binary_fill_holes(lawn_image)
    lawn_image = zplib.image.mask.get_largest_object(lawn_image).astype('bool')
    lawn_image = scipy.ndimage.morphology.binary_erosion(lawn_image, iterations = 3)
    lawn_image = zplib.image.mask.get_largest_object(lawn_image).astype('bool')
    
    if debug_mode: return (lawn_image, thr_image, worm_image < np.percentile(worm_image,5))
    return lawn_image

def parallel_expt_maskbuilder(expt, save_dir, debug_mode):
    num_workers = max(multiprocessing.cpu_count()-4,3)
    #num_workers=3
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        job_data = [executor.submit(expt.build_mask_bgsubtract_perwell_batch, well, save_dir,debug_mode) for well in expt.get_acquired_wells()]
    concurrent.futures.wait(job_data)
    if any([not job.done() for job in job_data]):
        print('ended in a bad state; check job_data')
    return job_data
    
    
'''
def calibrate_images(expt_path, data_match_str = '', cali_match_str = 'bf_flatfield', save_path = '', save_str = '_corrected'):
    if save_path == '': save_path = expt_path
    if not os.path.exists(save_path): 
        print('save path not found; making save path at:' + save_path)
        os.mkdir(save_path)
    
    expt_mdata = import_json_fromfn(expt_path+os.path.sep+'experiment_metadata.json')
    cali_path = expt_path+os.path.sep+'calibrations'
    median_ref_intensity = np.median([expt_mdata['brightfield_metering'][dataset]['ref_intensity'] for dataset in expt_mdata['brightfield metering'].keys()])
    
    # Extract the data files to be processed in the data directory and their corresponding calibration files
    organized_files = {}
    for pos_label in [expt_dir for expt_dir in os.listdir(expt_path) if os.path.isdir(expt_dir)]:   # Use directory names (i.e. worm numbers) as the labels to iterate over
        data_path = expt_path+os.path.sep+pos_label
        for d_file in sorted(os.listdir(data_path+os.path.sep+pos_label)):
            if data_match_str in d_file and os.path.isfile(data_path+os.path.sep+pos_label):
                organized_files[d_file[0:15]]={'d_file':d_file}
        for c_file in sorted(os.listdir(cali_path+os.path.sep+'calibrations')):
            if cali_match_str in c_file and os.path.isfile(cali_path+os.path.sep+c_file):
                organized_files[c_file[0:15]]['c_file'] = c_file
        master_v_mask = np.logical_and.reduce(np.array([freeimage.read(cali_path+os.path.sep+organized_files[dataset]['c_file'])>0 for dataset in organized_files.keys()]))
        
        # Note: Calibration image is written (float32) so that one MULTIPLIES it against the data image (vs. dividing it)
        for dataset in organized_files.keys():
            metering_data = expt_mdata['brightfield metering'][dataset]
            [d_file, c_file] = [organized_files[dataset]['d_file'], organized_files[metering_data['cal_image_prefix']]['c_file']]
            trans_img = freeimage.read(data_path+os.path.sep+d_file)*median_ref_intensity/metering_data['ref_intensity']* \
                freeimage.read(cali_path+os.path.sep+c_file)*master_v_mask
            save_fn = save_path+os.path.sep+d_file[:-4]+save_str+d_file[-4:]
            freeimage.write(trans_img.astype('uint16'), save_fn)

# TODO: Write functionality to calibrate single image
'''


#~ def make_mega_lawn(worm_subdirectory, super_vignette):
    #~ '''
    #~ Make a mega lawn mask for use with all images of a worm. This avoids problems with detecting the relatively faint edge of the lawn, which depends on the exact focus plane and the brightness of the lamp.
    #~ '''
    #~ my_bf_files = [worm_subdirectory + o s.path.sep + a_bf for a_bf in os.listdir(worm_subdirectory) if a_bf[-6:] == 'bf.png']

    #~ my_workers = min(multiprocessing.cpu_count() - 1, 60)

    #~ with concurrent.futures.ProcessPoolExecutor(max_workers = my_workers) as executor:
        #~ chunk_masks = [executor.submit(lawn_maker, a_bf_file, super_vignette.copy()) for a_bf_file in my_bf_files]

    #~ concurrent.futures.wait(chunk_masks)
    #~ chunk_masks = [a_job.result() for a_job in chunk_masks]

    #~ mega_lawn = np.max(np.array(chunk_masks), axis = 0)
    #~ mega_lawn = mega_lawn.astype('uint8')   
    #~ mega_lawn[mega_lawn > 0] = -1   
    #~ freeimage.write(mega_lawn, worm_subdirectory + os.path.sep + 'great_lawn.png')
    #~ return mega_lawn

#~ def lawn_maker(a_bf_file, super_vignette):
    #~ '''
    #~ #Find the bacterial lawn in one worm image.
    #~ '''
    #~ # Prepare a worm image for use in lawn-finding.
    #~ worm_image = freeimage.read(a_bf_file)
    #~ renormalized_image = imageOperations.renormalize_image(worm_image)
    #~ renormalized_image = skimage.filters.rank.median(renormalized_image, skimage.morphology.disk(1))
    
    #~ # Remove extraneous edges and out-of-lawn junk by finding the lawn and also applying an "ultra-vignette" mask.
    #~ ultra_vignette = scipy.ndimage.morphology.binary_erosion(super_vignette, iterations = 10)
    #~ my_edges = skimage.feature.canny(renormalized_image, sigma = 0.05)
    #~ my_edges[np.invert(ultra_vignette)] = False 
    #~ my_edges = scipy.ndimage.morphology.binary_dilation(my_edges, iterations = 10)
    #~ my_lawn = scipy.ndimage.morphology.binary_fill_holes(my_edges)
    #~ my_lawn = zplib_image_mask.get_largest_object(my_lawn).astype('bool')
    #~ my_lawn = scipy.ndimage.morphology.binary_erosion(my_lawn, iterations = 10)
    #~ my_lawn = zplib_image_mask.get_largest_object(my_lawn).astype('bool')
    #~ return my_lawn

def grab_images_for_dev_stages(expt_fps, ann_fps, expt_mds, save_path, match_string='corrected'):
    bad_worm_kws = ['FERTILITY', 'Nh', 'NO HATCH', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS']
    
    
    my_ann_files = [annotation_file.AnnotationFile(ann_fp,annotation_prefix='D') for ann_fp in ann_fps]
    timestamped_data = [[] for idx in range(len(expt_fps))]
    for expt_num, [expt_md_fps, ann_file] in enumerate(zip(expt_mds, my_ann_files)):
        timestamped_data[expt_num] = ann_file.data_as_timestamps(expt_md_fps)
    
    
    # Save json readme for each expt to save_path
    info_dict = {'ann_fps':ann_fps, 'expts_fps':expt_fps, 'expt_mds':expt_mds}
    with open(save_path+os.path.sep+'data_info.json','w') as info_file:
        json.dump(info_dict,info_file)
    
    for expt_num, [expt, expt_events,ann_file] in enumerate(zip([Experiment(expt_fp) for expt_fp in expt_fps], timestamped_data,my_ann_files)):
        if not os.path.exists(save_path+os.path.sep+str(expt_num).zfill(2)):
            os.mkdir(save_path+os.path.sep+str(expt_num).zfill(2))
        
        first_worm_num = int(expt_events['Worm'][0][-2:])   # Add in adjustment for the first index of the worms not being at 0
        viable_worm = (expt_events['Hatch']!=-1) \
            & (expt_events['Death']!=-1) \
            & np.array([not any([kw in note for kw in bad_worm_kws]) for note in expt_events['Notes']])
        hatched_on_corral = timestamped_data[expt_num]['Hatch'] != 0
        dev_labels = expt.label_time_perstage(ann_file,expt_md_fps)
        print(np.where(viable_worm & hatched_on_corral)[0])
        
        for worm_num, worm in enumerate(np.where(viable_worm & hatched_on_corral)[0]):
            worm_save_path = save_path+os.path.sep+str(expt_num).zfill(2)+os.path.sep+str(worm+first_worm_num).zfill(len(expt.acquired_wells[0]))
            worm_source_path = expt.expt_path+os.path.sep+str(worm+first_worm_num).zfill(len(expt.acquired_wells[0]))
            if not os.path.exists(worm_save_path): 
                os.mkdir(worm_save_path) # Make directory for worm
            
            # For each dev. stage, search through source_path for appropriate file and copy file to save path
            #for stg in range(1,5):
                 #Get the size at hatch
                #if stg == 1:
                    #img_fn_prefix_first = expt.expt_mdata['timepoints'][np.where(dev_labels[worm,:]==stg)[0][0]]
                    #img_fn = [img_fn for img_fn in list_image_files(worm_source_path) if img_fn_prefix_first in img_fn and match_string in img_fn][0]
                    #shutil.copyfile(worm_source_path+os.path.sep+img_fn, worm_save_path+os.path.sep+img_fn)
                
                 #Get the last frame
                #img_fn_prefix_last = expt.expt_mdata['timepoints'][np.where(dev_labels[worm,:]==stg)[0][-1]]
                #img_fn = [img_fn for img_fn in list_image_files(worm_source_path) if img_fn_prefix_last in img_fn and match_string in img_fn][0]
                #shutil.copyfile(worm_source_path+os.path.sep+img_fn, worm_save_path+os.path.sep+img_fn)
            
            # Grab only first stage (2nd image)
            #stg = 1
            #img_fn_prefix_first = expt.expt_mdata['timepoints'][np.where(dev_labels[worm,:]==stg)[0][2]]
            #img_fn = [img_fn for img_fn in list_image_files(worm_source_path) if img_fn_prefix_first in img_fn and match_string in img_fn][0]
            #shutil.copyfile(worm_source_path+os.path.sep+img_fn, worm_save_path+os.path.sep+img_fn)
            
            # Grab first egg lay
            eggtime_diffs = (np.array(expt.expt_mdata['timestamps'])-expt.expt_mdata['timestamps'][0])-expt_events['First Egg Laid'][worm]
            timestamp_idx = np.where(np.abs(eggtime_diffs) == np.abs(eggtime_diffs).min())[0][0]
            #print(expt.expt_mdata['timestamps'][timestamp_idx]-expt.expt_mdata['timestamps'][0])
            
            img_fn_prefix_first = expt.expt_mdata['timepoints'][timestamp_idx]
            #print(img_fn_prefix_first)
            #print(list_image_files(worm_source_path))
            img_fn = [img_fn for img_fn in list_image_files(worm_source_path) if img_fn_prefix_first in img_fn and match_string in img_fn][0]
            shutil.copyfile(worm_source_path+os.path.sep+img_fn, worm_save_path+os.path.sep+img_fn)

def colored_outlines(a_directory):
    '''
    For each outlined worm in a_directory, convert it to a grayscale 8-bit image with the outline in white and the rest of the image in black.
    '''
    for a_subdir in sorted(os.listdir(a_directory)):
        for an_image in sorted(os.listdir(a_directory + os.path.sep + a_subdir)):
            if 'bf_corrected' in an_image:
                filepath = a_directory + os.path.sep + a_subdir + os.path.sep + an_image
                print(filepath)
                my_image = freeimage.read(filepath)
                masked_conversion = fill_colored_outline(my_image)
                outline_conversion = fill_colored_outline(my_image, outline_only = True)
                freeimage.write(masked_conversion, filepath.replace('bf_corrected.tiff', 'hmask.png'))
                freeimage.write(outline_conversion, filepath.replace('bf_corrected.tiff', 'houtline.png'))
    return

def fill_colored_outline(an_image, outline_only = False):
    '''
    Fill out a colored outline and return a mask.
    '''
    # Get rid of the alpha channel if present.
    an_image = an_image.copy()
    if an_image.shape[2] == 4: an_image[:, :, 3] = 0
    
    # Get pixels with disproportionately more red.  
    my_mask = np.abs(np.max(an_image[:, :, 1:], axis = 2).astype('float64') - an_image[:, :, 0].astype('float64')).astype('uint16') > 0
    my_mask = zplib_image_mask.get_largest_object(my_mask)  
    
    # Fill holes.
    if not outline_only:        
        my_mask = zplib_image_mask.fill_small_area_holes(my_mask, 300000).astype('uint8')
    my_mask = my_mask.astype('uint8')
    my_mask[my_mask >= 1] = -1 
    return my_mask

def copy_humanmasks_for_dev_stages():
    # TODO
    pass
