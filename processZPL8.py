import matplotlib.pyplot as plt
#import MeasureFluorescence.wormFinding.backgroundSubtraction as bs
import wormFinding.backgroundSubtraction as bs

import os
import freeimage
import numpy as np
import datetime
import json
# import skimage.filters.rank as rank
import scipy.ndimage.filters as filters
from zplib.image import colorize
import zplib.image.mask

from collections import OrderedDict

import scipy.ndimage.measurements as measurements
import skimage.morphology
import pandas as pd
import zplib.curve.interpolate as zplib_interpolate

import re
import csv

def find_char(string,ch):
    return [idx for (idx, letter) in enumerate(string) if letter == ch]

def extract_datetime_fromstr(time_str):
    return datetime.datetime(int(time_str[0:4]), int(time_str[5:7]), int(time_str[8:10]), int(time_str[11:13]),
                             int(time_str[13:15]))

def extract_datetime_fromstr(time_str):
    return datetime.datetime.strptime(time_str,'%Y-%m-%dt%H%M')

def import_json_fromfn(fn):
    return json.loads(open(fn).read())
    
'''
Exporting to json file
    json.dump(variable,open('/media/Data/Work/ZPLab/WormImages/20151113_ZPL8Prelim3/07/info.json','w'))
    json.dump(metadata,open('./sam.json','w'),indent=4)
'''

def calibrate_images(data_path, cali_path, data_match_str='', cali_match_str='', save_path='', save_str=''):
    '''
    Calibrate images using flatfield calibration images from job runner
    data_path, cali_path - (string) paths to the relevant data or calibration images
    data_match_str, cali_match_str - (string) filter strings for selecting which data and flatfield images are used
    save_str - (string) tagged into the saved "corrected" image filename
    
    '''
    if len(save_path) == 0:
        save_path = data_path

    if len(save_str) == 0:
        save_str = 'corrected'
        
    # Extract the data files to be processed in the data directory and their corresponding calibration files
    organized_files = OrderedDict({})
    
    for d_file in sorted(os.listdir(data_path)):
        if data_match_str in d_file and os.path.isfile(data_path+os.path.sep+d_file):
            organized_files[d_file[0:15]]={'d_file':d_file}
    
    for c_file in os.listdir(cali_path):
        if cali_match_str in c_file and os.path.isfile(cali_path+os.path.sep+c_file):
            organized_files[c_file[0:15]]['c_file'] = c_file

    # Build a master vignette mask to screen out pixels excluded by vignetting over all time points
    master_v_mask = np.logical_and.reduce(np.array([freeimage.read(cali_path+os.path.sep+organized_files[dataset]['c_file'])>0 for dataset in organized_files.keys()]))

    # Grab the first image's intensity at a fixed pixel in the field
    [d_file_init, c_file_init] = [organized_files[next(iter(organized_files))]['d_file'],organized_files[next(iter(organized_files))]['c_file']]
    base_px_val = (freeimage.read(data_path+os.path.sep+d_file_init)*freeimage.read(cali_path+os.path.sep+c_file_init))[300,1200]
    print(d_file_init)
    print(c_file_init)
    print(base_px_val)

    # Note: Calibration image is written (float32) so that one MULTIPLIES it against the data image (vs. dividing it)
    for dataset in organized_files.keys():
        [d_file, c_file] = [organized_files[dataset]['d_file'], organized_files[dataset]['c_file']]
        trans_img = freeimage.read(data_path+os.path.sep+d_file)*freeimage.read(cali_path+os.path.sep+c_file)*master_v_mask
        if trans_img[300,1200]>base_px_val*1.8: 
            print(trans_img[300,1200])
            print(base_px_val*1.8)
            print(dataset)
            trans_img=trans_img/2
        save_fn = save_path+os.path.sep+d_file[:-4]+save_str+d_file[-4:]
        freeimage.write(trans_img.astype('uint16'), save_fn)

def make_composite_maskfile_batch(data_path, mask_path, save_path, data_str='', mask_str=''):
    data_fns = [data_f for data_f in sorted(os.listdir(data_path)) if data_str in data_f]
    
    print('importing data images')
    data_imgs = np.array([(freeimage.read(data_path + os.path.sep + data_f)) for data_f in data_fns])
    print('importing mask images')
    mask_imgs = np.array(
        [freeimage.read(mask_path + os.path.sep + mask_f) > 0 for data_f in data_fns for mask_f in sorted(os.listdir(mask_path)) 
         if data_f[0:15] in mask_f if mask_str in mask_f])
    try:
        os.stat(save_path)
    except:
        os.mkdir(save_path)
    
    print('generating and saving composites')
    comp = np.zeros(np.shape(data_imgs[[0]]))
    print('got here')
    for d_img,m_img,data_f in zip(data_imgs,mask_imgs,data_fns):
        comp=colorize.scale(np.copy(d_img),output_max=254)
        comp[m_img]=255
        #if data_f==data_fns[2]: return
        freeimage.write(comp.astype('uint8'), save_path+os.path.sep+data_f[:-4]+'composite'+data_f[-4:])

def extract_mask_fromcomposite_batch(comp_path, save_path, comp_str='', save_str=''):
    try:
        os.stat(save_path)
    except:
        os.mkdir(save_path)
    [freeimage.write(
        (freeimage.read(comp_path+os.path.sep+comp_f) == 255).astype('uint16')*65535,save_path+os.path.sep+comp_f[:-4]+save_str+comp_f[-4:]) \
        for comp_f in sorted(os.listdir(comp_path))]

def convert_to_8bit(data_path, save_path, filter_str = ''): #Scales based on the min and max OF EACH IMAGE!!! (NOT STANDARD ACROSS ALL); TODO EXTRA MODE FOR THIS)
    data_fns = [data_f for data_f in sorted(os.listdir(data_path)) if filter_str in data_f]
    data_imgs = np.array([freeimage.read(data_path+os.path.sep+data_f) for data_f in data_fns])
    for data_f, d_img in zip(data_fns, data_imgs):
        converted_img = colorize.scale(np.copy(d_img),min=d_img.min(),max=d_img.max(),output_max=255).astype('uint8')
        freeimage.write(converted_img, save_path+os.path.sep+data_f[:-4]+'_8bit'+data_f[-4:])

# From Willie's code
def fit_simple_skeleton(image_masks, save_dir = '', reverse_spine=np.array([])):
    def prune_skeleton(skeleton_graph, verbose_mode):
        '''
        Prune the graph of the skeleton so that only the single linear longest path remains.
        '''
        if verbose_mode:
            print('Pruning skeleton graph.')

        def farthest_node(my_graph, a_node, verbose_mode):
            '''
            Find the farthest node from a_node.
            '''
            reached_nodes = [a_node]
            distance_series = pd.Series(index = [str(a_node) for a_node in skeleton_graph.node_list])
            distance_series.loc[str(a_node)] = 0
            steps = 1
            current_circle = [a_node]
            next_circle = []
                            
            while len(reached_nodes) < len(my_graph.node_list):
                if verbose_mode:        
                    print('Reached nodes: ' + str(len(reached_nodes)))
                for current_node in current_circle:
                    next_steps = my_graph.edges(current_node)
                    if verbose_mode:
                        print('Current circle')
                        print(current_circle)
                        print('next_steps')
                        print(next_steps)
                    for one_step in next_steps:
                        other_node = [the_node for the_node in one_step if the_node != current_node][0]
                        if other_node not in reached_nodes:
                            distance_series.loc[str(other_node)] = steps 
                            next_circle.append(other_node)
                            reached_nodes.append(other_node)
                steps += 1
                current_circle = next_circle
                next_circle = []
            my_node = distance_series.argmax()
            my_node = (int(my_node[1:-1].split(',')[0]), int(my_node[1:-1].split(',')[1]))
            return my_node
        
        def find_minimal_path(my_graph, node_a, node_b, verbose_mode):
            '''
            Find the minimal path between node_a and node_b using my_graph.
            '''
            reached_nodes = [node_a]
            steps = 1
            current_circle = [[node_a]]
            next_circle = []
            got_to_b = False
            while not got_to_b:
                if verbose_mode:
                    print(len(reached_nodes))       
                    print([the_node for the_node in reached_nodes if the_node not in my_graph.node_list])       
                for current_path in current_circle:
                    current_node = current_path[-1]
                    next_steps = my_graph.edges(current_node)
                    for one_step in next_steps:
                        other_node = [the_node for the_node in one_step if the_node != current_node][0]
                        if other_node == node_b:
                            final_path = list(current_path)
                            final_path.append(other_node)
                            return (final_path, steps)
                            
                        elif other_node not in reached_nodes:
                            next_path = list(current_path)
                            next_path.append(other_node)
                            next_circle.append(next_path)
                            reached_nodes.append(other_node)
                steps += 1
                current_circle = next_circle
                next_circle = []    
            return
                            
        one_end = farthest_node(skeleton_graph, skeleton_graph.node_list[0], verbose_mode = verbose_mode)
        if verbose_mode:
            print('First end is: ' + str(one_end))
        other_end = farthest_node(skeleton_graph, one_end, verbose_mode = verbose_mode)
        if verbose_mode:
            print('Second end is: ' + str(other_end))
            
        (my_path, path_length) = find_minimal_path(skeleton_graph, one_end, other_end, verbose_mode = verbose_mode)
        my_path = np.array(my_path)
        return my_path

    def skeletonize_mask(raster_worm, verbose_mode):
        '''
        Given a masked worm in raster format, return a skeletonized version of it.
        '''
        if verbose_mode:
            print('Skeletonizing mask.')
        zero_one_mask = np.zeros(raster_worm.shape)
        zero_one_mask[raster_worm > 0] = 1
        zero_one_mask = zplib.image.mask.get_largest_object(zero_one_mask)
        my_skeleton = skimage.morphology.skeletonize(zero_one_mask)
        skeleton_mask = np.zeros(raster_worm.shape).astype('uint8')
        skeleton_mask[my_skeleton] = -1
        return skeleton_mask
            
    def skeleton_to_graph(skeleton_mask, verbose_mode):
        '''
        Converts a skeleton to a graph, which consists of a dictionary with a list of nodes (tuples containing the coordinates of each node) in 'nodes' and a list of edges (lists of two tuples containing coordinates of the nodes connected by the edge; all edges have length 1).
        '''
        if verbose_mode:
            print('Converting skeleton to graph.')
        node_list = [tuple(a_point) for a_point in np.transpose(np.array(np.where(skeleton_mask > 0)))]
        edge_list = []
        for point_a in node_list:
            for point_b in node_list:
                distance_vector = np.array(point_a) - np.array(point_b)
                check_distance = np.max(np.abs(distance_vector))
                my_edge = sorted([point_a, point_b])
                if check_distance == 1:
                    if my_edge not in edge_list:
                        edge_list.append(my_edge)
        
        class a_graph():
            def __init__(self, node_list, edge_list):
                self.node_list = node_list
                self.edge_list = edge_list
                return
            def edges(self, a_node):
                return [an_edge for an_edge in edge_list if a_node in an_edge]

        my_graph = a_graph(node_list, edge_list)
        return my_graph
    
    n_points = 10
    spine_out = np.zeros([(mask_imgs.shape)[0], n_points,2])
    
    if(reverse_spine.size==0):
        print('No orientation for spine found')
        reverse_spine=np.zeros([(mask_imgs.shape)[0]])
    else:
        print('Using supplied orientation')
    
    for (mask_idx, mask_img), reverse_this_spine in zip(enumerate(mask_imgs),reverse_spine):
        print('processing simple spine for mask img {:03}'.format(mask_idx))
        pruned_graph = prune_skeleton(skeleton_to_graph(skeletonize_mask(mask_img, True), True), True)  
        spine_tck = zplib_interpolate.fit_spline(pruned_graph)
        spine_out[mask_idx,:,:] = zplib_interpolate.spline_interpolate(spine_tck, n_points)
        if reverse_this_spine: spine_out[mask_idx,:,:] = np.flipud(spine_out[mask_idx,:,:])
        if save_dir != '':
            plt.figure(0)
            plt.gcf().clf()
            plt.imshow(mask_img.T)
            plt.scatter(spine_out[mask_idx,:,0],spine_out[mask_idx,:,1],c='b')
            plt.scatter(spine_out[mask_idx,0,0],spine_out[mask_idx,0,1],c='r')
            plt.figure(0).savefig(save_dir+os.path.sep+'mask_spine_{:03}.png'.format(mask_idx))
    
    return spine_out
    
def parse_acquisitions_log(expt_path):
    '''
        Simple format
    '''
    position_data = dict()
    timept_list = []
    
    with open(expt_path+os.path.sep+'acquisitions.log') as log_file:
        while True: # Look for next time point
            timept_str = log_file.readline()
            if timept_str == '': break  # EOF
            
            timept_id = (re.search(r'(?<=Starting timepoint ).{15}', timept_str)).group(0)
            timept_list.append(extract_datetime_fromstr(timept_id))
            position_data[timept_id] = dict()
            
            log_file.readline() # Discard next line
            while True: # Look for positions
                loc_line = log_file.readline()
                if 'ended' in loc_line: break
                position_line = log_file.readline()
                try:
                    loc_id = (re.search(r'(?<=Acquiring position )[0-9]*', loc_line)).group(0)
                    position_val = float((re.search(r'(?<=position: )[0-9]*[\.][0-9]*', position_line)).group(0))
                except AttributeError: # Can't find pattern (new version of acquisitions file)
                    #print(loc_line)
                    loc_id = (re.search(r'(?<=Acquiring Position: )[0-9]*', loc_line)).group(0)
                    position_val = float((re.search(r'(?<=z: )[0-9]*[\.][0-9]*', position_line)).group(0))
                
                position_data[timept_id][int(loc_id)] = position_val
                
    acquisition_times = np.array([(timept_event-timept_list[0]).total_seconds()/3600 for timept_event in timept_list])
    # Works for when every location was acquired at every time point
    #organized_z_loc = np.array([[position_data[timept_idx]['{:02}'.format(loc_idx)] for timept_idx in sorted(list(position_data.keys()))] for loc_idx in range(len(position_data[(list(position_data.keys()))[0]]))])
    temp_loc = [[] for idx in range(len(position_data[(list(position_data.keys()))[0]]))]
    for loc_idx in range(len(position_data[(list(position_data.keys()))[0]])):
        for timept_idx in sorted(list(position_data.keys())):
            if loc_idx in position_data[timept_idx].keys(): 
                temp_loc[loc_idx].append(position_data[timept_idx][int(loc_idx)])
    organized_z_loc = np.zeros([len(temp_loc),np.max(np.array([len(loc_data) for loc_data in temp_loc]))])
    for loc_idx in range(len(temp_loc)):
        organized_z_loc[loc_idx,0:len(temp_loc[loc_idx])] = temp_loc[loc_idx][:]
    return [acquisition_times, organized_z_loc]

def plot_autofocus_dryingcurve(expt_path, skip_positions=[]):
    [acq_times, z_locs] = parse_acquisitions_log(expt_path)
    
    print(acq_times)
    print(z_locs)
    print('Maximal drying between time points:{}'.format(np.max(np.abs(np.diff(z_locs,axis=1)))))
    print('Maximal drying during first day:{}'.format(np.max(np.abs(z_locs[:,(np.argwhere(acq_times>acq_times[0]+24.0))[0]] - z_locs[:,0]))))
    
    
    plt.gcf().clf()
    plt.ion()
    [plt.plot(acq_times[z_loc_set!=0],z_loc_set[z_loc_set!=0]) for z_loc_idx,z_loc_set in enumerate(z_locs) if z_loc_idx not in skip_positions]
    plt.xlabel('Time (from expt start; hr)')
    plt.ylabel('z_position (mm)')
    
    return (acq_times, z_locs)

'''
with open(annotation_file,newline='') as annotation_fp:
    bad_pos = csv.reader(annotation_fp,delimiter='\t')    [print(line) for line in bad_pos]
'''   

'''
with open(annotation_file,newline='') as annotation_fp:
    bad_pos = [idx-1 for idx,worm_info_line in enumerate(csv.reader(annotation_fp,delimiter='\t')) if 'OOF' in worm_info_line[5]]
'''

def process_directories(*data_dirs,mask_str='',fl_str='',save_dir=''):
    if(type(mask_str)==str): mask_str = [mask_str]*len(data_dirs)
    if(type(fl_str)==str): fl_str = [fl_str]*len(data_dirs)
    
    stat_fun = lambda data:np.percentile(data,95)
    fl_data = np.zeros([
        len(data_dirs),len([data_file for data_file in os.listdir(data_dirs[0]+os.path.sep+'corrected') if os.path.isfile(data_dirs[0]+os.path.sep+'corrected'+os.path.sep+data_file) and fl_str[0] in data_file])])
    life_events = []
    hatch_times = np.zeros([len(data_dirs)])
    
    time_strs = [time_str[0:15] for time_str in sorted(os.listdir(data_dirs[0]+os.path.sep+'masks')) if mask_str[0] in time_str and os.path.isfile(data_dirs[0]+os.path.sep+'masks'+os.path.sep+time_str)]
    times = [extract_datetime_fromstr(time_str) for time_str in time_strs]
    expt_start_t = times[0]
    times = np.array([(expt_time-expt_start_t).total_seconds() for expt_time in times])
    
    for dir_num,data_dir in enumerate(data_dirs):
        print('Processing directory '+str(dir_num))
        #Find time point where hatch occurs
        metadata_obj = import_json_fromfn(data_dir+os.path.sep+'metadata.json')
        life_events.append(metadata_obj['obs'])
        hatch_times[dir_num]=[(extract_datetime_fromstr(event_time)-expt_start_t).total_seconds() for [event_time,event_desc] in metadata_obj['obs'].items() if 'hatch' in event_desc][0]
        
        # Enumerate files and extract data
        fl_files = [data_dir+os.path.sep+'corrected'+os.path.sep+fl_file for fl_file in sorted(os.listdir(data_dir+os.path.sep+'corrected'+os.path.sep)) if fl_str[dir_num] in fl_file]
        m_files = [data_dir+os.path.sep+'masks'+os.path.sep+m_file for m_file in sorted(os.listdir(data_dir+os.path.sep+'masks'+os.path.sep)) if mask_str[dir_num] in m_file]
        
        for time_num,[m_file,fl_file] in enumerate(zip(m_files,fl_files)):
            fl_data[dir_num,time_num] = stat_fun((freeimage.read(fl_file))[freeimage.read(m_file)>0])
    
    # Subtract hatching time from each column of the matrix of accumulated developmental time points taken
    dev_times = ((np.tile(times,[len(data_dirs),1])).transpose()-np.array(hatch_times)).transpose()
    
    bins = np.arange(dev_times.min(),dev_times.max(),20*60) # ~20 mins/time pt
    binned_worm_data = [np.array([fl_data[dir_num,(dev_times[dir_num,:]>=bins[time_pt]) & (dev_times[dir_num,:]<bins[time_pt+1])] for time_pt in range(len(bins)-1)]) for dir_num in np.arange(len(data_dirs))]
    pop_stats_perbin = np.array([[np.mean(binned_worm_data[dir_num][time_pt]) for time_pt in range(len(bins)-1)] for dir_num in range(len(data_dirs))]) # num_worms x num_timepts
    pop_avg = np.nanmean(pop_stats_perbin,axis=0)   # Screen out nans due to empty bins
    
    pop_avg_ma = pop_avg-np.convolve(pop_avg,np.ones([5,])/5,mode='same')
    
    plt.figure(1)
    plt.gcf().clf()
    plt.ion()
    plt.show()

    # Plot stats per worm
    for dset_data in pop_stats_perbin:
        plt.plot((bins[:-1])[~np.isnan(dset_data)]/3600,dset_data[~np.isnan(dset_data)],'k')
    
    # Plot population average
    plt.plot((bins[:-1])[~np.isnan(pop_avg)]/3600,pop_avg[~np.isnan(pop_avg)],'r')

    plt.title('Aggregate raw fluorescence')
    plt.gca().set_xlabel('Time post-hatch (hr)')
    plt.gca().set_ylabel('Intensity (au)')
    
    if save_dir != '':
        plt.figure(1).savefig(save_dir+os.path.sep+'aggregate_fl_raw.png')
        
    plt.figure(2)
    plt.gcf().clf()
    plt.show()
    plt.plot((bins[:-1])[~np.isnan(pop_avg_ma)]/3600,pop_avg_ma[~np.isnan(pop_avg_ma)],'r')
    plt.title('Mean adjusted fluorescence')
    plt.gca().set_xlabel('Time post-hatch (hr)')
    plt.gca().set_ylabel('Intensity (au)')
    
    if save_dir != '':
        plt.figure(2).savefig(save_dir+os.path.sep+'aggregate_fl_ma.png')


if __name__ == "__main__":
    bs_match_string = 'focus-04'  # Fill this
    f_match_string = 'fluoro_z-10'  # Fill this
    data_path = '/media/Data/Work/ZPLab/WormImages/20151113_ZPL8Prelim3/07/corrected/'
    cali_path = '/media/Data/Work/ZPLab/WormImages/20151113_ZPL8Prelim3/calibrations/'
    do_calibration = False
    gen_masks = False
    save_figs = True
    save_notes = '20160106_devtrack_20151113_07'
    
    plot_raw_data = False
    plot_ma_data = False
    plot_centtrack = False
    plot_spinetrack = True
    gen_spine = False

    print("initializing")
    plt.ion()

    if do_calibration:
        # Assume that the raw data is in the data_path; save calibrated images in a separate directory in sbiling directory and change data_path accordingly
        try:
            os.stat(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'corrected'+os.path.sep)
        except:
            os.mkdir(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'corrected'+os.path.sep)

        calibrate_images(data_path, cali_path, bs_match_string, 'bf', data_path[:find_char(data_path,os.path.sep)[-2]+1]+'corrected'+os.path.sep)
        calibrate_images(data_path, cali_path, f_match_string, 'fl', data_path[:find_char(data_path,os.path.sep)[-2]+1]+'corrected'+os.path.sep)
        data_path = data_path[:find_char(data_path,os.path.sep)[-2]+1]+'corrected'+os.path.sep
        print('finished calibrating')
        print('changed data_path to:'+data_path)

    # masks in sibling directory 'masks'
    mask_path = data_path[:find_char(data_path,os.path.sep)[-2]+1]+'masks'+os.path.sep
    if gen_masks:
        print('generating masks in dir:'+mask_path)
        try:
            os.stat(mask_path)
        except:
            os.mkdir(mask_path)
        bs.overallBackgroundSubtract(data_path, bs_match_string, bs_match_string + 'mask', 10, mask_path)
        print('finished generating masks')
    
    print("loading image files")
    mask_imgs = np.asarray(
        [freeimage.read(mask_path + os.path.sep + mask_f) > 0 for mask_f in sorted(os.listdir(mask_path)) if \
         bs_match_string in mask_f if 'corrected' in mask_f])
    f_data_imgs = np.asarray(
        [freeimage.read(data_path + os.path.sep + data_f) for data_f in sorted(os.listdir(data_path)) if \
         f_match_string in data_f if 'corrected' in data_f])
    print("loaded "+str(np.shape(mask_imgs))+" image masks")

    marked_pixels = [f_img[m_img] for m_img, f_img in list(zip(mask_imgs, f_data_imgs))]

    stats = {0: np.mean, 1: lambda data: np.percentile(data, 95)}
    stats_strs = ['mean', '95th percentile']
    extracted_data = np.zeros((len(stats.keys()), np.shape(mask_imgs)[0]))
    extracted_data = np.asarray([[stats[row_key](pixels) for pixels in marked_pixels] for row_key in stats.keys()])

    # Get respective offsets for each timepoint relative to beginning of expt
    time_strs = [time_str[0:15] for time_str in sorted(os.listdir(mask_path)) if bs_match_string in time_str if
                 'corrected' in time_str]
    times = [extract_datetime_fromstr(time_str) for time_str in time_strs]
    time_offsets = [exp_time - times[0] for exp_time in times]
    #xax_ticklbls = ['{:3.1f}'.format(offset.total_seconds() / 3600) for offset in time_offsets]
    xax_times = [offset.total_seconds() / 3600 for offset in time_offsets]

    metadata_obj = import_json_fromfn(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'metadata.json')

    # Plot extracted data normally
    if plot_raw_data:
        plt.figure(1)
        plt.gcf().clf()
        stats_series = plt.plot(xax_times, extracted_data[0, :], '--r', xax_times, extracted_data[1, :], 'b')
            
        # Get metadata and grab observations
        for time_key, obs_value in metadata_obj['obs'].items(): # Plot events from recorded observations in metadata_obj
            plt.plot(np.array([(extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                (extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600]),
                np.array([extracted_data.min(), extracted_data.max()]),'k')
            plt.text((extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                extracted_data.max(), obs_value)

        plt.legend(stats_series, stats_strs,loc=2)
        plt.title('Raw tracked fluorescence')
        plt.gca().set_xlabel('Time from Expt Start (hr)')
        plt.gca().set_ylabel('Intensity (au)')
        plt.show()
        
        if save_figs:
            anal_path = data_path[:find_char(data_path,os.path.sep)[-2]+1]+'analysis'+os.path.sep
            try:
                os.stat(anal_path)
            except:
                os.mkdir(anal_path)
            plt.figure(1).savefig(anal_path+save_notes+'_fl_raw.png')
        
    
    # Plot signal after subtracting running mean (generated by convolving with a box filter)
    if plot_ma_data:
        plt.figure(2)
        plt.gcf().clf()
        stats_series = plt.plot(xax_times, extracted_data[0, :]-np.convolve(extracted_data[0,:],np.ones([5,])/5,mode='same'), '--r', \
            xax_times, extracted_data[1, :]-np.convolve(extracted_data[1,:],np.ones([5,])/5,mode='same'), 'b')
        metadata_obj = import_json_fromfn(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'metadata.json')
        for time_key, obs_value in metadata_obj['obs'].items(): # Plot events from recorded observations in metadata_obj
            plt.plot(np.array([(extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                (extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600]),
                np.array([extracted_data.min(), extracted_data.max()]),'k')
            plt.text((extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                extracted_data.max(), obs_value)
        plt.legend(stats_series, stats_strs,loc=2)
        plt.title('Mean subtracted fluorescence')
        plt.gca().set_xlabel('Time from Expt Start (hr)')
        plt.gca().set_ylabel('Intensity (au)')
        plt.show()
        if save_figs:
            anal_path = data_path[:find_char(data_path,os.path.sep)[-2]+1]+'analysis'+os.path.sep
            try:
                os.stat(anal_path)
            except:
                os.mkdir(anal_path)
            plt.figure(2).savefig(anal_path+save_notes+'_fl_ma.png')
        
    # Track centroid across image sequence
    if plot_centtrack:
        
        centroid_location = np.array([measurements.center_of_mass(m_img,m_img,index=1) for m_img in mask_imgs])
        plt.figure(3)
        plt.gcf().clf()
        plt.plot(xax_times[1:],np.linalg.norm(np.diff(centroid_location,axis=0),axis=1))      #Diff across time points, then take norms over columns
        
        metadata_obj = import_json_fromfn(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'metadata.json')
        for time_key, obs_value in metadata_obj['obs'].items(): # Plot events from recorded observations in metadata_obj
            plt.plot(np.array([(extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                (extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600]),
                np.array([extracted_data.min(), extracted_data.max()]),'k')
            plt.text((extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                extracted_data.max(), obs_value)
        plt.title('Tracked centroid')
        plt.gca().set_xlabel('Time from Expt Start (hr)')
        plt.gca().set_ylabel('Change in centroid position (px)')
        plt.show()    

        if save_figs:
            anal_path = data_path[:find_char(data_path,os.path.sep)[-2]+1]+'analysis'+os.path.sep
            try:
                os.stat(anal_path)
            except:
                os.mkdir(anal_path)
            plt.figure(3).savefig(anal_path+save_notes+'_centtrack.png')
    
    # Track spine location across image sequence
    if plot_spinetrack:
        # Import oriented spine if available by looking in metadata for orientation; otherwise, don't worry about it
        if gen_spine:
            try:    # Look for spine orientation
                spine_pts = fit_simple_skeleton(mask_imgs,save_dir=data_path[:find_char(data_path,os.path.sep)[-2]+1]+'spines'+os.path.sep,reverse_spine=np.array(metadata_obj['spine_needs_reversal']))
            except KeyError:
                spine_pts = fit_simple_skeleton(mask_imgs,save_dir=data_path[:find_char(data_path,os.path.sep)[-2]+1]+'spines'+os.path.sep)
            # Save spine in 'spines' directory
            spine_data = {variable:eval(variable) for variable in ['spine_pts']}
            json.dump(spine_data,open(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'spines'+os.path.sep+'spine.json','w'),indent=4)
        else:
            spine_pts = (import_json_fromfn(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'spines'+os.path.sep+'spine.json'))['spine_pts']
        spine_loc_diff = np.sum(np.sum(np.diff(spine_pts,axis=0)**2,axis=1)**0.5,axis=1) #sum RMSD (sum on distances b/t points) assumes num_txnum_ptsx2
        
        plt.figure(4)
        plt.gcf().clf()
        plt.plot(xax_times[1:],spine_loc_diff)
        
        metadata_obj = import_json_fromfn(data_path[:find_char(data_path,os.path.sep)[-2]+1]+'metadata.json')
        for time_key, obs_value in metadata_obj['obs'].items(): # Plot events from recorded observations in metadata_obj
            plt.plot(np.array   ([(extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                (extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600]),
                np.array([extracted_data.min(), extracted_data.max()]),'k')
            plt.text((extract_datetime_fromstr(time_key)-times[0]).total_seconds()/3600,
                extracted_data.max(), obs_value)
        plt.title('Tracked spine')
        plt.gca().set_xlabel('Time from Expt Start (hr)')
        plt.gca().set_ylabel('Change in centroid position (px)')
        plt.show()    

        if save_figs:
            anal_path = data_path[:find_char(data_path,os.path.sep)[-2]+1]+'analysis'+os.path.sep
            try:
                os.stat(anal_path)
            except:
                os.mkdir(anal_path)
            plt.figure(4).savefig(anal_path+save_notes+'_spinetrack.png')
