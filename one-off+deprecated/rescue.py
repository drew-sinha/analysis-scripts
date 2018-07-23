import shutil
import pathlib
import zplib.datafile as datafile
import json
import os

def backup_metadata(experiment_dir, backup_dir):
	experiment_dir = pathlib.Path(experiment_dir)
	backup_dir = pathlib.Path(backup_dir)
	for position_dir in sorted(list(experiment_dir.iterdir())):
		if position_dir.is_dir() and position_dir.parts[-1][0].isnumeric():
			if not(backup_dir/position_dir.parts[-1]).exists(): os.mkdir(str(backup_dir/position_dir.parts[-1]))
			shutil.copyfile(str(position_dir/'position_metadata.json'), str(backup_dir/position_dir.parts[-1]/'position_metadata.json'))

def prune_z(experiment_dir):
    '''
    For when I needed to rescue some scripts after Zach's update to new code
    '''
    
	experiment_dir = pathlib.Path(experiment_dir)
	for position_dir in sorted(list(experiment_dir.iterdir())):
		if position_dir.is_dir() and position_dir.parts[-1][0].isnumeric():
			shutil.copyfile(str(position_dir/'position_metadata.json'), str(position_dir/'position_metadata_old.json'))
			position_metadata = json.load((position_dir/'position_metadata.json').open('r'))
			if type(position_metadata[-1]['fine_z']) is list:
				position_metadata[-1]['fine_z'] = position_metadata[-1]['fine_z'][0]
			datafile.json_encode_atomic_legible_to_file(position_metadata, (position_dir/'position_metadata.json'))


# From debug_timepoints
#~ def reset_mds(expt_path):
    #~ if type(expt_path) is str: expt_path = pathlib.Path(expt_path)  
    #~ for sub_dir in expt_path.iterdir():
        #~ if sub_dir.is_dir() and (sub_dir/'position_metadata_old.json').exists():
            #~ (sub_dir/'position_metadata_old.json').rename(sub_dir/'position_metadata.json')
            #~ (sub_dir/'position_metadata_old.json').unlink()

"""
Some debugging for to check consistency of multiframe movement acquistion and imputation
"""

def check_expt_movement_acquisitions(expt_dir, mode='absolute'):
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)
    expt_pos_data, dirs_with_mdata, dirs_miss_mdata = load_movement_metadata(expt_dir)

    expt_mdata = json.load((expt_dir/'experiment_metadata.json').open())
    if mode == 'absolute':  # Look at the absolut time at which a frame was taken (with 0 being the first bright field image)
        unpacked_times = [[[item[i] for item in expt_pos_data['movement_frame_times'][timepoint]] for timepoint in expt_mdata['timepoints']] for i in range(5)] # unpacked_times [=] [pic_number, timepoint, well_number]

        fig_h,ax_h = plt.subplots(1,5)
        [[movement_ax.scatter(idx*np.ones([len(timepoint_times),1]),np.array(timepoint_times)) for idx,timepoint_times in enumerate(frame_times)] for movement_ax,frame_times in zip(ax_h,unpacked_times)] 
        [movement_ax.set_title('Movement Frame {}'.format(n+1)) for n,movement_ax in enumerate(ax_h)]
        [movement_ax.set_xlabel('Timepoint Index') for movement_ax in ax_h]
        [movement_ax.set_ylabel('Time Taken (s)') for movement_ax in ax_h]
    elif mode == 'deltas': # Look at time taken between each set of consecutive frames in an acquisition
        expt_pos_data_diff = {timepoint:[np.diff(item) for item in expt_pos_data['movement_frame_times'][timepoint]] for timepoint in expt_mdata['timepoints']} #{timepoint: list of lists}
        unpacked_diffs = [[[item[i] for item in expt_pos_data_diff[timepoint]] for timepoint in expt_mdata['timepoints']] for i in range(4)] # unpacked_times [=] [pic_number, timepoint, well_number]
        
        fig_h,ax_h = plt.subplots(1,4)
        [[movement_ax.scatter(idx*np.ones([len(timepoint_times),1]),np.array(timepoint_times)) for idx,timepoint_times in enumerate(frame_diffs)] for movement_ax,frame_diffs in zip(ax_h,unpacked_diffs)] 
        [movement_ax.set_title('Interval b/t Frames {} & {}'.format(n+1,n+2)) for n,movement_ax in enumerate(ax_h)]
        [movement_ax.set_xlabel('Timepoint Index') for movement_ax in ax_h]
        [movement_ax.set_ylabel('Time Taken (s)') for movement_ax in ax_h]

def copy_position_md(expt_dir, dest_dir):
    '''
        Copies position metadatas all to one destination folder ('dest_dir') for backing up
    '''
    
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)
    if type(dest_dir) is str: dest_dir = pathlib.Path(dest_dir)
    for sub_dir in expt_dir.iterdir():
        if sub_dir.is_dir() and str(sub_dir.parts[-1]) != 'calibrations':
            if not (dest_dir/sub_dir.parts[-1]).exists(): (dest_dir/sub_dir.parts[-1]).mkdir(mode=744)
            shutil.copyfile(str(sub_dir/'position_metadata.json'),str((dest_dir/sub_dir.parts[-1])/'position_metadata.json'))


def impute_missing_metadata_WZacquisition(expt_dir, dry_run = False):
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)
    expt_pos_data, dirs_with_mdata, dirs_miss_mdata = load_movement_metadata(expt_dir)  #(expt_pos_data = {data_key: timepoint: list of list of values})
    
    expt_mdata = json.load((expt_dir/'experiment_metadata.json').open())
    unpacked_times_tmpt = [[[item[pic_num] for item in expt_pos_data['movement_frame_times'][timepoint]] for pic_num in range(5)] for timepoint in expt_mdata['timepoints']]# unpacked_times_tmpt [=] [timepoint,pic_number,well_number]
    
    entry_dict = {
        'fine_z':None,  # Not needed
        'image_timestamps': {
            'bf.png': 0.0
        },
        'movement_frame_times': None,
        'timepoint':None,
        'timepoints':None,
        'timestamp':None,
        'notes':'added by automatic imputation (DS)',
    }
    # WARNING!!!!!! TODO Make this scalable to arbitrary position metadata....
    
    new_md_data = []
    for timepoint,timestamp,timepoint_times in zip(expt_mdata['timepoints'],expt_mdata['timestamps'],unpacked_times_tmpt):
        new_entry = entry_dict.copy()
        new_entry['timepoint'] = new_entry['timepoints'] = timepoint
        new_entry['timestamp'] = timestamp
        new_entry['movement_frame_times'] =[np.median(frame_times) for frame_times in timepoint_times]
        new_md_data.append(new_entry)
    
    if dry_run:
        print('New metadata contents:')
        print(new_md_data)
        
        print('Directories missing metadata')
        print(dirs_miss_mdata)
    else:
        for pos_dir in dirs_miss_mdata:
            print('Writing new metadata to directory: '+str(pos_dir))
            with (pos_dir/'position_metadata.json').open('w') as md_fp:
                encode_legible_to_file(new_md_data,md_fp)

def load_movement_metadata(expt_dir):
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)

    # Load experiment_metadata and grab timepoints to populate
    expt_mdata = json.load((expt_dir/'experiment_metadata.json').open())


    # Find dirs with position_metadata
    dirs_with_mdata = []
    dirs_miss_mdata = []
    for pos_dir in expt_dir.iterdir():
        if (pos_dir/'position_metadata.json').exists(): dirs_with_mdata.append(pos_dir)
        elif pos_dir.is_dir() and str(pos_dir.parts[-1]) != 'calibrations': dirs_miss_mdata.append(pos_dir)
    
    pos_mdata = json.load((dirs_with_mdata[0] / 'position_metadata.json').open())
    mdata_keys = [mdata_key for mdata_key in pos_mdata[0].keys()]
    expt_pos_data = {mdata_key:{timepoint:[] for timepoint in expt_mdata['timepoints']} \
        for mdata_key in mdata_keys}    #(expt_pos_data = {data_key: timepoint: list of values})

    for pos_dir in dirs_with_mdata:
        pos_mdata = json.load((pos_dir/'position_metadata.json').open())
        for item in pos_mdata:
            for mdata_key in mdata_keys:
                expt_pos_data[mdata_key][item['timepoint']].append(item[mdata_key])
    
    return (expt_pos_data, dirs_with_mdata, dirs_miss_mdata)

