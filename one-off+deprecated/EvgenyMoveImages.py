import pathlib, shutil
import numpy as np
import pandas as pd
from corral_annotations import annotation_file
import json

def move_pics(expt_dir,dest_parent_dir):
    expt_dir = pathlib.Path(expt_dir)
    dest_parent_dir = pathlib.Path(dest_parent_dir)
    
    if not (dest_parent_dir / expt_dir.parts[-1]).exists():
        print("Making directory: " + expt_dir.parts[-1])
        (dest_parent_dir / expt_dir.parts[-1]).mkdir()
    
    my_af = annotation_file.AnnotationFile(
        [filename 
        for filename in expt_dir.iterdir()
        if filename.is_file() and 'annotation' in str(filename).lower() 
        and '.tsv' in str(filename).lower()][0])
    mdata_fn = [filename for filename in expt_dir.iterdir()
        if 'experiment_metadata' in str(filename)][0]
    
    with mdata_fn.open('r') as mdata_fp:
        mdata = json.load(mdata_fp)
        mdata_timestamps = np.array(mdata['timestamps'])
    
    worm_dirs = sorted([dir_name for dir_name in expt_dir.iterdir()
        if dir_name.is_dir() and not dir_name.parts[-1][0].isalpha()])
    
    
    good_worms = my_af.data['Worm'][my_af.get_goodworms()]
    good_worm_idxs = {worm:worm_index for worm_index,worm in good_worms.iteritems()}
    print(good_worms)
    print(good_worm_idxs)
    
    # Patch to fix spe-9 datasets where development and adulthood were annotated separately
    if (~my_af.data['Death'].apply(lambda annotation: annotation.isnumeric())).all():
        my_af.data['Death'] = my_af.data['Death'].apply(
            lambda annotation: int(annotation[1:]) if len(annotation) > 0 else '')
        
        # Note: Need to mount to appropriate location and replace mount point here.
        dev_mddict = {
            '2016.02.16 spe-9 Run 9': '/media/drew/Untitled/20160216_spe9Acquisition',
            '2016.02.29 spe-9 Run 12A': '/media/drew/Untitled/20160229_spe9Acquisition_DevVarA',
            '2016.02.29 spe-9 Run 12B': '/media/drew/Untitled/20160229_spe9Acquisition_DevVarB',
            '2016.03.04 spe-9 Run 13B': '/media/drew/Untitled/20160304_spe9Acquisition_DevVarB',
            '2016.03.04 spe-9 Run 13C': '/media/drew/Untitled/20160304_spe9Acquisition_DevVarC',
        }
        
        # Estimate first egg laid timepoint by picking closest 3-hour resolution timepoint to that determined by the high-res developmental timepoints
        with (pathlib.Path(dev_mddict[expt_dir.parts[-1]]) / 'experiment_metadata.json').open('r') as mdata_fp:
            dev_mdata = json.load(mdata_fp)
            dev_mdata_timestamps = np.array(dev_mdata['timestamps'])
        #~ my_af.data['First Egg Laid'] = pd.Series([   # Works for all except the one time in 12A where animal didn't lay eggs
            #~ np.argmin(np.abs(mdata_timestamps - dev_mdata_timestamps[int(fe_timept)]))
            #~ if fe_timept is not '' else ''
            #~ for fe_timept in my_af.data['First Egg Laid']])
        def assign_egg_time(mdata_timestamps,dev_mdata_timestamps,fe_timept):
            if fe_timept is '': return ''
            elif fe_timept[0]=='W':
                return fe_timept[1:]
            else:
                return np.argmin(np.abs(mdata_timestamps - dev_mdata_timestamps[int(fe_timept)]))
        my_af.data['First Egg Laid'] = my_af.data['First Egg Laid'].apply(
            lambda fe_timept: assign_egg_time(mdata_timestamps,dev_mdata_timestamps,fe_timept))
            
    
    annotations_byframe = my_af.data_as_frames()

    for worm_dir in worm_dirs:
        frames_to_get = []
        print(('/'+worm_dir.parts[-1]))
        if ('/'+worm_dir.parts[-1]) in good_worm_idxs:
            worm_entry = annotations_byframe.iloc[good_worm_idxs[('/'+worm_dir.parts[-1])]]
            frames_to_get.append(worm_entry['Death'])
            
            adultspan = (mdata_timestamps[worm_entry['Death']]-mdata_timestamps[worm_entry['First Egg Laid']])/(3600*24) # Days
            for day_back in range(int(np.floor(adultspan))):
                frames_to_get.append(
                    np.argmin(np.abs(mdata_timestamps[worm_entry['Death']]-mdata_timestamps-(day_back+1)*(3600*24))))
            
            for quarter_back in range(3):   # Get the -6,-12,-18 frames
                frames_to_get.append(
                    np.argmin(np.abs(mdata_timestamps[worm_entry['Death']]-mdata_timestamps-(quarter_back+1)*(3600*6))))
            
            files_to_get = [
                worm_dir / (mdata['timepoints'][frame_num] + ' bf.png')
                for frame_num in frames_to_get]
            
            if not (dest_parent_dir / expt_dir.parts[-1] / worm_dir.parts[-1]).exists(): 
                (dest_parent_dir / expt_dir.parts[-1] / worm_dir.parts[-1]).mkdir()
            #~ [shutil.copy(str(my_file),       # Works in all cases except in a place where bad cleaning occurred (e.g. 12B '043')
                    #~ str(dest_parent_dir / expt_dir.parts[-1] / worm_dir.parts[-1] / my_file.parts[-1])) 
                    #~ for my_file in files_to_get]

            for my_file in files_to_get:
                if (dest_parent_dir / expt_dir.parts[-1] / worm_dir.parts[-1] / my_file.parts[-1]).exists():
                    continue
                
                if my_file.exists():
                    shutil.copy(str(my_file),
                        str(dest_parent_dir / expt_dir.parts[-1] / worm_dir.parts[-1] / my_file.parts[-1]))
                else:
                    print('Moved into life_after_death to get pic')
                    shutil.copy(str(my_file.parent / 'life_after_death' / my_file.parts[-1]), 
                        str(dest_parent_dir / expt_dir.parts[-1] / worm_dir.parts[-1] / my_file.parts[-1])) 
            
            #raise Exception() # For debugging.
            
            

'''
# Good for non-cleaned data.

good_worms = np.where(my_af.get_goodworms())[0]

annotations_byframe = my_af.data_as_frames()

for worm_dir,worm_entry in zip(worm_dirs,annotations_byframe.iterrows()):
    frames_to_get = []
    if int(worm_dir.parts[-1]) in good_worms:
        #print(worm_entry[1]['Death'])
        frames_to_get.append(worm_entry[1]['Death'])
        
        adultspan = (mdata_timestamps[worm_entry[1]['Death']]-mdata_timestamps[worm_entry[1]['First Egg Laid']])/(3600*24) # Days
        for day_back in range(int(np.floor(adultspan))):
            frames_to_get.append(
                np.argmin(np.abs(mdata_timestamps[worm_entry[1]['Death']]-mdata_timestamps-(day_back+1)*(3600*24))))
        
        for quarter_back in range(3):   # Get the -6,-12,-18 frames
            frames_to_get.append(
                np.argmin(np.abs(mdata_timestamps[worm_entry[1]['Death']]-mdata_timestamps-(quarter_back+1)*(3600*6))))
        
        files_to_get = [
            worm_dir / (mdata['timepoints'][frame_num] + ' bf.png')
            for frame_num in frames_to_get]
        z
        #print(files_to_get)
        
        if not (dest_parent_dir / worm_dir.parts[-1]).exists(): 
            (dest_parent_dir / worm_dir.parts[-1]).mkdir()
        [shutil.copy(str(my_file), 
            str(dest_parent_dir / worm_dir.parts[-1] / my_file.parts[-1])) 
            for my_file in files_to_get]
        
        #raise Exception() # For debugging.
'''


def move_masks(source_dir,dest_dir):
    
    '''
        source_dir - contains folders pooling all acquired wells
        dest_dir - directory where folder structure will be copied over
    '''
    
    source_dir = pathlib.Path(source_dir)
    dest_dir = pathlib.Path(dest_dir)
    
    if not (dest_dir).exists():
        print("Making directory: " + dest_dir)
        dest_dir.mkdir()
    
    mdata_jsons = []
    
    well_dirs = sorted([my_dir for my_dir in source_dir.iterdir()])
    num_wells = len(well_dirs)
        
    for well_num, well_dir in enumerate(well_dirs):
        
        print('On well {} of {}'.format(well_num,num_wells))
        with (well_dir / 'position_metadata_extended.json').open('r') as mdata_fp:
            mdata = json.load(mdata_fp)
        timestamps = np.array([mdata_item['timestamp'] for mdata_item in mdata])
        timepoints = [mdata_item['timepoint'] for mdata_item in mdata]
        egg_ages = [mdata_item['egg_age'] for mdata_item in mdata]
        ghost_ages = [mdata_item['ghost_age'] for mdata_item in mdata]
        
        # Find index corresponding to first egg laid and death
        egg_idx = np.argmin(np.abs(egg_ages))
        death_idx = np.argmin(np.abs(ghost_ages))
        
        frames_to_get = []
        frames_to_get.append(death_idx)
        
        adultspan = (timestamps[death_idx] - timestamps[egg_idx])/(3600*24)
        
        for day_back in range(int(np.floor(adultspan))):
            frames_to_get.append(
                np.argmin(np.abs(timestamps[death_idx] - timestamps - (day_back+1)*(3600*24))))
        
        for quarter_back in range(3):   # Get the -6,-12,-18 frames
            frames_to_get.append(
                np.argmin(np.abs(timestamps[death_idx] - timestamps - (quarter_back+1)*(3600*6))))
        
        files_to_get = [
            well_dir / (timepoints[frame_num] + ' mask.png')
            for frame_num in frames_to_get]
            
        files_to_get.extend([
            well_dir / (timepoints[frame_num] + ' bf.png')
            for frame_num in frames_to_get])
            
        if not (dest_dir / well_dir.parts[-1]).exists():
            (dest_dir / well_dir.parts[-1]).mkdir()
        
        for my_file in files_to_get:
            if (dest_dir / well_dir.parts[-1] / my_file.parts[-1]).exists():
                continue
            if my_file.exists():
                shutil.copy(str(my_file),
                    str(dest_dir / well_dir.parts[-1] / my_file.parts[-1]))
            else:
                print('Weird. Did not find file ' + str(my_file))
        
        
        ''' 
            This recreates the file structure similar to the original move_files
            TODO Finish this if desired
        #~ expt_prefix = ' '.join(well_dir.parts[-1].split()[:-1])
        #~ well_idx = well_dir.parts[-1].split()[-1]
        #~ if not (dest_dir / expt_prefix).exists():
            #~ (dest_dir / expt_prefix).mkdir()
        #~ if not (dest_dir / expt_prefix / well_idx).exists():
            #~ (dest_dir / expt_prefix / well_idx).mkdir()
        #~ for my_file in files_to_get:
            #~ if (dest_dir / expt_prefix
        '''

