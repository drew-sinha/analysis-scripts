import pathlib
import shutil
import json
import numpy as np

from zplib import datafile

def remove_offending_tp(expt_path,timept_str,dry_run=False):
    '''
        expt_path - str/pathlib.Path pointing to expt. dir
        timept_str - str in format yyyymmdd-thhmm
        dry_run - Toggles taking action (if False, will not delete but will verbosely specify where offending files are found
    '''
    
    if type(expt_path) is str: expt_path = pathlib.Path(expt_path)
    
    for sub_dir in expt_path.iterdir():
        if sub_dir.is_dir():
            
            # Move bad files
            offending_files = [im_file for im_file in sub_dir.iterdir() if timept_str in str(im_file)]
            if len(offending_files)>0:
                print('Found offending files in: '+str(sub_dir))
                
                if not dry_run:
                    if not (sub_dir/'offending_files').exists(): (sub_dir/'offending_files').mkdir(mode=755)
                    [im_file.rename(im_file.parent/'offending_files'/im_file.parts[-1]) for im_file in offending_files]
            
            # Check if position metadata is bad and handle 
            if sub_dir.parts[-1] != 'calibrations':
                md_file = (sub_dir/'position_metadata.json')
                if not md_file.exists(): print('No metadata file found at:'+str(sub_dir))
                else:
                    with md_file.open() as md_fp:
                        pos_md = json.load(md_fp)
                    has_bad_md_tmpt = any([tmpt['timepoints'] == timept_str for tmpt in pos_md])
                    if has_bad_md_tmpt:
                        print('Found offending entry in position_metadata in: '+str(sub_dir))
                        if not dry_run:
                            if not (sub_dir/'offending_files').exists(): (sub_dir/'offending_files').mkdir(mode=755)
                            md_file.rename(md_file.parent/'offending_files'/'position_metadata_old.json')     # Backup old
                            pos_md = [tp_data for tp_data in pos_md if tp_data['timepoints'] != timept_str]
                            with (sub_dir/'position_metadata.json').open('w') as md_fp:
                                encode_legible_to_file(pos_md,md_fp)   # Write out new position_metadata
            
            if (expt_dir / 'annotations' / sub_dir.name):
                position_annotations,timepoint_annotations = load_data.read_annotation_file(expt_dir / 'annotations' / sub_dir.name)

                if timept_str in timepoint_annotations:
                    print(f'Found annotation for position {sub_dir.name}')
                    if not dry_run:
                        timepoint_annotations.pop(timept_str)
                        load_data.write_annotation_file(expt_dir / 'annotations' / sub_dir.name, 
                            position_annotations, timepoint_annotations)
                
    
    md_file = (expt_path/'experiment_metadata.json')
    with md_file.open() as md_fp:
        expt_md = json.load(md_fp)
    
    try:
        tp_idx = np.where(np.array(expt_md['timepoints']) == timept_str)[0][0]
        del expt_md['timepoints'][tp_idx]
        del expt_md['timestamps'][tp_idx]
        del expt_md['durations'][tp_idx]
        expt_md['brightfield_metering']={key:val for key,val in expt_md['brightfield_metering'] if key != timept_str}
        
        md_file.rename(md_file.parent/'experiment_metadata_old.json') #Backup
        with (expt_path/'experiment_metadata.json').open('w') as md_fp:  # Write out new
            datafile.encode_atomic_legible_to_file(expt_md,md_fp)
    except:  # Offending timepoint didn't make it into metadata
        pass
