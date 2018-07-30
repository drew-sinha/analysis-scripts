import pathlib
import shutil
import json
import datetime
import numpy as np

from zplib import datafile
from elegant import load_data

def remove_offending_tp(experiment_root,timept_str,dry_run=False):
    '''Removes all of the files/entries for a given timepoint (backs them up in an associated "offending_files" directory)

    Parameters
        experiment_root - str/pathlib.Path pointing to experiment directory
        timept_str - str with standard format for timepoint acquisition (YYYY-MM-DDtHHMM)
        dry_run - bool flag that toggles taking action (if False, only specifies when offending files are found)
    '''

    if type(experiment_root) is str: experiment_root = pathlib.Path(experiment_root)

    for sub_dir in sorted(experiment_root.iterdir()):
        if sub_dir.is_dir():

            # Move bad files
            offending_files = [im_file for im_file in sub_dir.iterdir() if timept_str in str(im_file)]
            if len(offending_files)>0:
                print('Found offending files in: '+str(sub_dir))

                if not dry_run:
                    if not (sub_dir/'offending_files').exists(): (sub_dir/'offending_files').mkdir(mode=775)
                    [im_file.rename(im_file.parent/'offending_files'/im_file.name) for im_file in offending_files]

            # Check if position metadata is bad and handle
            if sub_dir.name.isnumeric():
                md_file = (sub_dir/'position_metadata.json')
                with md_file.open() as md_fp:
                    pos_md = json.load(md_fp)
                has_bad_md_tmpt = any([tmpt['timepoints'] == timept_str for tmpt in pos_md])
                if has_bad_md_tmpt:
                    print('Found offending entry in position_metadata in: '+str(sub_dir))
                    if not dry_run:
                        if not (sub_dir/'offending_files').exists(): (sub_dir/'offending_files').mkdir(mode=775)
                        md_file.rename(md_file.parent/'offending_files'/'position_metadata_old.json')     # Backup old
                        pos_md = [tp_data for tp_data in pos_md if tp_data['timepoints'] != timept_str]
                        datafile.json_encode_atomic_legible_to_file(pos_md,sub_dir/'position_metadata.json')   # Write out new position_metadata

            position_annotation_file = experiment_root / 'annotations' / f'{sub_dir.name}.pickle'
            if position_annotation_file.exists():
                general_annotations,timepoint_annotations = load_data.read_annotation_file(position_annotation_file)

                if timept_str in timepoint_annotations:
                    print(f'Found offending entry in annotation for position {sub_dir.name}')
                    if not dry_run:
                        del timepoint_annotations[timept_str]
                        load_data.write_annotation_file(position_annotation_file,
                            general_annotations, timepoint_annotations)


    md_file = (experiment_root/'experiment_metadata.json')
    with md_file.open() as md_fp:
        expt_md = json.load(md_fp)

    try:
        tp_idx = np.where(np.array(expt_md['timepoints']) == timept_str)[0][0]
        del expt_md['timepoints'][tp_idx]
        del expt_md['timestamps'][tp_idx]
        del expt_md['durations'][tp_idx]
        expt_md['brightfield_metering']={key:val for key,val in expt_md['brightfield_metering'] if key != timept_str}

        md_file.rename(md_file.parent/'experiment_metadata_old.json') #Backup
        datafile.json_encode_atomic_legible_to_file(expt_md,experiment_root/'experiment_metadata.json')
    except:  # Offending timepoint didn't make it into metadata
        pass

def remove_dead_timepoints(experiment_root, postmortem_time, delete_excluded=False):
    '''Deletes excess timepoints in an experiment where worms are dead

        Parameters
            experiment_root - str/pathlib.Path to experiment
            postmortem_time - Number of hours of data to keep past the annotated death timepoint;
                useful for keeping extra timepoints in case one ever wants to validate
                death for the previously made annotations
            delete_excluded - bool flag for whether to delete excluded positions; if True
                deletes folders for each position, but keeps the relevant annotation as a
                record of what happened at that position
    '''

    experiment_root = pathlib.Path(experiment_root)
    annotations = load_data.read_annotations(experiment_root)
    good_annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)

    if delete_excluded:
        excluded_positions = sorted(set(annotations.keys()).difference(set(good_annotations.keys())))
        for position in excluded_positions:
            if (experiment_root / position).exists():
                print(f'deleting excluded position {position}')
                shutil.rmtree(str(experiment_root / position))

    for position, position_annotations in good_annotations.items():
        general_annotations, timepoint_annotations = position_annotations
        timepoint_keys = list(timepoint_annotations.keys())
        timepoint_stages = [timepoint_info['stage'] for timepoint_info in timepoint_annotations.values()]
        death_timepoint = timepoint_keys[timepoint_stages.index('dead')]

        for timecourse_file in sorted((experiment_root / position).iterdir()):
            if timecourse_file.suffix[1:] not in ['png', 'tiff']: continue  # Also handles directories implicitly
            timepoint_label = timecourse_file.name.split(' ')[0]
            time_since_death = (_extract_datetime_fromstr(timepoint_label) - _extract_datetime_fromstr(death_timepoint))/3600
            if time_since_death > postmortem_time:
                timecourse_file.unlink()

        timepoints_to_delete = []
        for timepoint, timepoint_info in timepoint_annotations.items():
            if timepoint_info['stage'] == 'dead':
                time_since_death = (_extract_datetime_fromstr(timepoint) - _extract_datetime_fromstr(death_timepoint))/3600
                if time_since_death > postmortem_time:
                    timepoints_to_delete.append(timepoint)
        for dead_timepoint in timepoints_to_delete:
            del timepoint_annotations[dead_timepoint]
        load_data.write_annotation_file(experiment_root /'annotations' / f'{position}.pickle', general_annotations, timepoint_annotations)

        # TODO: Need to handle position_metadatas

def _extract_datetime_fromstr(time_str):
    '''Converts standard experimental timepoint string to time representation (seconds since epoch)'''
    return datetime.datetime.strptime(time_str,'%Y-%m-%dt%H%M').timestamp()
