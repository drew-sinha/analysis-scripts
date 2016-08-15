'''
    annotation_file.py
    Contains an interface for interacting with and extracting data out of annotation tsv files
'''

import numpy as np
import csv
import json
#from collections import OrderedDict
import os

def find_char(string,ch):
    '''
        Pulls out positions for a particular character in a string
    '''
    return [idx for (idx, letter) in enumerate(string) if letter == ch]

class AnnotationFile:
    
    def __init__(self, file_name, annotation_prefix=''):
        '''
            Arguments:
                file_name - file path to annotation file of interest
                annotation_prefix - character to add to front of annotation; useful for annotations coming from multiple experiments
        '''
        #self.data = OrderedDict()
        self.data = {}
        with open(file_name) as a_file:
            a_file_reader = csv.reader(a_file,delimiter='\t')
            a_tags = next(a_file_reader)
            for a_tag in a_tags:
                self.data[a_tag] = np.array([])
            for worm_data in a_file_reader:
                for a_tag,worm_tag_val in zip(a_tags,worm_data):
                    if annotation_prefix != '' and worm_tag_val != '' and a_tag not in ['Notes','Worm'] and not worm_tag_val[0].isalpha():   # If frame annotation field not labeled by experiment
                        self.data[a_tag] = np.append(self.data[a_tag],annotation_prefix+worm_tag_val)
                    else:
                        self.data[a_tag] = np.append(self.data[a_tag],worm_tag_val)
        
    def raw_data(self, expt_name='',restricted_list=None):
        raw_data = self.data.copy()
        if expt_name is not '':
            raw_data['Worm_FullName'] = np.array([expt_name+' '+worm_name[1:] for worm_name in raw_data['Worm']])
        if restricted_list is None:
            return raw_data
        else:
            return {a_tag:raw_data[a_tag][restricted_list] for a_tag in raw_data.keys()}
        
        
    def data_as_timestamps(self, metadata_list, expt_name='',restricted_list=None):
        '''
            metadata_list: (character delimiter: experiment metadata file paths)
            expt_name: Optional string identifying this experiment
            
            out_data: timestamps (-1 replaces any field that is left empty)
        '''
        out_data = self.data.copy()
        
        metadata_info = {}
        for metadata_pair in metadata_list.items():
            with open(metadata_pair[1]) as metadata_fp:
                metadata_info[metadata_pair[0]] = json.load(metadata_fp)
        
        for a_tag in out_data.keys():
            if a_tag not in ['Notes', 'Worm']:
                    out_data[a_tag] = \
                        np.array(
                        [metadata_info[time_idx[0]]['timestamps'][int(time_idx[1:])]-metadata_info[time_idx[0]]['timestamps'][0] if time_idx != '' else -1 for time_idx in out_data[a_tag]])
        
        if expt_name is not '':
            out_data['Worm_FullName'] = np.array([expt_name+' '+worm_name[1:] for worm_name in out_data['Worm']])
        
        if restricted_list is None:
            return out_data
        else:
            return {a_tag:out_data[a_tag][restricted_list] for a_tag in out_data.keys()}


    def data_as_timestamps_simple(self, metadata_file,expt_name='', restricted_list=None):
        '''
            Returns each set of data in the 'data' field as times relative to the start of the experiment (IN SECONDS)
            
            ARGUMENTS:
                metadata_file (string) - Corresponding metadata file to use that relates each frame number from the annotation to an actual data/time
                restricted_list (numpy array suitable for indexing) - list for grabbing data from a specific set of worms
        '''
        
        out_data = self.data.copy()
        
        with open(metadata_file) as metadata_fp:
            metadata_info = json.load(metadata_fp)
        
        # For each column in tsv,
        for a_tag in out_data.keys():
            out_data[a_tag][out_data[a_tag]=='']='-1'   # Empty field gets '-1'
            if a_tag != 'Notes' and a_tag != 'Worm':    # If not the non-time fields
                out_data[a_tag]=np.array([metadata_info['timestamps'][int(time_idx)]-metadata_info['timestamps'][0] for time_idx in out_data[a_tag]])
                out_data[a_tag][self.data[a_tag]=='']=-1
        
        if expt_name is not '':
            out_data['Worm_FullName'] = np.array([expt_name+' '+worm_name[1:] for worm_name in out_data['Worm']])

        if restricted_list is None:
            return out_data
        else:
            return {a_tag:out_data[a_tag][restricted_list] for a_tag in out_data.keys()}
    
    def get_goodworms(self, bad_worm_kws=[], restrict_to_hatched=False, expt_path=None):
        #if bad_worm_kws is []:
            #bad_worm_kws = ['FERTIL', 'Nh', 'NO HATCH', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS'] # Note: first item should be FERTIL!
        #viable_worm = (self.data['Hatch']!='') \
            #& (self.data['Death']!='') \
            #& np.array([not any([kw in note for kw in bad_worm_kws]) for note in self.data['Notes']])
        
        if len(bad_worm_kws) is 0:   # Use DEAD as marker
            viable_worm = (self.data['Hatch']!='') \
                & (self.data['Death']!='') \
                & np.array([('DEAD' in note and not 'NEVER LAID EGGS' in note) for note in self.data['Notes']])
        else:
            viable_worm = (self.data['Hatch']!='') \
                & (self.data['Death']!='') \
                & np.array([not any([kw in note for kw in bad_worm_kws]) for note in self.data['Notes']])
        
        goodworms = viable_worm #& worm_was_acquired
        if expt_path is not None:   # Need to screen for wells that weren't acquired (i.e. deleted)
            first_worm_num = int(self.data['Worm'][0][-2:])   # Add in adjustment for the first index of the worms not being at 0
            
            expt = Experiment(expt_path)
            worm_was_acquired = [str(worm_num+first_worm_num).zfill(len(expt.get_acquired_wells()[0])) in expt.get_acquired_wells() for worm_num in np.arange(len(self.data['Hatch']))]
        
        if restrict_to_hatched: goodworms = goodworms & (self.data['Hatch'] !=0)
        
        #if return_worm_data:
            #return (goodworms, {a_tag:self.data[a_tag][goodworms] for a_tag in data.keys()})
        #else:
        return goodworms
    
    def save_timestamp_tsv(self, metadata_file,output_file):
        with open(output_file,'w') as output_fp:
            output_writer = csv.writer(output_fp,delimiter='\t')
            output_writer.writerow([a_tag for a_tag in self.data.keys()])
            for worm_data in zip(*(self.data_as_timestamps(metadata_file).values())):
                output_writer.writerow(worm_data)

def compile_expt_timestamped_data(expt_dirs, md_dict=None):
    timestamped_data = {}
    if md_dict is None:
        for expt_dir in expt_dirs: 
            ann_file = AnnotationFile(
                [expt_dir+os.path.sep+my_file for my_file in sorted(os.listdir(expt_dir)) if '.tsv' in my_file][0])
            #if not any(timestamped_data):
                #[timestamped_data.setdefault(expt_key,np.array([])) for expt_key in ann_file.data.keys()]
            if expt_dir[-1] is not os.path.sep:
                expt_name = expt_dir[find_char(expt_dir,os.path.sep)[-1]:]
            else:
                expt_name = expt_dir[find_char(expt_dir,os.path.sep)[-2]:-1]
            ann_file_data = ann_file.data_as_timestamps_simple(expt_dir+os.path.sep+'experiment_metadata.json',restricted_list=ann_file.get_goodworms(),expt_name=expt_name)
            #print(timestamped_data)
            #print(ann_file.data['Notes'][ann_file.get_goodworms()])
            if not any(timestamped_data):
                [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in ann_file_data.keys()]
            timestamped_data = {expt_key:np.append(timestamped_data[expt_key],ann_file_data[expt_key]) if expt_key in ann_file_data.keys() else np.append(timestamped_data[expt_key],[-1]*np.count_nonzero(ann_file.get_goodworms())) for expt_key in timestamped_data.keys()}
    else:
        for expt_dir, md_map in zip(expt_dirs, md_dict):
            if expt_dir[-1] is not os.path.sep:
                expt_name = expt_dir[find_char(expt_dir,os.path.sep)[-1]:]
            else:
                expt_name = expt_dir[find_char(expt_dir,os.path.sep)[-2]:-1]

            if list(md_map.keys())[0] is '':
                ann_file = AnnotationFile(
                    [expt_dir+os.path.sep+my_file for my_file in sorted(os.listdir(expt_dir)) if '.tsv' in my_file and 'lock' not in my_file][0])
                ann_file_data = ann_file.data_as_timestamps_simple(expt_dir+os.path.sep+'experiment_metadata.json',restricted_list=ann_file.get_goodworms(),expt_name=expt_name)
            else:
                ann_file = AnnotationFile(
                    [expt_dir+os.path.sep+my_file for my_file in sorted(os.listdir(expt_dir)) if '.tsv' in my_file and 'lock' not in my_file][0],
                    #annotation_prefix=list(md_map.keys())[0])
                    annotation_prefix='D')
                #if list(md_map.values())[0] is '':
                    #ann_file_data = ann_file.data_as_timestamps({list(md_map.keys())[0]:expt_dir+os.path.sep+'experiment_metadata.json'},restricted_list=ann_file.get_goodworms())
                #else:
                ann_file_data = ann_file.data_as_timestamps(md_map,restricted_list=ann_file.get_goodworms(),expt_name=expt_name)

            if not any(timestamped_data):
                [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in ann_file_data.keys()]
            timestamped_data = {expt_key:np.append(timestamped_data[expt_key],ann_file_data[expt_key]) if expt_key in ann_file_data.keys() else np.append(timestamped_data[expt_key],[-1]*np.count_nonzero(ann_file.get_goodworms())) for expt_key in timestamped_data.keys()}
            #print(timestamped_data)
    
    return timestamped_data

def compile_expt_raw_data(expt_dirs):
    raw_data = {}
    for expt_dir in expt_dirs:
        ann_file = AnnotationFile(
            [expt_dir+os.path.sep+my_file for my_file in sorted(os.listdir(expt_dir)) if '.tsv' in my_file][0])
        if expt_dir[-1] is not os.path.sep:
            expt_name = expt_dir[find_char(expt_dir,os.path.sep)[-1]:]
        else:
            expt_name = expt_dir[find_char(expt_dir,os.path.sep)[-2]:-1]
        # Get data
        ann_file_data = ann_file.raw_data(expt_name=expt_name,restricted_list=ann_file.get_goodworms())
        if not any(raw_data):
            [raw_data.setdefault(expt_key,np.array([])) for expt_key in ann_file_data.keys()]
        raw_data = {expt_key:np.append(raw_data[expt_key],ann_file_data[expt_key]) if expt_key in ann_file_data.keys() else np.append(raw_data[expt_key],[-1]*np.count_nonzero(ann_file.get_goodworms())) for expt_key in raw_data.keys()}
    return raw_data
