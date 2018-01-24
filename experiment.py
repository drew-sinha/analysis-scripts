import freeimage
import numpy as np
import json
import pathlib
import pickle
import corral_annotations.annotation_file as annotation_file

def import_json_fromfn(fn):
    return json.loads(open(str(fn)).read())


class Experiment:
    BF_REF_INTENSITY = 11701.7207031
    
    def __init__(self, expt_path):
        self.expt_path = pathlib.Path(expt_path)
        self.expt_mdata = import_json_fromfn(self.expt_path / 'experiment_metadata.json')
        
        self.acquired_wells = [well_dir.name for well_dir in self.expt_path.iterdir()
            if well_dir.is_dir() and not well_dir.name[0].isalpha()]
        
        self.position_mdata = {}
        for well in self.acquired_wells:
            self.position_mdata[well] = import_json_fromfn(self.expt_path / well / 'position_metadata.json')
        tsv_files = [my_file for my_file in self.expt_path.iterdir() if my_file.suffix ==  '.tsv']
        if len(tsv_files) > 0:
            self.expt_af = annotation_file.AnnotationFile(tsv_files[0])
        else:
            self.expt_af = None
    
    def get_calibrated_image(self, well, timepoint, label, calibration_type, apply_vignette = True):
        '''
            well - string for well
            timepoint - string for timepoint per YYYY-MM-DDtHHMM format
            label - string suffix coming after timepoint
            calibration_type - ['bf,'fl']
        '''
        image_fn = self.expt_path / well / (
            timepoint + ' ' + label + '.png')
        calibration_fn = self.expt_path / 'calibrations' / (
            timepoint + ' ' + calibration_type + '_flatfield.tiff')
        
        img = freeimage.read(str(image_fn))
        calibration_img = freeimage.read(str(calibration_fn))
        ref_intensity = self.expt_mdata['brightfield metering'][timepoint]['ref_intensity']
        
        if (self.expt_path / 'super_vignette.pickle').exists() and apply_vignette:
            with (self.expt_path / 'super_vignette.pickle').open('rb') as sv_file:
                sv = pickle.load(sv_file)
            calibrated_img = (img * calibration_img / ref_intensity * self.BF_REF_INTENSITY).astype('uint16')
            calibrated_img[~sv] = 0
            return calibrated_img
        else:
            return (img * calibration_img / ref_intensity * self.BF_REF_INTENSITY).astype('uint16')
    
    def get_raw_image(self, well, timepoint, label):
        image_fn = self.expt_path / well / (
            timepoint + ' ' + label + '.png')
        return freeimage.read(str(image_fn))
    
    def get_super_vignette(self):
        if (self.expt_path / 'super_vignette.pickle').exists():
            with (self.expt_path / 'super_vignette.pickle').open('rb') as sv_file:
                sv = pickle.load(sv_file)
            return sv
        else:
            raise FileNotFoundError('No super vignette found for this experiment.')
    
    def get_af_data(self,restrict_good = True):
        return self.expt_af.data_as_timestamps_simple(
        self.expt_path / 'experiment_metadata.json',
            restricted_list = self.expt_af.get_goodworms(bad_worm_kws = ['BURST']))
