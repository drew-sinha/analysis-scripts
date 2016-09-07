import pathlib
import csv
import numpy as np

########
# Load tsv into dict
def load_tsv_as_dict(input_file, convert_strs_to_nums=True):
    data={}
    with input_file.open() as input_fp:
        input_reader = csv.reader(input_fp, delimiter='\t')
        fields = next(input_reader)
        data_field_types = []
        for a_field in fields:
            data[a_field] = []
        for worm_data in input_reader:
            if convert_strs_to_nums:
                if len(data_field_types) is 0:  # Run on first iteration
                    for data_val_str in worm_data:
                        if any([c.isalpha() for c in data_val_str]):
                            data_field_types.append(str)
                        elif '.' in data_val_str:
                            data_field_types.append(np.float64)
                        else:
                            data_field_types.append(np.int16)
                    print(data_field_types)
                for data_field, data_val,data_type in zip(fields,worm_data,data_field_types):
                    data[data_field].append(data_type(data_val))
            else:
                for data_field, data_val in zip(fields,worm_data):
                    data[data_field].append(data_val)
    return data

def sort_dict_by_field(my_dict, my_field):
    #def argsort(my_iter):
        #return sorted(range(len(my_iter)), key=my_iter.__getitem__)
    
    #field_idxs = argsort(my_dict[my_field])
    np_dict = {a_field:np.array(a_valueset) for a_field,a_valueset in my_dict.items()}
    field_idxs = np.argsort(np_dict[my_field])
    
    return {a_field:a_valueset[field_idxs] for a_field,a_valueset in np_dict.items()}

input_filepath = '/media/Data/Work/ZPLab/Analysis/WormDeath/deathdataForEmily_v2_allgood.tsv'   
image_dirpath = '/media/Data/Work/ZPLab/Analysis/WormDeath/EmilyDeathPics/'

input_file = pathlib.Path(input_filepath)
worm_deathdata = load_tsv_as_dict(input_file)       # Need to sort by lifespan!!!
worm_deathdata = sort_dict_by_field(worm_deathdata,'AnnotatedLifespan')

img_dir = pathlib.Path(image_dirpath)
#rw.flipbook.add_image_files([img_dir/(worm_data['Worm']+' '+worm_data['AnnotatedDeathTimepoint']+' bf.png') for worm_data in worm_deathdata[0:19]])
#image_list = [img_dir/(worm_name+' '+worm_deathtp+' bf.png') for worm_name,worm_deathtp in zip(worm_deathdata['Worm'],worm_deathdata['AnnotatedDeathTimepoint'])]
image_list = [img_dir / (worm_name[1:]+'_'+worm_deathtp+' bf.png') for worm_name,worm_deathtp in zip(worm_deathdata['Worm'],worm_deathdata['AnnotatedDeathTimepoint'])]
rw.flipbook.add_image_files(image_list[:20])
rw.flipbook.add_image_files(image_list[-20:])
