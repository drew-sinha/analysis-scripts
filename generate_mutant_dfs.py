from basicOperations import folderStuff
from analyzeHealth import characterizeTrajectories
import pickle
import os

working_directory = r'/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/work_dir/'
human_directory = r'/media/Data/Work/ZPLab/Analysis/MutantHealth/worm_health_data/utilities/'
#working_directory = '/home/dsinha/work_dir/'
#human_directory = '/home/dsinha/utilities/'

data_list = {   
    'spe-9': {
        'data_directories':[
            r'/mnt/iscopearray/Zhang_William/2016.02.16 spe-9 Run 9',      #0
            r'/mnt/iscopearray/Zhang_William/2016.02.20 spe-9 Run 10A',    #1
            r'/mnt/iscopearray/Zhang_William/2016.02.20 spe-9 Run 10B',    #2
            r'/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11A',    #3
            r'/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11B',    #4
            r'/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11C',    #5
            r'/mnt/iscopearray/Zhang_William/2016.02.26 spe-9 Run 11D',    #6
            r'/mnt/iscopearray/Zhang_William/2016.02.29 spe-9 Run 12A',    #7
            r'/mnt/iscopearray/Zhang_William/2016.02.29 spe-9 Run 12B',    #8
            r'/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13A',    #9
            r'/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13B',    #10
            r'/mnt/iscopearray/Zhang_William/2016.03.04 spe-9 Run 13C',    #11
            r'/mnt/iscopearray/Zhang_William/2016.03.14 spe-9 Run 14',     #12
            r'/mnt/iscopearray/Zhang_William/2016.03.25 spe-9 Run 15A',    #13
            r'/mnt/iscopearray/Zhang_William/2016.03.25 spe-9 Run 15B',    #14
            r'/mnt/iscopearray/Zhang_William/2016.03.31 spe-9 Run 16'],    #15
        'extra_directories':[
            {'W': r'/mnt/scopearray/ZhangWillie/2016.02.16 spe-9 Run 9'},      #0
            None,                                                                                     #1
            None,                                                                                     #2
            None,                                                                                     #3
            None,                                                                                     #4
            None,                                                                                     #5
            None,                                                                                     #6
            {'W': r'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12A'},    #7
            {'W': r'/mnt/scopearray/ZhangWillie/2016.02.29 spe-9 Run 12B'},    #8
            None,                                                                                     #9
            {'W': r'/mnt/scopearray/ZhangWillie/2016.03.04 spe-9 Run 13B'},    #10
            {'W': r'/mnt/scopearray/ZhangWillie/2016.03.04 spe-9 Run 13C'},    #11
            None,                                                                                     #12
            None,                                                                                     #13
            None,                                                                                     #14
            None],                                                                                     #15
        'experiment_directories':[
            'W',                                                                                      #0
            None,                                                                                     #1
            None,                                                                                     #2
            None,                                                                                     #3
            None,                                                                                     #4
            None,                                                                                     #5
            None,                                                                                     #6
            'W',                                                                                      #7
            'W',                                                                                      #8
            None,                                                                                     #9
            'W',                                                                                      #10
            'W',                                                                                      #11
            None,                                                                                     #12
            None,                                                                                     #13
            None,                                                                                     #14   
            None],
        'annotation_directories':[
            r'/mnt/purplearray/Sinha_Drew/20160216_spe9Acquisition',            #0
            None,                                                                                     #1
            None,                                                                                     #2
            None,                                                                                     #3
            None,                                                                                     #4
            None,                                                                                     #5
            None,                                                                                     #6
            r'/mnt/purplearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA',    #7
            r'/mnt/purplearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB',    #8
            None,                                                                                     #9
            r'/mnt/purplearray/Sinha_Drew/20160304_spe9Acquisition_DevVarB',    #10
            r'/mnt/purplearray/Sinha_Drew/20160304_spe9Acquisition_DevVarC',    #11
            None,                                                                                     #12
            None,                                                                                     #13
            None,                                                                                     #14
            None]
    },                                                                                     #15
    'age-1': {
        'data_directories':[
            r'/mnt/iscopearray/Zhang_William/2016.05.02 spe-9 age-1 osxi346 Run 17A', #16
            r'/mnt/iscopearray/Zhang_William/2016.05.02 spe-9 age-1 osxi346 Run 17B', #17
            r'/mnt/iscopearray/Zhang_William/2016.05.02 spe-9 age-1 osxi346 Run 17C', #18
            r'/mnt/iscopearray/Zhang_William/2016.05.02 spe-9 age-1 osxi346 Run 17D', #19
            r'/mnt/iscopearray/Zhang_William/2016.05.02 spe-9 age-1 osxi346 Run 17E', #20
            r'/mnt/scopearray/ZhangWillie/2016.05.12 spe-9 age-1 Run 18B',   #21
            r'/mnt/scopearray/ZhangWillie/2016.05.24 spe-9 age-1 Run 24A',  #22 7 worms alive as of 20160624
            r'/mnt/scopearray/ZhangWillie/2016.07.01 spe-9 age-1 Run 21',  #23
            r'/mnt/iscopearray/Sinha_Drew/2016.07.01 spe-9 age-1 Run 22',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
    'ife-2': {
        'data_directories':[
            r'/mnt/purplearray/Sinha_Drew/20160921 ife-2 spe-9 Run 23',
            r'/mnt/purplearray/Sinha_Drew/20160923 ife-2 spe-9 Run 24',
            r'/mnt/iscopearray/Sinha_Drew/20160927 ife-2 spe-9 Run 25A',
            r'/mnt/iscopearray/Sinha_Drew/20160927 ife-2 spe-9 Run 25B',
            r'/mnt/iscopearray/Sinha_Drew/20161206 ife-2 spe-9 Run 28',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
    'clk-1': {
        'data_directories':[
            r'/mnt/iscopearray/Sinha_Drew/20161025 clk-1 spe-9 Run 26A',
            r'/mnt/iscopearray/Sinha_Drew/20161025 clk-1 spe-9 Run 26B',
            r'/mnt/iscopearray/Sinha_Drew/20161027 clk-1 spe-9 Run 27B',
            r'/mnt/iscopearray/Sinha_Drew/20161209 clk-1 spe-9 Run 29',
            r'/mnt/iscopearray/Sinha_Drew/20161219 clk-1 spe-9 Run 30',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    }
}
for strain in data_list:
    for meta_dirs in ['extra_directories','experiment_directories','annotation_directories']:
        if data_list[strain][meta_dirs] is None:
            data_list[strain][meta_dirs] = [None for ii in range(len(data_list[strain]['data_directories']))]

def make_SVR(strains,svm_save_dir):
    save_directory = ''
    print('Work directory exists:'+str(os.path.isdir(working_directory)))
    print('Human directory exists:'+str(os.path.isdir(human_directory)))
    
    if strains == 'combined':
        compiled_list = {my_key:[] for my_key in ['data_directories','extra_directories','experiment_directories','annotation_directories']}
        for my_strain in data_list:
            for my_key in ['data_directories','extra_directories','experiment_directories','annotation_directories']:
                compiled_list[my_key].extend(data_list[my_strain][my_key])
        directory_bolus = folderStuff.DirectoryBolus(working_directory, human_directory, 
                *[compiled_list[dirs] for dirs in ['data_directories','extra_directories','experiment_directories','annotation_directories']], 
                ready = len(compiled_list['data_directories']))
        if svm_save_dir is not '':
            svm_dir_out = svm_save_dir+os.path.sep+'combined_health_SVR'+os.path.sep
            if not os.path.isdir(svm_dir_out): 
                print('(make_SVR) Making directory at: '+svm_dir_out)
                os.mkdir(svm_dir_out)
        adult_df = characterizeTrajectories.CompleteWormDF(directory_bolus, save_directory,
            {'adult_only': True, 'svm_dir_out':svm_dir_out})
    else:
        for strain in strains:
            directory_bolus = folderStuff.DirectoryBolus(working_directory, human_directory, 
                *[data_list[strain][dirs] for dirs in ['data_directories','extra_directories','experiment_directories','annotation_directories']], 
                ready = len(data_list[strain]['data_directories']))    
            
            if svm_save_dir is not '':
                svm_dir_out = svm_save_dir+os.path.sep+strain+'_health_SVR'+os.path.sep
                if not os.path.isdir(svm_dir_out): 
                    print('(make_SVR) Making directory at: '+svm_dir_out)
                    os.mkdir(svm_dir_out)
            adult_df = characterizeTrajectories.CompleteWormDF(directory_bolus, save_directory,
                {'adult_only': True, 'svm_dir_out':svm_dir_out})

def make_df(strains,df_savedir, svm_directory=''):
    if strains is 'combined':
        compiled_list = {my_key:[] for my_key in ['data_directories','extra_directories','experiment_directories','annotation_directories']}
        for my_strain in data_list:
            for my_key in ['data_directories','extra_directories','experiment_directories','annotation_directories']:
                compiled_list[my_key].extend(data_list[my_strain][my_key])
        directory_bolus = folderStuff.DirectoryBolus(working_directory, human_directory, 
                *[compiled_list[dirs] for dirs in ['data_directories','extra_directories','experiment_directories','annotation_directories']], 
                ready = len(compiled_list['data_directories']))
        if df_savedir is not '':
            full_savedir = df_savedir+os.path.sep+'combined_health'+os.path.sep
            if not os.path.isdir(full_savedir): 
                print('(make_df) Making directory at:'+full_savedir)
                os.mkdir(full_savedir)
        adult_df = characterizeTrajectories.CompleteWormDF(directory_bolus, full_savedir,
            {'adult_only': True, 'svm_directory':svm_directory})
        #~ life_df = characterizeTrajectories.CompleteWormDF(directory_bolus, full_savedir,
            #~ {'adult_only': False, 'svm_directory':svm_directory})
        with open(full_savedir+'df_combined.pickle','wb') as my_file:
            pickle.dump({'adult_df':adult_df,'svm_directory':svm_directory,'data_list':data_list},my_file)
    else:
        for strain in strains:
            directory_bolus = folderStuff.DirectoryBolus(working_directory, human_directory, 
                *[data_list[strain][dirs] for dirs in ['data_directories','extra_directories','experiment_directories','annotation_directories']], 
                ready = len(data_list[strain]['data_directories']))    
            
            if df_savedir is not '':
                full_savedir = df_savedir+os.path.sep+strain+'_health'+os.path.sep
                if not os.path.isdir(full_savedir): 
                    print('(make_df) Making directory at: '+full_savedir)
                    os.mkdir(full_savedir)
            adult_df = characterizeTrajectories.CompleteWormDF(directory_bolus, full_savedir,
                {'adult_only': True, 'svm_directory':svm_directory})
            #~ life_df = characterizeTrajectories.CompleteWormDF(directory_bolus, full_savedir,
                #~ {'adult_only': False, 'svm_directory':svm_directory})
            with open(full_savedir+'df_'+strain+'.pickle','wb') as my_file:
                pickle.dump({'adult_df':adult_df,'svm_directory':svm_directory,'data_list':data_list},my_file)
