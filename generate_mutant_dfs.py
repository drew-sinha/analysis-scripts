from analyzeHealth import characterizeTrajectoriesMinimal
import pickle
import pathlib

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
            r'/mnt/iscopearray/Sinha_Drew/2016.05.12 spe-9 age-1 Run 18B',   #21
            r'/mnt/iscopearray/Sinha_Drew/2016.05.24 spe-9 age-1 Run 24A',  #22 7 worms alive as of 20160624
            r'/mnt/iscopearray/Sinha_Drew/2016.07.01 spe-9 age-1 Run 21',  #23
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
    },
    'daf-16': {
        'data_directories':[
            r'/mnt/purplearray/Sinha_Drew/20170206 daf-16 spe-9 Run 32A',
            r'/mnt/purplearray/Sinha_Drew/20170407 daf-16 spe-9 Run 33A',
            r'/mnt/purplearray/Sinha_Drew/20170407 daf-16 spe-9 Run 33B',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
    'pqm-1': {
        'data_directories':[
            r'/mnt/purplearray/Sinha_Drew/20170706 pqm-1 spe-9 Run 1A',
            r'/mnt/purplearray/Sinha_Drew/20170706 pqm-1 spe-9 Run 1B',
            r'/mnt/purplearray/Sinha_Drew/20170717 pqm-1 spe-9 Run 2',
            r'/mnt/purplearray/Sinha_Drew/20170721 pqm-1 spe-9 Run 3',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
    'mev-1': {
        'data_directories':[
            r'/mnt/purplearray/Sinha_Drew/20170825 mev-1 spe-9 Run 4',
            r'/mnt/purplearray/Sinha_Drew/20170831 mev-1 spe-9 Run 6A',
            r'/mnt/purplearray/Sinha_Drew/20170908 mev-1 spe-9 Run 7A',
            r'/mnt/purplearray/Sinha_Drew/20170908 mev-1 spe-9 Run 7B',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
    'rsks-1': {
        'data_directories':[
            r'/mnt/scopearray/Sinha_Drew/20171030 rsks-1 spe-9 Run 1A',
            r'/mnt/scopearray/Sinha_Drew/20171030 rsks-1 spe-9 Run 1B',
            r'/mnt/scopearray/Sinha_Drew/20171114 rsks-1 spe-9 Run 2A',
            r'/mnt/scopearray/Sinha_Drew/20171114 rsks-1 spe-9 Run 2B',
            r'/mnt/scopearray/Sinha_Drew/20171213 rsks-1 spe-9 Run 3',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
    'age-1;daf-16': {
        'data_directories':[
            r'/mnt/scopearray/Sinha_Drew/20180108 age-1 daf-16 spe-9 Run 3A',
            r'/mnt/scopearray/Sinha_Drew/20180108 age-1 daf-16 spe-9 Run 3B',
            r'/mnt/scopearray/Sinha_Drew/20180130 age-1 daf-16 spe-9 Run 4A',
            r'/mnt/scopearray/Sinha_Drew/20180130 age-1 daf-16 spe-9 Run 4B',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },    
    'spe-9new': {
        'data_directories':[
            r'/mnt/scopearray/Sinha_Drew/20180108 spe-9 Ctrl',
            r'/mnt/scopearray/Sinha_Drew/20180130 spe-9 Ctrl',
            r'/mnt/scopearray/Sinha_Drew/20180220 spe-9 Ctrl',
            r'/mnt/scopearray/Sinha_Drew/20180309 spe-9 Ctrl',
        ],
        'extra_directories':None,
        'experiment_directories':None,
        'annotation_directories':None
    },
}
for strain in data_list:
    for meta_dirs in ['extra_directories','experiment_directories','annotation_directories']:
        if data_list[strain][meta_dirs] is None:
            data_list[strain][meta_dirs] = [None for ii in range(len(data_list[strain]['data_directories']))]

def make_SVR(strains,svm_save_dir,**df_extra_args):
    #~ save_directory = ''
    #~ print('Work directory exists:'+working_directory.is_dir())
    #~ print('Utility directory exists:'+working_directory.is_dir())
    #~ df_extra_args.setdefault('adult_only',True)
    #~ df_extra_args.setdefault('add_health',True)
    #~ df_extra_args['svm_save_dir'] = svm_save_dir
    
    #~ if strains == 'combined':
        #~ compiled_list = {my_key:[] for my_key in ['data_directories','extra_directories','experiment_directories','annotation_directories']}
        #~ for my_strain in data_list:
            #~ for my_key in ['data_directories','extra_directories','experiment_directories','annotation_directories']:
                #~ compiled_list[my_key].extend(data_list[my_strain][my_key])
        #~ directory_bolus = folderStuff.DirectoryBolus(working_directory, human_directory, 
                #~ *[compiled_list[dirs] for dirs in ['data_directories','extra_directories','experiment_directories','annotation_directories']], 
                #~ ready = len(compiled_list['data_directories']))
        #~ if svm_save_dir is not '':
            #~ svm_dir_out = svm_save_dir+os.path.sep+'combined_health_SVR'+os.path.sep
            #~ if not os.path.isdir(svm_dir_out): 
                #~ print('(make_SVR) Making directory at: '+svm_dir_out)
                #~ os.mkdir(svm_dir_out)
        #~ adult_df = characterizeTrajectories.CompleteWormDF(directory_bolus, save_directory,
            #~ df_extra_args)
    #~ else:
        #~ for strain in strains:
            #~ directory_bolus = folderStuff.DirectoryBolus(working_directory, human_directory, 
                #~ *[data_list[strain][dirs] for dirs in ['data_directories','extra_directories','experiment_directories','annotation_directories']], 
                #~ ready = len(data_list[strain]['data_directories']))    
            
            #~ if svm_save_dir is not '':
                #~ svm_dir_out = svm_save_dir+os.path.sep+strain+'_health_SVR'+os.path.sep
                #~ if not os.path.isdir(svm_dir_out): 
                    #~ print('(make_SVR) Making directory at: '+svm_dir_out)
                    #~ os.mkdir(svm_dir_out)
            #~ adult_df = characterizeTrajectories.CompleteWormDF(directory_bolus, save_directory,
                #~ df_extra_args)
                
    '''
        TODO - get rid of this or reimplement such that it utilizes the refactored fit_regressor func of the DF
    '''
    raise NotImplementedError

def make_df(strains,df_savedir='',custom_savename='',**df_extra_args):
    '''
        Important df arguments
            regressor_fp - (str/pathlib.Path) filepath to health score regressor
    '''
    
    df_extra_args.setdefault('adult_only',True)
    df_extra_args.setdefault('bad_worm_kws',['PVL','BURST'])
    df_savedir = pathlib.Path(df_savedir)
    
    dfs = []
    for strain in strains:  
        adult_df = characterizeTrajectoriesMinimal.BasicWormDF(
            data_list[strain]['data_directories'],
            **df_extra_args)
        
        if str(df_savedir) is not '':
            data_to_save = {'adult_df':adult_df,'data_list':data_list[strain]}
            data_to_save.update(df_extra_args)
            with (df_savedir/('df_'+strain+custom_savename+'.pickle')).open('wb') as my_file:
                pickle.dump(data_to_save,my_file)
        dfs.append(adult_df)
    return dfs
