import scipy.stats as stats
import sklearn.svm as svm
import sklearn.model_selection as model_selection

import concurrent.futures
import multiprocessing
import numpy as np

import pathlib

import time

import argparse
import pickle
import string
import platform

import analyzeHealth.selectData as selectData


def randomsearch_fitparams_epsSVR(X,y,n_iter=20,n_jobs=1):
    regressor = svm.SVR(kernel='rbf')
    param_dist = {'C': stats.expon(scale=50), 'gamma': stats.expon(scale=1), 'epsilon':stats.expon(scale=10)}
    cv = model_selection.ShuffleSplit(n_splits=5, test_size=0.7)
    random_search = model_selection.RandomizedSearchCV(regressor, param_distributions=param_dist, n_iter=n_iter, verbose=True, cv=cv, n_jobs=n_jobs)
    random_search.fit(X, y)
    return random_search

def gridsearch_fitparams_epsSVR(X,y,n_iter=20,n_jobs=1):
    regressor = svm.SVR(kernel='rbf')
    param_dist = {'C': stats.expon(scale=50), 'gamma': stats.expon(scale=1), 'epsilon': stats.expon(scale=10)}
    cv = model_selection.ShuffleSplit(n_splits=5, test_size=0.7)
    random_search = model_selection.RandomizedSearchCV(regressor, param_distributions=param_dist, n_iter=n_iter, verbose=True, cv=cv, n_jobs=n_jobs)
    #~ print('Finished searching; proceeding to fit model.')
    #~ random_search.fit(X, y)
    return random_search
    
func_lookup_table = {
    'randomsearch_fitparams_epsSVR': randomsearch_fitparams_epsSVR,
    'gridsearch_fitparams_epsSVR': gridsearch_fitparams_epsSVR
}

def paramsearch_df_parallel(df,param_search_func, num_workers=4):
    
    '''
        Still in beta
    '''
    
    biomarker_data = df.mloc(measures=raw_health_vars) #[animals, measurements, time]
    biomarker_data_flat = biomarker_data.swapaxes(1,2).reshape(-1,len(raw_health_vars)) #[animals x time,measurements]

    remaining_life = df.mloc(measures=['ghost_age'])[:,0,:] #[animals,time]
    remaining_life_flat = remaining_life.flatten() #[animals xtime]

    nan_mask = np.isnan(biomarker_data_flat).any(axis=1) | np.isnan(remaining_life_flat)
    with concurrent.futures.ProcessPoolExecutor(max_workers = num_workers) as executor:
        search_results = [
            executor.submit(param_search_func, 
                biomarker_data_flat[~nan_mask], 
                remaining_life_flat[~nan_mask])]
    concurrent_futures.wait(param_search_results)
    
    compiled_search_results = []
    for results in search_results:
        compiled_search_results.extend(results)
    return compiled_search_results
    
def paramsearch_df(df,param_search_func,selected_worms=None, **func_kwargs):
    if selected_worms is not None:
        print('Using custom-selected worms; length = {}'.format(len(selected_worms)))
    
    raw_health_vars = ['bulk_movement','stimulated_rate_a', 'stimulated_rate_b', 'unstimulated_rate', 'intensity_80', 'life_texture','adjusted_size','adjusted_size_rate', 'cumulative_eggs','cumulative_eggs_rate']
    biomarker_data = df.mloc(worms=selected_worms,measures=raw_health_vars) #[animals, measurements, time]
    biomarker_data_flat = biomarker_data.swapaxes(1,2).reshape(-1,len(raw_health_vars)) #[animals x time,measurements]

    remaining_life = df.mloc(worms=selected_worms,measures=['ghost_age'])[:,0,:] #[animals,time]
    remaining_life_flat = remaining_life.flatten() #[animals xtime]

    nan_mask = np.isnan(biomarker_data_flat).any(axis=1) | np.isnan(remaining_life_flat)
    search_results = param_search_func(biomarker_data_flat[~nan_mask], 
                        remaining_life_flat[~nan_mask],
                        **func_kwargs)
    
    return search_results

def write_cluster_jobfile(save_directory, **template_kws):
    jobscript_template = string.Template(
        '''
        #PBS -N $job_name
        #PBS -t 1-$num_runs
        #PBS -l nodes=1:ppn=1
        #PBS -l walltime=5:00:00
        #PBS -M drew.sinha@wustl.edu
        #PBS -d /scratch/sinhad/
        #PBS -e /scratch/sinhad/job_log_files/
        #PBS -o /scratch/sinhad/job_log_files/
        #PBS -q old

        export ANACONDA_PATH="/act/Anaconda3-2.3.0/bin"
        export PATH=$$ANACONDA_PATH:$$PATH

        source activate zplab_cluster
        cd $$PBS_O_WORKDIR
        python /home/sinhad/Scripts/df_param_search.py /scratch/sinhad/work_dir/param_search/dfs/WT_classifier_newcode_20171011/spe-9_health/df_spe-9.pickle $param_search_func $${PBS_ARRAYID} n_jobs=1, n_iter=$n_iter
        ''')
    param_search_func = template_kws['param_search_func']
    assert param_search_func in func_lookup_table.keys()
    
    template_kws.setdefault('n_iter',30)
    template_kws.setdefault('num_runs',5)
        
    code = jobscript_template.substitute(job_name='param_search_'+param_search_func,
        **template_kws)
    with (pathlib.Path(save_directory) / 'job_script_paramsearch.sh').open('w') as job_script:
        job_script.write(code)
        
def compile_paramsearch_results(result_directory,save_directory=None):
    '''
        Compile data from multiple separate parameter searches into a single serializable object (structure of sklearn.model_selection.RandomizedSearchCV) with additional attributes storing the top clasifier from each parameter search and its score (top_classifiers_ and top_scores_, resp.)

        Arguments:
            result_directory - directory with saved parameter search data
            save_directory - directory to save pickled results
    '''

    result_directory = pathlib.Path(result_directory)
    result_data = None
    for i,result_fn in enumerate(result_directory.iterdir()):
        if ((result_fn.parts[-1].split('.'))[0] == 'compiled_results') or not (result_fn.suffix[1:] == 'pickle'): 
            continue # suffix starts with '.'
        
        #~ print('Working on file {}/{}'.format(i,len(list(result_directory.iterdir()))),sep='',end='',flush=True)
        print('Working on file {}/{}'.format(i,len(list(result_directory.iterdir()))),sep='',end='\r',flush=True)
        with result_fn.open('rb') as result_file:
            load_data = pickle.load(result_file) # RandomizedSearchCV
            #print(load_data.cv_results_)
            if result_data is None: 
                result_data = load_data
                result_data.top_classifiers_ = [load_data.best_estimator_]
                result_data.top_scores_ = [load_data.best_score_]
            else:
                result_data.top_classifiers_.extend([load_data.best_estimator_])
                result_data.top_scores_.extend([load_data.best_score_])
                if result_data.best_score_ < load_data.best_score_:
                    result_data.best_score_ = load_data.best_score_
                    result_data.best_params_ = load_data.best_params_
                    result_data.best_estimator_ = load_data.best_estimator_
                    result_data.best_index_ = len(result_data.cv_results_['params']) + load_data.best_index_
                result_data.cv_results_ = {
                    attr:np.append(result_data.cv_results_[attr],load_data.cv_results_[attr])
                    for attr in result_data.cv_results_.keys()}
    result_data.cv_results_['rank_test_score'] = np.array([]) # Eliminate this since these ranks are based on the individual searches
    
    if save_directory is not None:
        with (save_directory / 'compiled_results.pickle').open('wb') as result_file:
            pickle.dump(result_data,result_file)
    return result_data


    
if __name__ == "__main__":
    ''' Call with signature
    python df_param_search.py DF_PATH [SAVE_DIR] PARAM_SEARCH_FUNC_STR ARRAYID [LS_PERCENTILE=LOW,HIGH] [ARG1=VAL1 [ARG2=VAL2...]]
    '''
    #~ multiprocessing.set_start_method('spawn',force=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('df_path',type=str)
    parser.add_argument('save_dir',type=str,nargs='?',default=None)
    parser.add_argument('param_search_func',type=str)
    parser.add_argument('job_id',type=int)
    parser.add_argument('kws',nargs='*')

    args = parser.parse_args()    

    kwargs = {}
    for kw_str in args.kws:
        kw_parsed = kw_str.split(sep='=')
        filtered_kw = kw_parsed[1].replace('-','')
        if kw_parsed[0].lower() == 'ls_percentile':
            kwargs[kw_parsed[0].lower()] = list(map(float,kw_parsed[1][1:-1].split(',')))
        else:
            if filtered_kw.isnumeric(): 
                kw_parsed[1] = int(kw_parsed[1])
            kwargs[kw_parsed[0]] = kw_parsed[1]    

    # Baby-sit loading for simultaneous access attempts on cluster
    READ_SUCCESS = False
    while not READ_SUCCESS:
        try:
            with pathlib.Path(args.df_path).open('rb') as df_file:
                my_df = pickle.load(df_file)['adult_df']
            READ_SUCCESS = True
        except PermissionError:
            pass

    if kwargs.get('ls_percentile',None) is None:
        selected_worms = None
    else:
        low_percentile, high_percentile = kwargs.pop('ls_percentile')
        lifespans = selectData.get_adultspans(my_df)
        low_cutoff, high_cutoff = np.percentile(lifespans,[low_percentile, high_percentile])
        selected_worms = np.array(my_df.worms)[
            (lifespans >= low_cutoff) & (lifespans <= high_cutoff)]
    
    timestamp = time.strftime('%Y-%m-%d_%H%M')
    
    search_results = paramsearch_df(my_df,func_lookup_table[args.param_search_func],selected_worms, **kwargs)
    
    if args.save_dir is None:
        if 'centos' in platform.platform(): # Environment @ CHPC
            save_dir = pathlib.Path('/scratch/sinhad/work_dir/param_search/results')
        else:
            save_dir = pathlib.Path('/home/drew/temp')
    else:
        save_dir = pathlib.Path(args.save_dir)

    with (save_dir / (timestamp + '_' + str(args.param_search_func) + '_results_' + str(args.job_id) +'.pickle')).open('wb') as result_file:
        pickle.dump(search_results, result_file)
    
