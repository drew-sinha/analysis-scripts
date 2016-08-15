
def bootstrap_long(data,stat,num_samples, num_reps):
    resamp_data = np.array([resample(my_time, num_samples, num_reps) for my_time in data.T]) # resamp_data.shape = [num_cols, num_rep, num_samps]
    return [stat(rep) for rep in np.swapaxes(np.swapaxes(resamp_data,0,1),1,2]  # resamp_data.shape = [num_reps, num_samps, num_cols] after swapping

# Resample from 1D array
def resample(data, num_samples, num_reps):
    return data[np.randint(0,np.len(data),[num_reps, num_samples]).flatten()].reshape([num_reps, num_samples])

# def bootstrap(data, stat, num_samples, num_reps):
    # '''
        # data - Data to be resampled in a flattened array
        # stat - Function handle to be applied to data; signature: stat:=stat(resampled) where resampled has a dimension of 2
        # num_samples - number of samples for a given repetition/resampling
        # num_reps - number of reptitions/resamplings (an array giving the axes of the passed back data)
        
        # Note: stat should be consistent with number of dimensions of data; namely, stat should act with a num_dim(data)+1 matrix
            # e.g. if data is 1D, stat should act over the columns of a 2D matrix with shape (num_rep, num_samp)
    # '''
    # if(np.len(data.shape) == 1):    # Sampling from 1D array
        # resamp_data = resample(data,num_samples,num_reps)
        # return np.apply_along_axis(stat, 1, resamp_data)    # 1 - apply over columns of data
    # elif(np.len(data.shape) ==2):   # 2D
        # resamp_data = np.apply_along_axis(resample, 0, data)    # resamp_data.shape = [num_cols, num_rep, num_samps]
        # return np.array([stat(resampl_data[my_slice,...]) for my_slice in resamp_data])
        
# def bootstrap_long(data, stat, num_samples, num_reps)
    # '''
        # data - Longitudinal data to bootstrap in a 2D matrix with shape (num_indiv, num_times)
        # stat - Function handle to statistic operating over a vector of observed values (e.g. one variable tracked for an individual over time)
        # num_samples
        # num_reps
    # '''
    # resamp_data = np.array([resample(my_time, num_samples, num_reps) for my_time in data.T]) # resamp_data.shape = [num_cols, num_rep, num_samps]
    # return [[stat(samp) for samp in rep] for rep in np.swapaxes(np.swapaxes(resamp_data,0,1),1,2)]
    

    
    
'''
    Stats
    ==============
    Cumulative absolute deviation from the group mean/median
    Cumulative absolute z-scored/percentile
    # of sign changes
'''