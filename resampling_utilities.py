import numpy as np
import scipy.stats


'''
Permutation test stuff
Taken fron 20170511_daf-16 notebook (and that from 20170217_..._FinalFollowUp)
'''

def permute_data(data,num_samples):
    return np.array([
        np.random.permutation(data) for i in range(num_samples)
    ])

def resample_data(data,num_samples):
    data = np.array(data)
    #~ return np.array(
            #~ data[[np.random.randint(len(data),(len(data),))]] for i in num_samples)
    return data[[np.random.randint(len(data),(num_samples,len(data)))]]

def permutetest_meandiff(samp1,samp2,num_samples=10000): # samp2-samp1, default to two-sided test
    resampled_data = permute_data(np.append(samp1,samp2),num_samples=num_samples)
    resampled_diffs = np.mean(resampled_data[:,:len(samp2)],axis=1)-np.mean(resampled_data[:,len(samp2):],axis=1)
    obs_diff = np.mean(samp2)-np.mean(samp1)
    return np.count_nonzero(
        np.abs(resampled_diffs)>np.abs(obs_diff))/num_samples

def permutetest_mediandiff(samp1,samp2,num_samples=10000): # samp2-samp1, default to two-sided test
    resampled_data = permute_data(np.append(samp1,samp2),num_samples=num_samples)
    resampled_diffs = np.median(resampled_data[:,:len(samp2)],axis=1)-np.median(resampled_data[:,len(samp2):],axis=1)
    obs_diff = np.median(samp2)-np.median(samp1)
    return np.count_nonzero(
        np.abs(resampled_diffs)>np.abs(obs_diff))/num_samples

def permutetest_levene(samp1,samp2,num_samples=10000):
    resampled_data = permute_data(np.append(samp1,samp2),num_samples=num_samples)
    resampled_L = [scipy.stats.levene(sample[:len(samp1)],sample[len(samp1):])[0] for sample in resampled_data]
    obs_L = scipy.stats.levene(samp1,samp2)[0]
    return np.count_nonzero(resampled_L>obs_L)/num_samples

def permutetest_tstat(samp1,samp2,num_samples=10000):
    resampled_data = permute_data(np.append(samp1,samp2),num_samples=num_samples)
    resampled_t = scipy.stats.ttest_ind(resampled_data[:,:len(samp1)],resampled_data[:,len(samp1):],axis=1,equal_var=False)[0]
    obs_t = scipy.stats.ttest_ind(samp1,samp2,equal_var=False)[0]
    return np.count_nonzero(np.abs(resampled_t)>np.abs(obs_t))/num_samples
    
def bootstrap_CI(samp,num_samples=10000,alpha=0.95):
    resampled_data = resample_data(samp)
    resampled_data = resampled_data - resampled_data.mean(axis=0)
    
    return samp.mean() - np.percentile(resampled_data,[alpha/2, (1-alpha)/2])
