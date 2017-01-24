import annotation_file
import matplotlib.pyplot as plt
import numpy as np

def get_lsdata(expt_dirs,md_list=None):
    timestamped_data = annotation_file.compile_expt_timestamped_data(expt_dirs,md_list)
    lifespan = (timestamped_data['Death']-timestamped_data['Hatch'])/(3600*24) # Days
    return (lifespan, timestamped_data['Worm_FullName'])
    

expt_dirs = [
    '/mnt/iscopearray/zhangw/2016.02.16 spe-9 Run 9/',      #0
    '/mnt/iscopearray/zhangw/2016.02.20 spe-9 Run 10A/',    #1
    '/mnt/iscopearray/zhangw/2016.02.20 spe-9 Run 10B/',    #2
    '/mnt/iscopearray/zhangw/2016.02.26 spe-9 Run 11A/',    #3
    '/mnt/iscopearray/zhangw/2016.02.26 spe-9 Run 11B/',    #4
    '/mnt/iscopearray/zhangw/2016.02.26 spe-9 Run 11C/',    #5
    '/mnt/iscopearray/zhangw/2016.02.26 spe-9 Run 11D/',    #6
    '/mnt/iscopearray/zhangw/2016.02.29 spe-9 Run 12A/',    #7
    '/mnt/iscopearray/zhangw/2016.02.29 spe-9 Run 12B/',    #8
    '/mnt/iscopearray/zhangw/2016.03.04 spe-9 Run 13A/',    #9
    '/mnt/iscopearray/zhangw/2016.03.04 spe-9 Run 13B/',    #10
    '/mnt/iscopearray/zhangw/2016.03.04 spe-9 Run 13C/',    #11
    '/mnt/iscopearray/zhangw/2016.03.14 spe-9 Run 14/',     #12
    '/mnt/iscopearray/zhangw/2016.03.25 spe-9 Run 15A/',    #13
    '/mnt/iscopearray/zhangw/2016.03.25 spe-9 Run 15B/',    #14
    '/mnt/iscopearray/zhangw/2016.03.31 spe-9 Run 16/']     #15
md_list = [
    {'W':'/mnt/iscopearray/zhangw/2016.02.16 spe-9 Run 9/experiment_metadata.json','D':'/mnt/scopearray/Sinha_Drew/20160216_spe9Acquisition/experiment_metadata.json'},   #0
    {'':''},    #1
    {'':''},    #2
    {'':''},    #3
    {'':''},    #4
    {'':''},    #5
    {'':''},    #6
    {'W':'/mnt/iscopearray/zhangw/2016.02.29 spe-9 Run 12A/experiment_metadata.json', 'D':'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarA/experiment_metadata.json'},    #7
    {'W':'/mnt/iscopearray/zhangw/2016.02.29 spe-9 Run 12B/experiment_metadata.json','D':'/mnt/scopearray/Sinha_Drew/20160229_spe9Acquisition_DevVarB/experiment_metadata.json'},    #8
    {'':''},    #9
    {'W':'/mnt/iscopearray/zhangw/2016.03.04 spe-9 Run 13B/experiment_metadata.json', 'D':'/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarB/experiment_metadata.json'},    #10
    {'W':'/mnt/iscopearray/zhangw/2016.03.04 spe-9 Run 13B/experiment_metadata.json', 'D':'/mnt/scopearray/Sinha_Drew/20160304_spe9Acquisition_DevVarC/experiment_metadata.json'},    #11
    {'':''},    #12
    {'':''},    #13
    {'':''},    #14
    {'':''},    #15
]

plt.close('all')
plt.ion()
plt.show()
    
(spe_ls, spe_worms) = get_lsdata(expt_dirs,md_list)

# Derive survival curve implicitly by treating deaths event-wise
sorted_ls = np.concatenate(([0],sorted(spe_ls)))
prop_alive = np.concatenate((1 -(np.arange(start=0,stop=np.size(spe_ls)))/np.size(spe_ls),[0]))

# Mortality rate
bins = np.arange(start=0, stop=int(spe_ls.max()//1 + 1),step=0.5)
(binned_deaths,trash) = np.histogram(spe_ls, bins)
remaining_pop = np.concatenate(([np.sum(binned_deaths)],(np.sum(binned_deaths)-np.cumsum(binned_deaths))[:-1]))
raw_rate=binned_deaths/remaining_pop

# Plot survival curve, lifespans, mortality rate
fig_h, ax_h = plt.subplots(3,1,sharex=True)
ax_h[0].plot(sorted_ls,prop_alive)
ax_h[1].hist(spe_ls, np.arange(0,(np.floor_divide(spe_ls.max(),2)+1)*2))
ax_h[2].plot(bins[:-1],np.log10(raw_rate))

import scipy.stats
no_mortality = raw_rate==0
lin_fit = scipy.stats.linregress(bins[:-1][~no_mortality],np.log10(raw_rate)[~no_mortality])
