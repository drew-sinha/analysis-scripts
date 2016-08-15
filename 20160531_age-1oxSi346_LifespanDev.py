import plotting_tools

if __name__=="__main__":
    plotting_tools.quick_plot_dev(
    ['/mnt/scopearray/Sinha_Drew/20160502_spe9age1oxsi346_10ms/20160502_age1-spe9-oxSi346_Development - 10msExposure_WormNotes.tsv',
    '/mnt/scopearray/Sinha_Drew/20160502_spe9age1oxsi346_20ms/20160502_age1-spe9-oxSi346_Development - 20msExposure_WormNotes.tsv'],
    ['/mnt/scopearray/Sinha_Drew/20160502_spe9age1oxsi346_10ms/experiment_metadata.json','/mnt/scopearray/Sinha_Drew/20160502_spe9age1oxsi346_20ms/experiment_metadata.json'],
    bad_worm_kws = ['FERTILITY', 'Nh', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS','Nw','DAUER','BURST']
    )

    plotting_tools.quick_plot_lifespan(
    ['/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17A/2016.05.02_age-1-spe-9-oxSi346_Lifespan - RunA_WormNotes.tsv',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17B/2016.05.02_age-1-spe-9-oxSi346_Lifespan - RunB_WormNotes.tsv',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17C/2016.05.02_age-1-spe-9-oxSi346_Lifespan - RunC_WormNotes.tsv',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17D/2016.05.02_age-1-spe-9-oxSi346_Lifespan - RunD_WormNotes.tsv',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17E/2016.05.02_age-1-spe-9-oxSi346_Lifespan - RunE_WormNotes.tsv'],
    ['/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17A/experiment_metadata.json',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17B/experiment_metadata.json',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17C/experiment_metadata.json',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17D/experiment_metadata.json',
    '/mnt/scopearray/ZhangWillie/2016.05.02 spe-9 age-1 osxi346 Run 17E/experiment_metadata.json'],
    bad_worm_kws = ['FERTILITY', 'Nh', 'DOUBLE WORM', 'OUT OF FOCUS', 'NO EGGS','NO WORM', 'RAN LOW ON FOOD', 'NEVER LAID EGGS','Nw','DAUER','BURST']
    )
