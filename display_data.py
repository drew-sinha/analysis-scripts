import numpy
import matplotlib.pyplot as plt

from elegant import worm_data

def display_worm_subset(worms, feature, specified_subset=None, **plot_timecourse_kws):
    '''
        specified_subset - either number of worms, names, or list of indices
    '''

    try:
        if type(specified_subset[0]) is str:
            subset = [worm for worm in worms if worm.name in specified_subset]
        else:
            subset = numpy.array(worms)[specified_subset]
    except TypeError: # Not subscriptable, so int
        subset = numpy.random.permutation(worms)[:specified_subset]

    subset = worm_data.Worms(subset)
    fig_h, ax_h = plt.subplots()
    subset.plot_timecourse(feature, **plot_timecourse_kws)

    return fig_h, ax_h
