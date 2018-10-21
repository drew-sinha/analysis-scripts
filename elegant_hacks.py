
import numpy
import matplotlib.pyplot as plt; plt.ion(); plt.show()
import matplotlib.cm

from elegant import worm_data
from zplib.image import colorize

import plotting_tools

#==================================
# Plotting timecourses
#==================================

def _timecourse_plot_data(worms, feature, min_age=-numpy.inf, max_age=numpy.inf, age_feature='age', color_by='lifespan',
    color_map=None):

    if color_map is None:
        color_map = lambda array: colorize.color_map(array, uint8=False)

    time_ranges = worms.get_time_range(feature, min_age, max_age, age_feature)
    if color_by is None:
        color_vals = [0]*len(times_range)
    else:
        color_vals = colorize.scale(worms.get_feature(color_by), output_max=1)
    colors = color_map(color_vals)
    out = []
    for time_range, color in zip(time_ranges, colors):
        x, y = time_range.T
        out.append((x, y, color))
    return out

def plot_timecourse(worms, feature, min_age=-numpy.inf, max_age=numpy.inf,
    age_feature='age', time_units='hours', color_by='lifespan',
    color_map=None, **plot_kws):
    """Plot values of a given feature for each worm, colored by a given
    worm feature (defaults to lifespan).

    Parameters:
        feature: the name of a timecourse feature available in worm.td to
            retrieve (such as "gfp_95th", say), or a function that will be
            called as feature(worm) that will calculate a timecourse feature
            (see examples in Worm.get_time_range).
        min_age, max_age: the beginning and end of the window in which to
            retrieve features. If not specified, features from the very
            beginning and/or to the very end of the timecourse will be
            retrieved.
        age_feature: the name of the feature to use to determine the age
            window and the plot axes. Defaults to 'age', but other
            monotonic features, such as a time-left-alive 'ghost_age' could
            also be used. If this is a function, it will be called as
            age_feature(worm) to generate the ages to examine (see
            examples in Worm.get_time_range).
        time_units: one of "days", "hours", "minutes", or "seconds",
            representing the units in which the time scale should be plotted.
        color_by: worm feature to use for color scale of each timecourse.
    """
    if time_units not in worm_data.TIME_UNITS:
        raise ValueError(f"'time_units' must be one of: {list(TIME_UNITS)}")

    time_scale = worm_data.TIME_UNITS[time_units]
    fig_h, ax_h = plt.subplots()
    ax_h.set_xlabel(f'{age_feature} ({time_units})')
    ax_h.set_ylabel(f'{feature}')
    for x, y, c in _timecourse_plot_data(worms, feature, min_age, max_age, age_feature, color_by,color_map):
        ax_h.plot(x/time_scale, y, color=c, **plot_kws)
    return fig_h, ax_h

'''
    # Code for monochromatic colormap....

    base_color = [] # FILL in.
    cmap = cm.LinearlySegmentedColormap('Monochrome_lab',
        plotting_tools.build_gradient_palette(base_color,256))  # Assume palette is fine enough to linearly segment
    plot_timecourse(worms, feature, min_age=0, age_feature='adult_age', time_units='days', color_map=cmap)
'''
