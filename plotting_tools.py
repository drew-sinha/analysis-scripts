import numpy as np
import matplotlib.pyplot as plt
import os
import zplib.scalar_stats.kde
import skimage.color

from utilities import utilities

#~ import annotation_file

# I Want Hue color map
#~ qual_colors = np.array([[0,0,0],
    #~ [211,66,126],
    #~ [104,116,201],
    #~ [127,111,47],
    #~ [172,88,197],
    #~ [215,144,72],
    #~ [178,176,68],
    #~ [194,99,109],
    #~ [190,117,177],
    #~ [93,174,70],
    #~ [86,165,116],
    #~ [81,173,208],
    #~ [203,81,54]])/255
qual_colors = np.array([[0,0,0],    # Switched the 2nd and 9th colors
    [93,174,70],
    [104,116,201],
    [127,111,47],
    [172,88,197],
    [215,144,72],
    [178,176,68],
    [194,99,109],
    [81,173,208],   # Moved this guy up from 2nd to last
    [190,117,177],
    [211,66,126],
    [86,165,116],
    [203,81,54]])/255

def build_gradient_palette(base_color, num_colors):
    '''
        Build a gradient color palette using a provided color as the base and modifying the intensity in LAB space

        base_color - 3-item list or np.array containing the base color (no alpha channel)
        num_colors - (int) Number of desired colors
    '''

    # Handle black (really want white here...)
    if (np.array(base_color)==np.array([0,0,0])).all():
        base_color = np.array([1,1,1])*.99 # can't do pure white for stability

    base_color_lab = skimage.color.rgb2lab(np.array([[base_color]]))

    return ([skimage.color.lab2rgb(
        np.array([[[
            100-(base_color_lab[0,0,0]*i/num_colors),
            base_color_lab[0,0,1],
            base_color_lab[0,0,2]]]]))[0,0,:]
        for i in range(num_colors+1)][1:])

#~ def quick_plot_dev(ann_fps, expt_mds,bad_worm_kws=[]):
    #~ my_ann_files = [annotation_file.AnnotationFile(ann_fp) for ann_fp in ann_fps]

    #~ timestamped_data = {}
    #~ [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    #~ for [expt_md_fp, ann_file] in zip(expt_mds, my_ann_files):
        #~ ann_file_data = ann_file.data_as_timestamps_simple(expt_md_fp)
        #~ for expt_key in timestamped_data.keys():
            #~ timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])
    #~ print(timestamped_data)

    #~ viable_worm = (timestamped_data['Hatch']!=-1) \
        #~ & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    #~ hatched_on_corral = timestamped_data['Hatch'] != 0

    #~ L1_durations = (timestamped_data['L1 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    #~ L2_durations = (timestamped_data['L2 ecdysis']-timestamped_data['L1 ecdysis'])[viable_worm]/3600
    #~ L3_durations = (timestamped_data['L3 ecdysis']-timestamped_data['L2 ecdysis'])[viable_worm]/3600
    #~ L4_durations = (timestamped_data['L4 ecdysis']-timestamped_data['L3 ecdysis'])[viable_worm]/3600
    #~ larval_durations = (timestamped_data['L4 ecdysis']-timestamped_data['Hatch'])[viable_worm & hatched_on_corral]/3600
    #~ print(L1_durations)

    #~ plt.show()
    #~ plt.ion()

    #~ plt.gcf().clf()
    #~ fig_h,ax_h = plt.subplots(5,1,sharex=True)

    #~ ax_h[0].hist(L1_durations)
    #~ ax_h[0].set_xlabel('Duration (hr)')
    #~ ax_h[0].set_ylabel('Frequency')
    #~ ax_h[0].set_title('L1 (n={}) - Mean: {}, Std:{}'.format(np.size(L1_durations),np.mean(L1_durations),np.std(L1_durations)))

    #~ ax_h[1].hist(L2_durations)
    #~ ax_h[1].set_title('L2 (n={}) - Mean: {}, Std:{}'.format(np.size(L2_durations),np.mean(L2_durations),np.std(L2_durations)))

    #~ ax_h[2].hist(L3_durations)
    #~ ax_h[2].set_title('L3 - Mean: {}, Std:{}'.format(np.mean(L3_durations),np.std(L3_durations)))

    #~ ax_h[3].hist(L4_durations)
    #~ ax_h[3].set_title('L4 - Mean: {}, Std:{}'.format(np.mean(L4_durations),np.std(L4_durations)))

    #~ ax_h[4].hist(larval_durations)
    #~ ax_h[4].set_xlabel('Time to maturity (hr)')
    #~ ax_h[4].set_title('Maturity (n={}) - Mean: {}, Std:{}'.format(np.size(larval_durations),np.mean(larval_durations),np.std(larval_durations)))

#~ def quick_plot_lifespan(ann_fps, expt_mds, bad_worm_kws=[],debug_mode=False):
    #~ my_ann_files = [annotation_file.AnnotationFile(ann_fp) for ann_fp in ann_fps]

    #~ timestamped_data = {}
    #~ [timestamped_data.setdefault(expt_key,np.array([])) for expt_key in my_ann_files[0].get_data_keys()]
    #~ for [expt_md_fp, ann_file] in zip(expt_mds, my_ann_files):
        #~ if debug_mode: print(expt_md_fp)
        #~ ann_file_data = ann_file.data_as_timestamps_simple(expt_md_fp)
        #~ for expt_key in timestamped_data.keys():
            #~ timestamped_data[expt_key] = np.append(timestamped_data[expt_key],ann_file_data[expt_key])

    #~ viable_worm = (timestamped_data['Hatch']!=-1) \
        #~ & (timestamped_data['Death']!=-1) \
        #~ & np.array([not any([kw in note for kw in bad_worm_kws]) for note in timestamped_data['Notes']])
    #~ print(timestamped_data['Worm'][viable_worm])

    #~ lifespan = (timestamped_data['Death']-timestamped_data['Hatch'])[viable_worm]/(3600*24) # Days
    #~ sorted_ls = sorted(lifespan)
    #~ prop_alive = 1 - (np.arange(start=0,stop=np.size(lifespan)))/np.size(lifespan)

    #~ plt.show()
    #~ plt.ion()

    #~ plt.gcf().clf()
    #~ fig_h, ax_h = plt.subplots(2,1,sharex=True)
    #~ ax_h[0].plot(np.append([0],sorted_ls),np.append([1],prop_alive))
    #~ ax_h[0].set_xlabel('Time since expt. start (d)')
    #~ ax_h[0].set_ylabel('Proportion alive')
    #~ ax_h[0].set_title('Survival curve - n = {}'.format(np.size(lifespan)))

    #~ ax_h[1].hist(lifespan)
    #~ ax_h[1].set_xlabel('Time to death (d)')
    #~ ax_h[1].set_ylabel('Frequency')
    #~ ax_h[1].set_title('Mean+/-STD: {:.2f}+/-{:.2f}d\nMedian:{}d'.format(np.mean(lifespan),np.std(lifespan),np.median(lifespan)))

def clean_plot(my_plot, cleaning_mode=None,**kws):
    if cleaning_mode is None:
        make_labels = kws.get('make_labels',False)
        suppress_ticklabels = kws.get('suppress_ticklabels',False)
    spines_off = kws.get('spines_off',['right','top'])
    num_xticks = kws.get('num_xticks', None)
    num_yticks = kws.get('num_yticks', None)
    square_aspect = kws.get('square_aspect',False)
    ylim = kws.get('ylim', None)
    xlim = kws.get('xlim', None)

    if cleaning_mode == 'visualize':
        make_labels=False
        suppress_ticklabels = False
    elif cleaning_mode == 'PPT':
        make_labels = False
        suppress_ticklabels = True
    elif cleaning_mode == 'verbose':
        make_labels = True
        suppress_ticklabels = False

    [my_plot.spines[spine_line].set_visible(False) for spine_line in spines_off]

    my_plot.tick_params(axis='both',which='both', top='off', bottom='off', left='off', right='off')

    if ylim is not None:
        my_plot.set_ylim(ylim)
    if xlim is not None:
        my_plot.set_xlim(xlim)

    if num_xticks is not None and len(my_plot.get_xticks())>0:
        full_xticks = my_plot.get_xticks().copy()
        full_xticklabels = my_plot.get_xticklabels().copy()
        xlim = my_plot.get_xlim()
        first_xtick = full_xticks[full_xticks>=xlim[0]][0]  # Sometimes there are hidden ticks not visible within limits...
        last_xtick = full_xticks[full_xticks<=xlim[1]][-1]
        if first_xtick*last_xtick < 0:
            my_plot.set_xticks([first_xtick, 0, last_xtick])
        else:
            my_plot.set_xticks([first_xtick, last_xtick])
        if num_xticks > 0:
            my_plot.set_xticks(np.linspace(first_xtick, last_xtick, num_xticks))
        my_plot.xaxis.set_ticks_position('none')

    if num_yticks is not None and len(my_plot.get_yticks())>0:
        full_yticks = my_plot.get_yticks().copy()
        full_yticklabels = my_plot.get_yticklabels().copy()
        ylim = my_plot.get_ylim()
        first_ytick = full_yticks[full_yticks>=ylim[0]][0]
        last_ytick = full_yticks[full_yticks<=ylim[1]][-1]
        if first_ytick*last_ytick < 0:
            my_plot.set_yticks([first_ytick, 0, last_ytick])
        else:
            my_plot.set_yticks([first_ytick, last_ytick])
        if num_yticks > 0:
            my_plot.set_yticks(np.linspace(first_ytick, last_ytick, num_yticks))
        my_plot.yaxis.set_ticks_position('none')

    if suppress_ticklabels:
        my_plot.set_xticks([])
        my_plot.set_yticks([])

    if not make_labels:
        my_plot.set_xlabel('')
        my_plot.set_ylabel('')
        my_plot.set_title('')

    if square_aspect:
        square_plot_aspect(my_plot)

def square_plot_aspect(my_plot):
    x0,x1 = my_plot.get_xlim()
    y0,y1 = my_plot.get_ylim()
    my_plot.set_aspect(abs(x1-x0)/abs(y1-y0))

def save_cleaned_fig(fig_h,ax_h,out_fn,**kws):
    if 'cleaning_mode' not in kws:
        kws['cleaning_mode'] = 'verbose'
    if type(ax_h) in [list, np.ndarray]:
        [clean_plot(my_ax,**kws) for my_ax in utilities.flatten_list(ax_h)]
    else:
        clean_plot(ax_h,**kws)
    fig_h.savefig(out_fn)

def force_same_plot_attributes(my_axes, *args):
    for attr in args:
        if attr=='xlim':
            min_x = min([ax.get_xlim()[0] for ax in my_axes])
            max_x = max([ax.get_xlim()[1] for ax in my_axes])
            [ax.set_xlim([min_x,max_x]) for ax in my_axes]
        elif attr=='ylim':
            min_y = min([ax.get_ylim()[0] for ax in my_axes])
            max_y = max([ax.get_ylim()[1] for ax in my_axes])
            [ax.set_ylim([min_y,max_y]) for ax in my_axes]
