import pathlib

import numpy
import matplotlib as mpl
import matplotlib.colors as mpl_colors

from zplib.image import colorize
from elegant import load_data

import plotting_tools
from utilities import utilities


# def draw_position_map(worms, *features, experiment_dir=None):
#     '''
#     feature - callable returning true/false
#     '''
#     import matplotlib.pyplot as plt; plt.ion(); plt.show()
#
#     if experiment_dir:
#         annotations = load_data.read_annotations(experiment_dir)
#         included_worms = [worm.name.split()[-1] for worm in worms]
#         excluded = [position for position in annotations.keys() if position not in included_worms]
#
#     fig_h, ax_h = plt.subplots()
#     grouped_worms = []
#     for feature in feature:
#         grouped_worms.append(worms.filter(feature))
#     categorized_worms
#
#     big_worms = worms.filter(lambda worm: life_quality_stats.peak_size('length')(worm) > 1200)
#     small_worms = worms.filter(lambda worm: worm not in big_worms)
#     ref_x = min([worm.stage_x for worm in worms])
#
#     for worm in big_worms:
#         ax_h.add_artist(plt.Circle((worm.stage_x-ref_x, worm.stage_y), 1, color = color_cycle[1]))
#
#     for worm in small_worms:
#         ax_h.add_artist(plt.Circle((worm.stage_x-ref_x, worm.stage_y), 1, color = color_cycle[2]))
#     ax_h.set_xlim([-5,20])
#     ax_h.set_ylim([0,65])
#     ax_h.set_aspect('equal')
#     ax_h.set_title(f'{worms[0].name.split()[0]}')

def get_worms_without_feature(worms, feature):
    missing_feature = [worm for worm in worms if feature not in dir(worm.td)]
    null_feature = [worm for worm in worms
        if (worm not in missing_feature and ~numpy.isnan(worm.get_feature(feature)))
    ]

    return missing_feature, null_feature

def enumerate_common_annotations(experiment_dir,bad_kws=None,verbose=True, filter_good=True):
    if not bad_kws:
        bad_kws = [
            ['Nw','Nh'],
            ['REFERENCE'],
            ['CONTAMINATION'],
            ['DOUBLE WORM', 'TRIPLE WORM'],
            ['ESCAPE', 'VISITED'],
            ['LOST'],
            ['LARVAL', 'DELAYED'],
            ['FERTILE'],
            ['PVL', 'BURST',],
            ['bag\'d'],
            ['small','sickly','scrawny','mottled']
        ]

    experiment_dir = pathlib.Path(experiment_dir)

    inventory_bad_worms = []
    annotations = load_data.read_annotations(experiment_dir)
    if filter_good:
        original_worms = annotations.keys()
        annotations = load_data.filter_annotations(annotations, load_data.filter_excluded)
        excluded_worms = sorted(list(set(original_worms) - set(annotations.keys())))

    for kw_group in bad_kws:
        group_list = []
        for kw in kw_group:
            group_list.extend([worm for worm, worm_annotations in annotations.items() if kw in worm_annotations[0]['notes']])

        group_list = utilities.unique_items(group_list)
        inventory_bad_worms.append(group_list)

    if verbose:
        print(f'\n{experiment_dir.name} (n = {len(annotations)})')
        if filter_good:
            print(f'(excluded): {len(excluded_worms)} ({excluded_worms})')

        for kw_group, bad_worms in zip(bad_kws, inventory_bad_worms):
            print(f'{"/".join(kw_group)}: {len(bad_worms)} ({bad_worms})')
    utilities.print_table(
        [[len(bad_worms) for bad_worms in inventory_bad_worms] + [f'{len(annotations)}']],
        column_names = [f'{"/".join(kw_group)}' for kw_group in bad_kws] + ['Total'],
        row_names = [experiment_dir.name]
    )

def spatial_distribution(worms, feature, ax_h=None, fig_h=None, source='gradient'):
    import matplotlib.pyplot as plt

    INTERSLOT_SPACING = 30 # spacing between slots on stage insert (mm)

    ax_provided = ax_h is not None
    if not ax_provided:
        fig_h, ax_h = plt.subplots()
        ax_h.set_title(f'{worms[0].name.split()[0]}')
    else:
        assert fig_h is not None

    feature_vals = worms.get_feature(feature)
    nan_mask = numpy.isnan(feature_vals)
    if all(nan_mask):
        print(f'feature not found for any individuals')
        return
    elif any(nan_mask):
        print(f'Warning: nan\'s found for {numpy.array(worms)[nan_mask]}')

    good_worms = numpy.array(worms)[~nan_mask]
    feature_vals = feature_vals[~nan_mask]

    # color_vals = colorize.scale(feature_vals,output_max=1,gamma=0.5)
    # color_vals = colorize.scale(feature_vals,output_max=1,gamma=0.5)
    # color_vals = plotting_tools.build_gradient_palette(
    #     base_color = numpy.array(plotting_tools.qual_colors[numpy.random.randint(len(plotting_tools.qual_colors))]),
    #     num_colors=len(feature_vals)
    # )
    if source == 'gradient':
        color_vals = plotting_tools.build_gradient_palette(
            base_color = numpy.array(plotting_tools.qual_colors[numpy.random.randint(len(plotting_tools.qual_colors))]),
            num_colors=256
        )
    elif source == 'zplib':
        color_vals = colorize.color_map(numpy.linspace(0,1,257))

    median_x = numpy.median([worm.stage_x for worm in good_worms])
    ref_x = (median_x // INTERSLOT_SPACING) * INTERSLOT_SPACING # Adjusts to the nearest slide position

    cmap = mpl_colors.ListedColormap([numpy.append(val, 1) for val in color_vals])
    norm = mpl.colors.Normalize(vmin=feature_vals.min(),vmax=feature_vals.max())
    sm = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
    sm.set_array([])

    rect = ax_h.add_artist(plt.Rectangle((-5,-15),35,90,linewidth=1,edgecolor='k',facecolor='none'))
    # for worm, color in zip(good_worms,color_vals):
    #     ax_h.add_artist(
    #         plt.Circle((worm.stage_x-ref_x+(median_x-ref_x), worm.stage_y), 1, edgecolor = [0,0,0], facecolor=color))

    for worm, feature_val in zip(good_worms, feature_vals):
        ax_h.add_artist(
            plt.Circle((worm.stage_x-ref_x+(median_x-ref_x), worm.stage_y), 1, edgecolor = [0,0,0], facecolor=sm.to_rgba(feature_val)))


    for worm in worms.filter(lambda worm: worm not in good_worms):
        ax_h.add_artist(
            plt.Circle((worm.stage_x-ref_x+(median_x-ref_x), worm.stage_y), 1, edgecolor = [0,0,0], facecolor=[0,0,0]))


    ax_h.set_xlim([-5,30])
    ax_h.set_ylim([-15,75])
    ax_h.set_xlim([-10,35])
    ax_h.set_ylim([-20,80])
    ax_h.spines['bottom'].set_visible(False)
    ax_h.spines['left'].set_visible(False)

    ax_h.set_aspect('equal')

    # Show the colorbar for these values
    fig_h.colorbar(sm,ax=ax_h)

    if not ax_provided:
        return fig_h, ax_h
