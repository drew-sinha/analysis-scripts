import numpy as np
import matplotlib.pyplot as plt

def plot_trajectories(strain_df, health_var='health', worms_to_plot=None, num_worms=20, worm_subset = None):
    if worms_to_plot is None:
        if worm_subset is not None:
            worms_to_plot = np.random.permutation(strain_df.worms)[worm_subset][:num_worms]
        else:
            worms_to_plot = np.random.permutation(strain_df.worms)[:num_worms]
    elif type(worms_to_plot[0]) is not str:
        worms_to_plot = np.array(strain_df.worms)[worms_to_plot]

    fig_h, ax_h = plt.subplots()
    for worm in worms_to_plot:
        worm_health = strain_df.display_variables(
            strain_df.mloc(worms=[worm],measures=[health_var])[0,0,:],
            my_var=health_var)[0]
        worm_health = worm_health[~np.isnan(worm_health)]

        ax_h.plot(np.array(strain_df.ages)[:len(worm_health)], worm_health)
    if health_var in ['movement','texture','autofluorescence','eggs','size']:
        ax_h.set_ylabel('{} health across life'.format(health_var))
    elif health_var == 'health':
        ax_h.set_ylabel('Overall health across life')
    else:
        ax_h.set_ylabel('{} across life'.format(health_var))
    ax_h.set_xlabel('Time Post-Adulthood (days)')
    return fig_h, ax_h
