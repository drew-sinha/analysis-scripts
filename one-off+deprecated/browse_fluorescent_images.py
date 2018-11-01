
import pathlib

import numpy as np
import matplotlib.pyplot as plt
plt.ion(); plt.show()

import freeimage
from zplib.scalar_stats import kde
from ris_widget import ris_widget
from elegant import load_data

import plotting_tools

if __name__ == "__main__":
    view_files = False

    experiment_dir = pathlib.Path('/mnt/9karray/Sinha_Drew/20180518_spe-9_Run_3/')
    timepoint = '2018-05-22t1716' # Early timepoint
    timepoint = '2018-05-28t1710' # Late timepoint ~8d - eggs are super autofluorescent!!!!! churned up food doesn't help either....

    experiment_dir = pathlib.Path('/mnt/9karray/Sinha_Drew/20180524_spe-9_Run_4/')
    timepoint = '2018-05-28t1526' # Early timepoint
    timepoint = '2018-06-03t1714' # Late timepoint ~8d ..

    if view_files:
        try:
            rw
        except NameError:
            rw = ris_widget.RisWidget()
            rw.show()

    def timepoint_filter(position_name, timepoint_name):
        return timepoint_name == timepoint
    experiment_images = load_data.scan_experiment_dir(experiment_dir,timepoint_filter=timepoint_filter, channels = ['bf','green_yellow_excitation_autofluorescence'])

    image_names, images = [], []
    fl_image_modes = [] # Measurement of the gel
    fl_image_center_medians = [] # Measurement of the food (1100,900) to (1300, 1100)
    for position_name, position_images in experiment_images.items():
        if position_images[timepoint]:
            processed_position_images = []
            for image_name in position_images[timepoint]:
                raw_image = freeimage.read(str(image_name))
                if 'fl' in image_name.name:
                    fl_image_modes.append(np.bincount(raw_image.flat)[1:].argmax()+1)
                    fl_image_center_medians.append(np.median(raw_image[1100:1300,900:1100]))

                if view_files:
                    if 'fl' in image_name.name.split()[1]:
                        flatfield_name = experiment_dir / 'calibrations' / (timepoint + ' fl_flatfield.tiff')
                    elif 'bf.png' in image_name.name:
                        flatfield_name = experiment_dir / 'calibrations' / (timepoint + ' bf_flatfield.tiff')
                    flatfield = freeimage.read(str(flatfield_name))
                    image = raw_image.astype('float32') * flatfield
                    processed_position_images.append(image)

                    # Mode-centering
                    # noise_floor = np.mean(raw_image[0:51,0:51])
                    # mode = numpy.bincount(bf.flat)[1:].argmax()+1
                    # image -= noise_floor
                    # image *= (300-noise_floor) / (mode-noise_floor)

                    # Consider dF/F


            if view_files: images.append(processed_position_images)
    if view_files: rw.flipbook_pages = images

    fig_h, ax_h = plt.subplots()
    ax_h.set_title(f'''{experiment_dir.name}: {timepoint}
        Mode (Background) = {np.median(fl_image_modes)} +/- {np.std(fl_image_modes)}
        Center Median (Pad) = {np.median(fl_image_center_medians)} +/- {np.diff(np.percentile(fl_image_center_medians,[25, 75]))}''')
    for i, vals in enumerate([fl_image_modes, fl_image_center_medians]):
        x,y,e = kde.kd_distribution(vals)
        ax_h.plot(x,y,color=plotting_tools.qual_colors[i])
    ax_h.legend(['mode', 'center_medians'])
