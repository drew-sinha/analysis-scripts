import pathlib
from concurrent import futures

import numpy
from PyQt5 import Qt

import freeimage
from ris_widget import ris_widget
from ris_widget import image as rw_image
from ris_widget.qwidgets import flipbook
from zplib.image import mask as zplib_image_mask

class RWPainter:
    """A painter that loads images from a directory into a RisWidget object and draw overlays that can be saved"""

    def __init__(self, rw, image_dir, out_dir):
        """
        Parameters
            rw - RisWidget object to load images into
            image_dir - str/pathlib.Path to directory images for loading
            out_dir - str/pathlib.Path to directory for saving overlays
        """

        self.rw = rw
        self.image_dir = pathlib.Path(image_dir)
        self.out_dir = pathlib.Path(out_dir)
        assert self.image_dir.exists(); assert self.out_dir.exists()

        self.editing=False

        # To be low-maintenance, just add everything to the flipbook layout
        layout = Qt.QFormLayout()
        widget = Qt.QWidget()
        widget.setLayout(layout)
        self.rw.flipbook.layout().addWidget(widget)

        self.rw.add_painter()
        self.rw.flipbook_pages.clear()

        self.image_list = sorted(list(self.image_dir.glob('*.png')))
        for image_path in self.image_list:
            image = freeimage.read(image_path)
            image_list = flipbook.ImageList()
            image_list.append(rw_image.Image(data=image,name=image_path.name))
            if (self.out_dir / image_path.name).exists():
                image_list.append(freeimage.read((self.out_dir / image_path.name)))
            else:
                image_shape = image.shape
                if len(image_shape) == 3: 
                    image_list.append(numpy.zeros_like(image))
                else:
                    new_shape = list(image_shape) + [3]
                    image_list.append(numpy.zeros(shape=new_shape).astype('uint8'))
            rw.flipbook_pages.append(image_list)

        self.rw.painter.brush_size.value = 13

        self.clear = self._add_button(layout, 'Clear All', self._on_clear_clicked)
        self.reload = self._add_button(layout, 'Reload Overlay', self._on_reload_clicked)
        self.save = self._add_button(layout, 'Save Overlay', self._on_save_clicked)

    def _add_button(self, layout, title, callback):
        button = Qt.QPushButton(title)
        button.clicked.connect(callback)
        layout.addWidget(button)
        return button

    def _on_clear_clicked(self):
        image_data = rw.layers[1].image.data
        image_data *= 0
        rw.layers[1].image.refresh()

    def _on_reload_clicked(self):
        if (self.out_dir / rw.layers[0].image.name).exists():
            original_overlay = freeimage.read(self.out_dir / image_path.name)
            rw.layers[1].image.data = original_overlay
        else:
            print('No original overlay exists')

    def _on_save_clicked(self):
        overlay = self.rw.layers[1].image.data
        if overlay.any():
            working_file = self.out_dir / pathlib.Path(self.rw.layers[0].image.name).name
            freeimage.write(overlay, working_file)

if __name__ == "__main__":
    image_dir = pathlib.Path('/Volumes/zplab/acquired_data/20190901_cholesterol_supp_expt/20190907_alternatelighting/DS/Fridge_food+chol')
    image_dir = pathlib.Path('/Volumes/zplab/acquired_data/20190904_freshfood_expt/20190910_alternatelighting/DS/Fresh_food_DS_annotations')
    if not image_dir.exists():
        raise Exception('image directory doesn\'t exist!')
    # out_dir = image_dir.parent / (image_dir.name + '_centerlines')
    suffix = 'WP'
    out_dir = image_dir / ('_centerlines' + suffix)

    out_dir.mkdir(exist_ok=True)

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()

    rwp = RWPainter(rw, image_dir,out_dir)
