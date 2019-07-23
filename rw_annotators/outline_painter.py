import pathlib
from concurrent import futures

import numpy
from PyQt5 import Qt

import freeimage
from ris_widget import ris_widget
from ris_widget import image as rw_image
from ris_widget.qwidgets import annotator, flipbook
from elegant import worm_spline
from zplib.image import mask as zplib_image_mask
from zplib.curve import spline_geometry

class OutlinePainter:
    def __init__(self, rw, image_dir, out_dir):
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
        self.rw.painter.hide()
        self.rw.painter.setVisible(False)
        self.rw.painter.painter_item.hide()

        self.rw.flipbook_pages.clear()

        self.image_list = sorted(list(self.image_dir.glob('*.png')))
        for image_path in self.image_list:
            image = freeimage.read(image_path)
            image_list = flipbook.ImageList()
            image_list.append(rw_image.Image(data=image,name=image_path.name))
            if (self.out_dir / image_path.name).exists():
                image_list.append(freeimage.read((self.out_dir / image_path.name)))
            else:
                image_list.append(numpy.zeros_like(image))
            rw.flipbook_pages.append(image_list)

        self.edit = self._add_button(layout, 'Edit Outlines', self._on_edit_clicked)

    def _add_button(self, layout, title, callback):
        button = Qt.QPushButton(title)
        button.clicked.connect(callback)
        layout.addWidget(button)
        return button

    def _on_edit_clicked(self):
        if self.editing:
            self.stop_editing()
        else:
            self.start_editing()

    def start_editing(self):
        self.editing = True
        self.edit.setText('Save Edits')
        self.rw.painter.show()
        self.rw.painter.painter_item.show()

        # Bring focus to mask layer
        sm = self.rw.qt_object.layer_stack._selection_model
        m = sm.model()
        sm.setCurrentIndex(m.index(0,0),
            Qt.QItemSelectionModel.SelectCurrent|Qt.QItemSelectionModel.Rows)

        self.rw.painter.brush_size.value = 10
        self.rw.painter.brush_val.value = (255,0,0)
        self.rw.layers[1].opacity = 1

    def stop_editing(self):
        outline = self.rw.layers[1].image.data
        if outline.any():
            working_file = self.out_dir / pathlib.Path(self.rw.layers[0].image.name).name
            freeimage.write(outline, working_file)

        self.rw.layers[1].opacity = 0.75

        self.rw.painter.hide()
        self.rw.painter.painter_item.hide()
        self.editing = False
        self.edit.setText('Edit Outlines')

def parse_outline(outline_mask):
    masks = []
    for layer in numpy.moveaxis(outline_mask,-1,0):
        while layer.any():
            outline = zplib_image_mask.get_largest_object(layer>0)
            new_mask = zplib_image_mask.fill_small_area_holes(outline,300000).astype('uint8')
            new_mask[new_mask>0] = -1
            masks.append(new_mask)
    return masks

def process_outline_dir(source_dir, microns_per_pixel, out_dir=None):
    if out_dir is None:
        out_dir = source_dir

    mask_data = {}
    mask_data_entries = ['mask_names', 'mask_areas', 'mask_lengths']
    for outline_image_path in pathlib.Path(source_dir).iterdir():
        if outline_image_path.suffix[1:] != 'png': continue
        outline_image = freeimage.load(outline_image_path)
        masks = parse_outline(outline_image)
        for mask_num, mask in enumerate(masks):
            freeimage.write(mask, out_dir / (outline_image_path.stem + f'_{mask_num}.png'))
            center_tck, width_tck = worm_spline.pose_from_mask(mask)
            length = spline_geometry.arc_length(center_tck) * self.microns_per_pixel

            mask_data.setdefault('mask_names',[]).append(outline_image_path.stem + f'_{mask_num}')
            mask_data.setdefault('mask_areas',[]).append(mask.sum()*(microns_per_pixel/1000)**2)
            mask_data.setdefault('mask_lengths',[]).append(length*microns_per_pixel/1000)

        with (self.out_dir / 'measurements.txt').open('w+') as mf_pointer:
            mf_pointer.write('\t'.join(mask_data_entries))
            # for

if __name__ == "__main__":
    image_dir = '/home/drew/20190705/concentrated_OP50'
    out_dir = pathlib.Path('/home/drew/20190705/concentrated_OP50_outlines')

    out_dir.mkdir(exist_ok=True)

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()

    op = OutlinePainter(rw, image_dir,out_dir)
