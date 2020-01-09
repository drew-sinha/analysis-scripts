import pathlib
import pickle

from PyQt5 import Qt

from ris_widget.qwidgets import annotator

import rw_funs
import spline_set

"""
This module defines an annotator for drawing one or more splines on multiple sets of images for a given directory
and saving that spline data in an associated pickle file.

Instructions:
    Replace the __main__ clause at the bottom with the path of interest for the image directory.
    Run this module as a script using python/ipython.
    Press "D" or the button in the annotator to start drawing a spline. Mouse up finishes drawing the spline.
    Clicking the "Save Annotations" button saves splines drawn **on all images** to a pickle file in the image directory.
    Clicking the "Reload Annotations" button restores the splines from the original pickle file used to load them.
    Press "T"  to initiate a zoom-to-fiT and reset the image zoom
"""


'''
    TODO:
        Accumulation of annotator fields in background (propensity for segfaults?)
        Getting rid of images (duplicates/bad images) after the fact
'''

def _add_button(layout, title, callback):
    button = Qt.QPushButton(title)
    button.clicked.connect(callback)
    layout.addWidget(button)
    return button

class MultisplineAnnotationField(annotator.AnnotationField):
    def __init__(self, ris_widget, default=[], color=(255, 0, 0)):
        self.ris_widget = ris_widget
        pen = Qt.QPen(Qt.QColor(*color))
        self.spline_set = spline_set.SplineSet(ris_widget, pen)
        self.spline_set.geometry_change_callbacks.append(self.on_geometry_change)
        super().__init__(name='MultisplineAnnotation', default=default)

    def init_widget(self):
        self.widget = Qt.QGroupBox(self.name)
        layout = Qt.QHBoxLayout()
        self.widget.setLayout(layout)
        self.label = Qt.QLabel()
        layout.addWidget(self.label)
        self.draw = _add_button(layout, 'Draw Spline', self._handle_spline_drawing)
        self.draw.setEnabled(True)
        Qt.QShortcut(Qt.Qt.Key_D, self.widget, self._handle_spline_drawing, context=Qt.Qt.ApplicationShortcut)

        self.show_spline = Qt.QCheckBox('Show Splines')
        self.show_spline.setChecked(True)
        self.show_spline.toggled.connect(self.toggle_show_spline)
        layout.addWidget(self.show_spline)

    def toggle_show_spline(self, show):
        [spline.setVisible(show) for spline in self.spline_set.splines]
        self.draw.setEnabled(show) # Assume we're not trying to draw a spline if we're hiding all of them..

    def _handle_spline_drawing(self):
        self.spline_set.setup_next_spline()
        self.draw.setText('Drawing...')
        self.draw.setEnabled(False)
        self.show_spline.setEnabled(False) # Don't allow disappearing splines while drawing.

    def on_geometry_change(self, spline_tcks):
        self.update_annotation(spline_tcks)
        if not self.spline_set.splines or not (self.spline_set.splines[-1].warping):
            self.draw.setText('Draw Spline')
            self.draw.setEnabled(True)
            self.show_spline.setEnabled(True)

    def update_widget(self, spline_tcks):
        self.spline_set.geometry = spline_tcks
        [spline.setVisible(self.show_spline.isChecked()) for spline in self.spline_set.splines]

class MultisplineAnnotator:
    """A spline annotator that loads images from a directory and draw splines that can be saved"""

    def __init__(self, rw, fields, image_dir, annotation_file='annotations.pickle'):
        self.rw = rw
        self.rw.add_annotator(fields)
        self.image_dir = pathlib.Path(image_dir)
        self.annotation_file = annotation_file
        assert self.image_dir.exists()

        self.image_paths = sorted(self.image_dir.glob('*.png'))
        self.rw.add_image_files_to_flipbook(
            self.image_paths,
            page_names=[path.stem for path in self.image_paths]
        )

        self.load_annotations()

        nb_widget = Qt.QGroupBox()
        button_layout = Qt.QHBoxLayout(nb_widget)
        button_layout.setSpacing(11)
        up_button = _add_button(button_layout, '\N{UPWARDS ARROW}', lambda: self.rw.flipbook.focus_prev_page())
        down_button = _add_button(button_layout, '\N{DOWNWARDS ARROW}', lambda: self.rw.flipbook.focus_next_page())
        self.rw.annotator.layout().insertRow(0, nb_widget)

        db_widget = Qt.QGroupBox()
        data_buttons = Qt.QVBoxLayout(db_widget)
        save = _add_button(data_buttons, 'Save Annotations', self.save_annotations)
        reload = _add_button(data_buttons, 'Reload Annotations', self.load_annotations)
        self.rw.annotator.layout().insertRow(1, db_widget)

        Qt.QShortcut(Qt.Qt.Key_T, self.rw.annotator, self.toggle_zoom_to_fit, context=Qt.Qt.ApplicationShortcut)
        Qt.QShortcut(Qt.Qt.Key_Q, self.rw.annotator, lambda: self.incremental_zoom(False), context=Qt.Qt.ApplicationShortcut)
        Qt.QShortcut(Qt.Qt.Key_W, self.rw.annotator, lambda: self.incremental_zoom(True), context=Qt.Qt.ApplicationShortcut)

    def toggle_zoom_to_fit(self):
        self.rw.qt_object.image_view.zoom_to_fit = True

    def incremental_zoom(self, zoom_in):
        self.rw.image_view.change_zoom(zoom_in)
        self.rw.image_view.zoom_to_fit = False

    def load_annotations(self):
        try:
            with (self.image_dir / self.annotation_file).open('rb') as annotation_fp:
                self.rw.annotator.all_annotations = pickle.load(annotation_fp)
        except FileNotFoundError:
            self.all_annotations = {}

    def save_annotations(self):
        # Sometimes, accidentally swiping and creates a bad invisible tck (with a None value).
        # Quickly check the annotations to remove any.
        all_annotations = self.rw.annotator.all_annotations
        replacement_annotations = []
        for page_annotations in all_annotations:
            replacement_annotations.append(page_annotations)
            good_annotations = [tck for tck in page_annotations['MultisplineAnnotation'] if tck is not None]
            replacement_annotations[-1]['MultisplineAnnotation'] = good_annotations
        all_annotations = replacement_annotations

        with (self.image_dir / self.annotation_file).open('wb') as annotation_fp:
            pickle.dump(all_annotations, annotation_fp)

if __name__ == "__main__":
    image_dir = pathlib.Path('/mnt/fluoro-scope/acquired_data/20200102_2dayRNAseqRun/20200106/Corral_Food_1')

    if not image_dir.exists():
        raise Exception('image directory doesn\'t exist!')
    annotation_file = 'annotations_DS.pickle'

    rw = rw_funs.make_global_riswidget()
    if hasattr(rw, 'annotator'):
        rw.annotator.close()
        del(rw.annotator)
    rw.show()
    rw.flipbook_pages.clear()
    fields = [MultisplineAnnotationField(rw)]
    msa = MultisplineAnnotator(rw, fields, image_dir,annotation_file=annotation_file)
