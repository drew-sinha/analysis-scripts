import pathlib
import pickle

from PyQt5 import Qt

from ris_widget.qwidgets import annotator

import rw_funs
import spline_set

'''
    TODO:
        Need tool to update legacy annotations to this new format.
        Functionality for hiding centerlines.
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


    def _handle_spline_drawing(self):
        self.spline_set.setup_next_spline()
        self.draw.setText('Drawing...')
        self.draw.setEnabled(False)

    def on_geometry_change(self, spline_tcks):
        self.update_annotation(spline_tcks)
        if not self.spline_set.splines or not (self.spline_set.splines[-1].warping):
            self.draw.setText('Draw Spline')
            self.draw.setEnabled(True)

    def update_widget(self, spline_tcks):
        self.spline_set.geometry = spline_tcks

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

        widget = Qt.QGroupBox()
        annotation_layout = Qt.QVBoxLayout(widget)
        save = _add_button(annotation_layout, 'Save Annotations', self.save_annotations)
        reload = _add_button(annotation_layout, 'Reload Annotations', self.load_annotations)
        self.rw.annotator.layout().insertRow(0,widget)

    def load_annotations(self):
        try:
            with (self.image_dir / self.annotation_file).open('rb') as annotation_fp:
                self.rw.annotator.all_annotations = pickle.load(annotation_fp)
        except FileNotFoundError:
            self.all_annotations = {}

    def save_annotations(self):
        with (self.image_dir / self.annotation_file).open('wb') as annotation_fp:
            pickle.dump(self.rw.annotator.all_annotations, annotation_fp)

if __name__ == "__main__":
    image_dir = pathlib.Path('/mnt/fluoro-scope/acquired_data/20191010_25ul_FridgeDropExperiment/20191010/Lawn_Ctrl_3')

    if not image_dir.exists():
        raise Exception('image directory doesn\'t exist!')
    annotation_file = 'annotations_DS.pickle'

    rw = rw_funs.make_global_riswidget()
    if hasattr(rw, 'annotator'):
        rw.annotator.close()
        del(rw.annotator)
    rw.show()
    rw.flipbook_pages.clear()
    fields = [MultisplineAnnotationField(rw,annotations_file=annotations_file)]
    msa = MultisplineAnnotator(rw, fields, image_dir)
