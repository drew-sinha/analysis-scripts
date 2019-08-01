from PyQt5 import Qt
from ris_widget.qwidgets import annotator
import pathlib
import sys
from ris_widget import ris_widget
import time
from concurrent import futures
import collections


"""
READ ME
 T/F Annotator for eggs given a file annotations.txt and folder of egg images
"""

class BasicEggAnnotationField(annotator.AnnotationField):

    COLORS = [(184, 255, 184), (255, 184, 184), (255, 255, 184)]
    CATEGORIES = ['egg','without_egg', 'questionable']
    SHORTCUTS = [category[0] for category in CATEGORIES]
    COLOR_MAP = {}; COLOR_MAP.update(zip(CATEGORIES, COLORS))
    assert COLOR_MAP

    def __init__(self):
        name = 'is_egg'
        super().__init__(name)

    def init_widget(self):
        self.widget = Qt.QGroupBox(self.name)
        layout = Qt.QHBoxLayout()
        self.widget.setLayout(layout)

        self.label = Qt.QLabel()
        layout.addWidget(self.label)

        for category, key in zip(self.CATEGORIES, self.SHORTCUTS):
            button = Qt.QPushButton(category)
            callback = self._make_annotation_callback(category)
            button.clicked.connect(callback)
            layout.addWidget(button)
            Qt.QShortcut(key, self.widget, callback, context=Qt.Qt.ApplicationShortcut)

    def _make_annotation_callback(self, category):
        def callback():
            self.update_annotation(category)
            self.update_widget(category)
            self.flipbook.focus_next_page()
        return callback

    def update_widget(self, value):
        if value is None:
            self.label.setText('')
        elif value not in self.CATEGORIES:
            raise ValueError(f'Value {value} not in defined categories')
        else:
            self.label.setText(value)

        for page in self.flipbook.pages:
            egg_annotation = self.get_annotation(page)
            page.color = self.COLOR_MAP[egg_annotation] if egg_annotation else None 

class BasicEggAnnotator(annotator.Annotator):

    def __init__(self, ris_widget, fields, image_directory, batch_size=200,annotation_filename='annotations'):
        super().__init__(ris_widget, fields)
        widget = Qt.QGroupBox()
        annotation_layout = Qt.QVBoxLayout(widget)

        ris_widget.add_annotator(fields)
        self.ris_widget = ris_widget

        self.image_directory = pathlib.Path(image_directory)
        # image_list = list(self.image_directory.glob('*.png'))
        image_list = list(self.image_directory.glob('*_*.png'))
        self.images = sorted(
            image_list, 
            key=lambda entry:int(entry.name.split('_')[0])) # Sort by first field (alt. just rename to have constant number of zeros at beginning)
        assert self.images

        self.annotation_file = self.image_directory / f'{annotation_filename}.txt'
        self.annotations = collections.defaultdict(dict)
        self.load_annotations()

        self.batch_size = batch_size
        self.batch_num = 0

        info_layout = Qt.QHBoxLayout()
        info_layout.setSpacing(11)
        save = self._add_button(info_layout, 'Save', self.save_annotations)
        annotation_layout.addLayout(info_layout)

        nav_buttons = Qt.QHBoxLayout()
        nav_buttons.setSpacing(11)
        self._prev_button = self._add_button(nav_buttons, '\N{LEFTWARDS ARROW TO BAR}', lambda: self.load_next_batch(-1))
        self._add_button(nav_buttons, '\N{UPWARDS ARROW}', self.prev_timepoint)
        self._add_button(nav_buttons, '\N{DOWNWARDS ARROW}', self.next_timepoint)
        self._next_button = self._add_button(nav_buttons, '\N{RIGHTWARDS ARROW TO BAR}', lambda: self.load_next_batch(1))
        annotation_layout.addLayout(nav_buttons)
        ris_widget.annotator.layout().insertRow(0,widget)
        self.load_next_batch(0)


    def load_annotations(self):
        if self.annotation_file.exists():
            with self.annotation_file.open('r') as annotation_file_pointer:
                header = annotation_file_pointer.__next__()
                for line in annotation_file_pointer:
                    tokens = line.split('\n')[0].split('\t')
                    self.annotations[tokens[0]]['is_egg'] = tokens[1] if tokens[1] else None

    def _add_button(self, layout, title, callback):
        button = Qt.QPushButton(title)
        button.clicked.connect(callback)
        layout.addWidget(button)
        return button

    def load_next_batch(self, increment):
        self.pull_annotations_from_flipbook()
        if increment != 0:
            self.save_annotations()

        next_batch_num = self.batch_num + increment
        start_idx = self.batch_size * next_batch_num
        stop_idx = self.batch_size * (next_batch_num+1)
        if stop_idx > len(self.images):
            stop_idx = len(self.images)

        self.ris_widget.flipbook_pages.clear()
        self.ris_widget.add_image_files_to_flipbook(self.images[start_idx:stop_idx])
        self.add_annotations_to_flipbook()

        self.batch_num = next_batch_num
        self._prev_button.setEnabled(self.batch_num != 0)
        self._next_button.setEnabled((self.batch_num + 1) * self.batch_size < len(self.images))

    def prev_timepoint(self):
        self.flipbook.focus_prev_page()

    def next_timepoint(self):
        self.flipbook.focus_next_page()

    def pull_annotations_from_flipbook(self):
        for page in self.flipbook.pages:
            self.annotations[page._abbreviated_name]['is_egg'] = self.fields[0].get_annotation(page)

    def add_annotations_to_flipbook(self):
        for page in self.ris_widget.flipbook_pages:
            abbreviated_page_name = page.name.split('.')[0]
            page._abbreviated_name = abbreviated_page_name
            # self.flipbook.pages.annotations = self.annotations.get(abbreviated_page_name,{}).copy()
            page.annotations = self.annotations.get(abbreviated_page_name,{}).copy()
        self.update_fields()

    def save_annotations(self):
        self.pull_annotations_from_flipbook()

        # names = self.annotations.keys()
        # sorted_names = sorted(
        #     names, 
        #     key=lambda entry:int(entry.name.split('_')[0]))
        # self.annotations = {name: self.annotations[name] for name in sorted_names}

        with self.annotation_file.open('w+') as annotation_file_pointer:
            annotation_file_pointer.write('\t'.join(['Image Name', 'is_egg']) + '\n')
            for image_name, annotation in self.annotations.items():
                if annotation['is_egg']:
                    annotation_file_pointer.write('\t'.join([image_name, annotation['is_egg']]) + '\n')

if __name__ == "__main__":
    '''Call as python egg_annotator.py IMAGE_DIRECTORY'''
    image_directory = pathlib.Path(sys.argv[1])
    annotations = {}

    try:
        rw
    except NameError:
        rw = ris_widget.RisWidget()
    egg_annotator = BasicEggAnnotator(rw, [BasicEggAnnotationField()], image_directory)
