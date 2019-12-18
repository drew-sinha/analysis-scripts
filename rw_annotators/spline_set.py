
from PyQt5 import Qt

from ris_widget.overlay import base
from elegant.gui.spline_overlay import center_spline

class SplineSet(base.RWGeometryItemMixin, Qt.QGraphicsPathItem):
    """Defines a container for a set of (center)splines"""

    def __init__(self, ris_widget, pen=None, geometry=None):
        self.splines = []
        self._last_click_deselected=False
        self.pen = pen
        super().__init__(ris_widget,geometry=geometry)

    @property
    def geometry(self):
        if len(self.splines) == 0:
            geometry = None
        else:
            geometry = [spline.geometry for spline in self.splines]
        return geometry

    @geometry.setter
    def geometry(self, geometry):
        [spline.remove() for spline in self.splines]
        self.splines = []
        if geometry:
            for tck in geometry:
                self.splines.append(
                    center_spline.CenterSpline(self.ris_widget, geometry=tck)
                )
                self.splines[-1].geometry_change_callbacks.append(self._update_spline_in_set)

    def remove(self):
        if self.splines:
            [spline.remove() for spline in self.splines]
        super().remove()

    def setup_next_spline(self):
        next_cs = center_spline.CenterSpline(self.ris_widget)
        next_cs.geometry_change_callbacks.append(self._update_spline_in_set)
        self.splines.append(next_cs)
        next_cs.start_drawing()

    def _update_spline_in_set(self, tck):
        if tck:
            self._geometry_changed()
        else:
            # Handle a deletion
            new_splines = []
            deleted = False
            for spline in self.splines:
                if spline.geometry is not None:
                    new_splines.append(spline)
                else:
                    spline.remove()
                    deleted = True
            if deleted:
                self.splines = new_splines
                self._geometry_changed()

    def sceneEventFilter(self, watched, event):
        if event.type() == Qt.QEvent.KeyPress and event.key() == Qt.Qt.Key_D:
            # Kept around mostly for debugging.
            self.setup_next_spline()
        return False
