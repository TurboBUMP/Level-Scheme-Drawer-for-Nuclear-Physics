import sys
import os

try:
    from PyQt5.QtWidgets import (
        QApplication, QMainWindow, QWidget, QLabel, QLineEdit,
        QPushButton, QCheckBox, QGridLayout, QMessageBox,
        QSizePolicy, QVBoxLayout, QMenu, QAction
    )
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont
except ImportError:
    print("Unable to import PyQt5, please install it before proceeding.")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Unable to import Numpy, please install it before proceeding.")
    sys.exit(1)

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except ImportError:
    print("Unable to import Matplotlib, please install it before proceeding.")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("Unable to import Pandas, please install it before proceeding.")
    sys.exit(1)

try:
    from lvlDrawFunc import drawLevelScheme
except ImportError:
    print("Error occurred while importing lvlDrawFunc\n")
    sys.exit(1)


FONT_FAMILY = "Noto Sans"
FONT_SIZE   = 13
WIDTH_CONST = 0.05   # fraction of x_max used for label margins

# Folder scanned by the right-click context menu on the file entries.
# Resolved relative to the directory that contains draw.py.
_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FILES_FOLDER = os.path.join(_SCRIPT_DIR, "..", "files")


def make_font(size=None):
    return QFont(FONT_FAMILY, size or FONT_SIZE)


# ─────────────────────────────────────────────────────────────────────────────
#  Label-collision solver
# ─────────────────────────────────────────────────────────────────────────────

def resolve_label_positions(energies, min_gap, draw_gs=True, spread=False):
    """
    Place label y-positions in one of two modes, selected by the
    "Spread labels evenly" checkbox.

    spread=False (default, box unchecked):
      Cluster-based bidirectional algorithm. Consecutive levels closer than
      min_gap form a cluster; each cluster's labels are spread symmetrically
      around the cluster mean. Isolated levels stay exactly at their energy.
      Only the minimum displacement needed to avoid overlap is applied.

    spread=True (box checked):
      Simple linspace from the GS level (y=0 when draw_gs=True, else e_min)
      to the top level energy. No clustering — every label gets an equal share
      of the full vertical range regardless of where the levels actually are.
    """
    n = len(energies)
    if n == 0:
        return np.array([])

    e = np.array(energies, dtype=float)

    if n == 1:
        return e.copy()

    # ── spread=True: evenly spaced from GS to top level ──────────────────────
    if spread:
        y_min = 0.0 + min_gap if draw_gs else e[0]
        y_max = e[-1]
        return np.linspace(y_min, y_max, n)

    # ── spread=False: cluster-based bidirectional spreading ───────────────────
    gs_floor = 0.0 if draw_gs else -np.inf

    # Step 1: build clusters
    clusters = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and (e[j + 1] - e[j]) < min_gap:
            j += 1
        clusters.append((i, j))
        i = j + 1

    # Step 2: spread each cluster symmetrically around its mean energy
    pos = np.empty(n, dtype=float)
    for (ci, cj) in clusters:
        k       = cj - ci + 1
        offsets = (np.arange(k) - (k - 1) / 2.0) * min_gap
        pos[ci:cj + 1] = e[ci:cj + 1].mean() + offsets

    # Step 3: clamp bottom against the GS floor
    if pos[0] < gs_floor:
        pos[clusters[0][0]:clusters[0][1] + 1] += gs_floor - pos[0]

    # Step 4: upward sweep to fix residual inter-cluster overlaps
    for i in range(1, n):
        needed = pos[i - 1] + min_gap
        if pos[i] < needed:
            owner = next((ci, cj) for (ci, cj) in clusters if ci <= i <= cj)
            pos[owner[0]:owner[1] + 1] += needed - pos[i]

    return pos


# ─────────────────────────────────────────────────────────────────────────────
#  Drawing window  (resizable, canvas fills all space)
# ─────────────────────────────────────────────────────────────────────────────

class DrawingWindow(QMainWindow):
    """Separate window that displays the matplotlib level scheme."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Level Scheme")
        self.resize(960, 860)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Figure with tight layout so axes fill the canvas
        self.fig = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()

        toolbar = NavigationToolbar(self.canvas, self)

        layout.addWidget(toolbar)
        layout.addWidget(self.canvas)

    def render(self, params: dict):
        """
        Re-draw the level scheme.  `params` is the dict built by SetParameters.
        Called both on first show and on every resize.
        """
        self.fig.clear()
        ax = self.fig.add_subplot(111)

        # Derive x_max from canvas width in inches so arrows/labels
        # are proportional to the actual rendered width.
        dpi    = self.fig.get_dpi()
        w_inch = self.canvas.width()  / dpi
        h_inch = self.canvas.height() / dpi

        n_tr   = params['stop_tr'] - params['start_tr']
        x_max  = max(w_inch * 0.8, 1.0)   # use 80% of figure width in data units
        delta_x = x_max / n_tr if n_tr > 0 else x_max

        wc = WIDTH_CONST
        x_fig_start       = -x_max * wc
        x_fig_end         =  x_max + x_max * wc
        x_left_label_dist =  x_fig_start
        x_right_label_dist = x_fig_end

        # Scale arrow geometry with x_max so it looks right at any size
        scale          = x_max / 5.0          # 5.0 is the historical default
        arrow_width    = params['arrow_width']    * scale
        arrow_hw       = params['arrow_head_width'] * scale
        arrow_hl       = params['arrow_head_length']   # keep in data-y units

        drawLevelScheme(
            params['nucleus_name'],
            params['nucleus_name_fontsize'],
            self.fig, ax,
            params['levels_pandas'],
            params['transitions_pandas'],
            delta_x, x_max,
            x_fig_start, x_fig_end,
            x_right_label_dist, x_left_label_dist,
            params['fontsize'],
            params['start_lvl'], params['stop_lvl'],
            params['start_tr'],  params['stop_tr'],
            arrow_width, arrow_hw, arrow_hl,
            params['arrow_color'],
            params['draw_gs'],
            params['draw_all_aligned'],
            params['energy_label_rotation'],
        )

        self.canvas.draw()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Re-render when the window is resized, if we have params
        if hasattr(self, '_params') and self._params is not None:
            self.render(self._params)

    def show_scheme(self, params: dict):
        self._params = params
        self.render(params)
        self.show()


# ─────────────────────────────────────────────────────────────────────────────
#  Parameter window
# ─────────────────────────────────────────────────────────────────────────────

class SetParameters(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Nuclear Physics Level Scheme")
        self.resize(1200, 740)

        self.uploaded_flag      = False
        self.levels_pandas      = None
        self.transitions_pandas = None
        self._drawing_win       = None

        self.defaults = {
            "transition_file":       "../files/transitions.csv",
            "level_file":            "../files/levels.csv",
            "arrow_width":           0.001,
            "arrow_head_width":      0.005,
            "arrow_head_length":     40.0,
            "energy_label_rotation": 90.0,
            "min_vert_label_dist":   100.0,
            "fontsize":              10,
            "arrow_color":           "black",
            "start_level":           0,
            "stop_level":            -1,
            "start_transition":      0,
            "stop_transition":       -1,
            "nucleus_name":          "42Ca",
            "nucleus_name_fontsize": 20,
        }

        self._build_ui()

    # ── UI ───────────────────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        grid = QGridLayout(central)
        grid.setContentsMargins(10, 10, 10, 10)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(6)

        def lbl(text):
            w = QLabel(text)
            w.setFont(make_font())
            return w

        def entry(key):
            w = QLineEdit(str(self.defaults[key]))
            w.setFont(make_font())
            w.setMinimumWidth(160)
            return w

        def chk(text, checked=False):
            w = QCheckBox(text)
            w.setFont(make_font())
            w.setChecked(checked)
            return w

        # widgets
        self.lbl_trans_file  = lbl("Transition Filename:")
        self.lbl_level_file  = lbl("Level Filename:")
        self.ent_trans_file  = entry("transition_file")
        self.ent_level_file  = entry("level_file")

        self.lbl_arrow_w     = lbl("Arrow Width:")
        self.lbl_arrow_hw    = lbl("Arrow Head Width:")
        self.lbl_arrow_hl    = lbl("Arrow Head Length:")
        self.lbl_elabel_rot  = lbl("Energy Label Rotation:")
        self.lbl_min_vdist   = lbl("Minimum Vertical Label Separation:")
        self.lbl_fontsize    = lbl("Fontsize:")
        self.lbl_arrow_col   = lbl("Arrow Color:")
        self.lbl_start_lvl   = lbl("Start Level:")
        self.lbl_stop_lvl    = lbl("Stop Level:")
        self.lbl_start_tr    = lbl("Start Transition:")
        self.lbl_stop_tr     = lbl("Stop Transition:")
        self.lbl_nuc_name    = lbl("Nucleus Name:")
        self.lbl_nuc_fs      = lbl("Name Fontsize:")

        self.ent_arrow_w     = entry("arrow_width")
        self.ent_arrow_hw    = entry("arrow_head_width")
        self.ent_arrow_hl    = entry("arrow_head_length")
        self.ent_elabel_rot  = entry("energy_label_rotation")
        self.ent_min_vdist   = entry("min_vert_label_dist")
        self.ent_fontsize    = entry("fontsize")
        self.ent_arrow_col   = entry("arrow_color")
        self.ent_start_lvl   = entry("start_level")
        self.ent_stop_lvl    = entry("stop_level")
        self.ent_start_tr    = entry("start_transition")
        self.ent_stop_tr     = entry("stop_transition")
        self.ent_nuc_name    = entry("nucleus_name")
        self.ent_nuc_fs      = entry("nucleus_name_fontsize")

        self.btn_upload      = QPushButton("Upload Files / Update Values")
        self.btn_upload.setFont(make_font())
        self.btn_draw        = QPushButton("Draw")
        self.btn_draw.setFont(make_font())

        self.chk_draw_gs        = chk("Draw G.S.", checked=True)
        self.chk_all_aligned    = chk("Draw Transitions Vertically Aligned")
        self.chk_spread_labels  = chk("Spread labels evenly")

        # Dropdown arrow buttons next to file entries
        self.btn_trans_browse = QPushButton("▼")
        self.btn_trans_browse.setFont(make_font())
        self.btn_trans_browse.setFixedWidth(32)
        self.btn_trans_browse.setToolTip("Pick a file from ../files/")
        self.btn_level_browse = QPushButton("▼")
        self.btn_level_browse.setFont(make_font())
        self.btn_level_browse.setFixedWidth(32)
        self.btn_level_browse.setToolTip("Pick a file from ../files/")

        # grid layout
        r = 0
        grid.addWidget(self.lbl_trans_file,  r, 0); grid.addWidget(self.ent_trans_file, r, 1)
        grid.addWidget(self.btn_trans_browse, r, 2)
        grid.addWidget(self.btn_upload,       r, 3)
        r += 1
        grid.addWidget(self.lbl_level_file,  r, 0); grid.addWidget(self.ent_level_file, r, 1)
        grid.addWidget(self.btn_level_browse, r, 2)
        r += 2
        grid.addWidget(self.lbl_arrow_w,     r, 0); grid.addWidget(self.ent_arrow_w,    r, 1)
        r += 1
        grid.addWidget(self.lbl_arrow_hw,    r, 0); grid.addWidget(self.ent_arrow_hw,   r, 1)
        r += 1
        grid.addWidget(self.lbl_arrow_hl,    r, 0); grid.addWidget(self.ent_arrow_hl,   r, 1)
        grid.addWidget(self.btn_draw,        r, 2, 1, 2)
        r += 1
        grid.addWidget(self.lbl_elabel_rot,  r, 0); grid.addWidget(self.ent_elabel_rot, r, 1)
        r += 1
        grid.addWidget(self.lbl_min_vdist,   r, 0); grid.addWidget(self.ent_min_vdist,  r, 1)
        r += 1
        grid.addWidget(self.lbl_fontsize,    r, 0); grid.addWidget(self.ent_fontsize,   r, 1)
        r += 1
        grid.addWidget(self.lbl_arrow_col,   r, 0); grid.addWidget(self.ent_arrow_col,  r, 1)
        r += 1
        grid.addWidget(self.lbl_start_lvl,   r, 0); grid.addWidget(self.ent_start_lvl,  r, 1)
        r += 1
        grid.addWidget(self.lbl_stop_lvl,    r, 0); grid.addWidget(self.ent_stop_lvl,   r, 1)
        r += 1
        grid.addWidget(self.lbl_start_tr,    r, 0); grid.addWidget(self.ent_start_tr,   r, 1)
        r += 1
        grid.addWidget(self.lbl_stop_tr,     r, 0); grid.addWidget(self.ent_stop_tr,    r, 1)
        r += 1
        grid.addWidget(self.lbl_nuc_name,    r, 0); grid.addWidget(self.ent_nuc_name,   r, 1)
        grid.addWidget(self.lbl_nuc_fs,      r, 2); grid.addWidget(self.ent_nuc_fs,     r, 3)
        r += 1
        grid.addWidget(self.chk_draw_gs,       r, 0)
        grid.addWidget(self.chk_all_aligned,   r, 1, 1, 2)
        r += 1
        grid.addWidget(self.chk_spread_labels, r, 0)

        self.btn_upload.clicked.connect(self.upload_files)
        self.btn_draw.clicked.connect(self.draw)
        self.btn_trans_browse.clicked.connect(
            lambda: self._show_file_menu(self.ent_trans_file,
                                         self.btn_trans_browse))
        self.btn_level_browse.clicked.connect(
            lambda: self._show_file_menu(self.ent_level_file,
                                         self.btn_level_browse))


    # ── helpers ───────────────────────────────────────────────────────────────

    def _show_file_menu(self, entry_widget, anchor_widget):
        """
        Dropdown menu anchored below anchor_widget (the ▼ button).
        Lists every .csv file found in FILES_FOLDER; clicking one writes
        its relative path (../files/<name>) into entry_widget.
        """
        folder = os.path.normpath(FILES_FOLDER)
        try:
            csv_files = sorted(
                f for f in os.listdir(folder)
                if f.lower().endswith('.csv')
            )
        except FileNotFoundError:
            csv_files = []

        menu = QMenu(self)
        menu.setFont(make_font())

        if csv_files:
            for name in csv_files:
                # Store as a ../files/-relative path so it matches the
                # default convention already used throughout the code.
                rel_path = os.path.join("..", "files", name)
                action = QAction(name, self)
                action.triggered.connect(
                    lambda checked, p=rel_path: entry_widget.setText(p))
                menu.addAction(action)
        else:
            no_files = QAction("(no .csv files found)", self)
            no_files.setEnabled(False)
            menu.addAction(no_files)

        # Pop up just below the bottom-left corner of the arrow button
        btn_rect = anchor_widget.rect()
        popup_pos = anchor_widget.mapToGlobal(btn_rect.bottomLeft())
        menu.exec_(popup_pos)

    def _set_file_style(self, ok):
        if ok is None:
            for w in (self.ent_trans_file, self.ent_level_file):
                w.setStyleSheet("")
            for w in (self.lbl_trans_file, self.lbl_level_file):
                w.setStyleSheet("")
            return
        colour = "green" if ok else "red"
        for w in (self.ent_trans_file, self.ent_level_file):
            w.setStyleSheet(f"background-color: {colour}; color: white;")
        for w in (self.lbl_trans_file, self.lbl_level_file):
            w.setStyleSheet(f"color: {colour};")

    def _float(self, w, name):
        try:
            return float(w.text())
        except ValueError:
            raise ValueError(f"Invalid value for '{name}': {w.text()!r}")

    def _int(self, w, name):
        try:
            return int(w.text())
        except ValueError:
            raise ValueError(f"Invalid value for '{name}': {w.text()!r}")

    def _error(self, msg):
        QMessageBox.critical(self, "Error", msg)

    # ── upload ────────────────────────────────────────────────────────────────

    def upload_files(self):
        try:
            self.uploaded_flag = False
            self._set_file_style(None)

            min_vdist   = self._float(self.ent_min_vdist, "Min Vert Label Dist")
            draw_gs        = self.chk_draw_gs.isChecked()
            spread_labels  = self.chk_spread_labels.isChecked()

            lp = pd.read_csv(
                self.ent_level_file.text(),
                names=['Level energy', 'Spin-Parity',
                       'Energy Label Position', 'Level color'],
                dtype={'Level energy': float, 'Spin-Parity': str,
                       'Energy Label Position': str, 'Level color': str},
                comment='#'
            )
            lp = lp.sort_values(by=['Level energy']).reset_index(drop=True)

            tp = pd.read_csv(
                self.ent_trans_file.text(),
                names=["Transition energy", "Initial level", "Final level",
                       "Multipolarity", "Transition color"],
                dtype={"Transition energy": float, "Initial level": float,
                       "Final level": float, "Multipolarity": str,
                       "Transition color": str},
                comment='#'
            )
            tp = tp.sort_values(
                by=["Initial level", "Transition energy"],
                ascending=True).reset_index(drop=True)

            # Resolve stop sentinels
            stop_level = lp.shape[0]
            stop_trans = tp.shape[0]
            self.ent_stop_lvl.setText(str(stop_level))
            self.ent_stop_tr.setText(str(stop_trans))
            self.ent_start_lvl.setText("0")
            self.ent_start_tr.setText("0")

            # ── improved label placement ──────────────────────────────────
            energies = lp['Level energy'].values.astype(float)
            label_pos = resolve_label_positions(
                energies, min_vdist,
                draw_gs=draw_gs, spread=spread_labels)
            lp['Energy Label Position'] = label_pos
            # ─────────────────────────────────────────────────────────────

            lp.fillna('', inplace=True)

            self.levels_pandas      = lp
            self.transitions_pandas = tp

            self.uploaded_flag = True
            self._set_file_style(True)
            print("Files uploaded successfully!")

        except Exception as e:
            self.uploaded_flag = False
            self._set_file_style(False)
            self._error(f"Something went wrong while uploading:\n{e}")
            print(f"Upload error: {e}")

    # ── draw ──────────────────────────────────────────────────────────────────

    def draw(self):
        if not self.uploaded_flag:
            self._error("No file to draw. Please upload valid files first!")
            return
        try:
            params = dict(
                levels_pandas        = self.levels_pandas,
                transitions_pandas   = self.transitions_pandas,
                nucleus_name         = self.ent_nuc_name.text(),
                nucleus_name_fontsize= self._int(self.ent_nuc_fs,      "Name Fontsize"),
                fontsize             = self._int(self.ent_fontsize,     "Fontsize"),
                start_lvl            = self._int(self.ent_start_lvl,   "Start Level"),
                stop_lvl             = self._int(self.ent_stop_lvl,    "Stop Level"),
                start_tr             = self._int(self.ent_start_tr,    "Start Transition"),
                stop_tr              = self._int(self.ent_stop_tr,     "Stop Transition"),
                arrow_width          = self._float(self.ent_arrow_w,   "Arrow Width"),
                arrow_head_width     = self._float(self.ent_arrow_hw,  "Arrow Head Width"),
                arrow_head_length    = self._float(self.ent_arrow_hl,  "Arrow Head Length"),
                arrow_color          = self.ent_arrow_col.text(),
                draw_gs              = int(self.chk_draw_gs.isChecked()),
                draw_all_aligned     = int(self.chk_all_aligned.isChecked()),
                energy_label_rotation= self._float(self.ent_elabel_rot,"Energy Label Rotation"),
            )

            if self._drawing_win is None:
                self._drawing_win = DrawingWindow()
            self._drawing_win.show_scheme(params)

        except Exception as e:
            self._error(f"Something went wrong while drawing:\n{e}")
            print(f"Draw error: {e}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SetParameters()
    window.show()
    sys.exit(app.exec_())
