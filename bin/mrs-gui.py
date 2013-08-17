#! /usr/bin/env ipython

"""
GUI for the analysis of data from MRS experiments
-------------------------------------------------

This GUI (graphical user interface) will do the following:

- Load data from a p-file
- Display spectra (on, off and diff), while allowing dynamic changing of:
   - Line-widening 
   - phase correction (zero order and first order)

- Output spectra to a csv file

"""
import sys, os, csv

import numpy as np
## Using PySide
from PySide import QtCore
from PySide import QtGui

import matplotlib

## Added for PySide
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

import MRS.analysis as ana
import MRS.utils as ut
import nitime as nt

class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.setWindowTitle('MRS')
        self.data = DataHolder()
        self.series_list_model = QtGui.QStandardItemModel()
        
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        
        self.load_file()
        self.update_ui()
        self.on_show()

    def create_main_frame(self):
        self.main_frame = QtGui.QWidget()
        plot_frame = QtGui.QWidget()
        
        self.dpi = 100
        self.fig = Figure((6.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        self.axes = self.fig.add_subplot(111)
        self.axes.invert_xaxis()
        self.axes.set_xlabel('ppm')
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        log_label = QtGui.QLabel("Data series:")
        self.series_list_view = QtGui.QListView()
        self.series_list_view.setModel(self.series_list_model)
        
        spin_label1 = QtGui.QLabel('Phase:')
        self.phase_spin = QtGui.QDoubleSpinBox(singleStep=0.1)
        spin_label2 = QtGui.QLabel('Linewidth:')
        self.line_spin = QtGui.QDoubleSpinBox(singleStep=0.1)
        spin_label3 = QtGui.QLabel('Cutoff:')
        self.cutoff_spin = QtGui.QDoubleSpinBox(singleStep=0.5)
        
        spins_hbox = QtGui.QHBoxLayout()
        spins_hbox.addWidget(spin_label1)
        spins_hbox.addWidget(self.phase_spin)
        spins_hbox.addWidget(spin_label2)
        spins_hbox.addWidget(self.line_spin)
        spins_hbox.addWidget(spin_label3)
        spins_hbox.addWidget(self.cutoff_spin)
        spins_hbox.addStretch(1)

        ticks_hbox = QtGui.QHBoxLayout()
        self.legend_cb = QtGui.QCheckBox("Show L&egend")
        ticks_hbox.addWidget(self.legend_cb) 

        buttons_hbox = QtGui.QHBoxLayout()
        self.show_button = QtGui.QPushButton("&Show")
        self.connect(self.show_button, QtCore.SIGNAL('clicked()'), self.on_show)
        buttons_hbox.addWidget(self.show_button)

        self.fit_button = QtGui.QPushButton("&Fit")
        self.connect(self.fit_button, QtCore.SIGNAL('clicked()'), self.on_fit)
        buttons_hbox.addWidget(self.fit_button)

        self.auc_button = QtGui.QPushButton("&AUC")
        self.connect(self.auc_button, QtCore.SIGNAL('clicked()'), self.on_auc)
        buttons_hbox.addWidget(self.auc_button)
        
        left_vbox = QtGui.QVBoxLayout()
        left_vbox.addWidget(self.canvas)
        left_vbox.addWidget(self.mpl_toolbar)

        right_vbox = QtGui.QVBoxLayout()
        right_vbox.addWidget(log_label)
        right_vbox.addWidget(self.series_list_view)
        right_vbox.addLayout(spins_hbox)
        right_vbox.addLayout(ticks_hbox)
        right_vbox.addLayout(buttons_hbox)
        right_vbox.addStretch(1)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addLayout(left_vbox)
        hbox.addLayout(right_vbox)
        self.main_frame.setLayout(hbox)

        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QtGui.QLabel("Please load a data file")
        self.statusBar().addWidget(self.status_text, 1)

    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_action = self.create_action("&Load file",
            shortcut="Ctrl+L", slot=self.load_file, tip="Load a file")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_action, None, quit_action))
            
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
        self.add_actions(self.help_menu, (about_action,))
        
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QtGui.QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, QtCore.SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


    def load_file(self, filename=None):
        filename = QtGui.QFileDialog.getOpenFileName(self,
            'Open a data file', '.', 'P files (*.7)')
        
        if filename:
            self.data.load_from_file(filename)
            self.data.get_spectra(self.line_spin.value(),
                                  self.cutoff_spin.value())
            self.fill_series_list(self.data.series_names())
            self.status_text.setText("Loaded " + filename[0])
            self.update_ui()
    
    def update_ui(self):
        if self.data.series_count() > 0 and self.data.series_len() > 0:
            self.phase_spin.setValue(0)
            self.phase_spin.setRange(-np.pi, np.pi)
            self.phase_spin.setEnabled(True)
            
            self.line_spin.setRange(0, 100)
            self.line_spin.setEnabled(True)
            self.line_spin.setValue(5.0)

            self.cutoff_spin.setRange(0, 1000)
            self.cutoff_spin.setEnabled(True)
            self.cutoff_spin.setValue(1)

        else:
            for w in [self.phase_spin, self.first_spin]:
                w.setEnabled(False)
    
    def on_show(self):
        self.axes.clear()
        self.axes.set_xlabel('ppm')
        self.axes.grid(True)
        has_series = False
        self.data.get_spectra(self.line_spin.value(),
                              self.cutoff_spin.value())
        
        for row in range(self.series_list_model.rowCount()):
            model_index = self.series_list_model.index(row, 0)
            checked = self.series_list_model.data(model_index,
                QtCore.Qt.CheckStateRole) == QtCore.Qt.Checked 
            name = str(self.series_list_model.data(model_index))
            
            if checked:
                has_series = True
                phase_zero = self.phase_spin.value()
                series = self.data.get_series_data(name)
                series = ut.phase_correct_zero(series, phase_zero)
                self.axes.plot(self.data.f_ppm, series, label=name)
        
        if has_series and self.legend_cb.isChecked():
            self.axes.legend()
        self.canvas.draw()

    def on_fit(self):
        """

        """
        print "This is where fitting data will go"
        # XXX Need to continue here:

        # Need to update the series list in the GUI to show the fit: 
        ## self.fill_series_list(self.data.series_names())
        ## self.update_ui()


    def on_auc(self):
        """

        """
        xlim = self.axes.get_xlim()
        line_list = self.axes.get_lines()
        txt = ''
        for row in range(self.series_list_model.rowCount()):
            model_index = self.series_list_model.index(row, 0)
            checked = self.series_list_model.data(model_index,
                QtCore.Qt.CheckStateRole) == QtCore.Qt.Checked 
            name = str(self.series_list_model.data(model_index))

            if checked:
                xy = line_list[row].get_xydata()
                idx = np.where(np.logical_and(xy[:,0]>=np.min(xlim),
                                              xy[:,0]<=np.max(xlim)))
                this_y = xy[idx, 1].squeeze()
                # Use both sides to find the baseline:
                baseline = np.mean([this_y[0], this_y[1]])
                this_y -= baseline
                # Integrate as sum * d(f_ppm)
                txt = txt + "%s AUC: %2.2f \n"%(name, np.sum(this_y))

        self.axes.text(0.05, 0.95, txt, horizontalalignment='left',
                       verticalalignment='top',
                       transform = self.axes.transAxes)
        self.canvas.draw()

    def on_about(self):
        msg = __doc__
        QMessageBox.about(self, "About the demo", msg.strip())

    def fill_series_list(self, names):
        self.series_list_model.clear()
        
        for name in names:
            item = QtGui.QStandardItem(name)
            item.setCheckState(QtCore.Qt.Unchecked)
            item.setCheckable(True)
            self.series_list_model.appendRow(item)
            
        


class DataHolder(object):
    """ Just a thin wrapper over a dictionary that holds integer 
        data series. Each series has a name and a list of numbers 
        as its data. The length of all series is assumed to be
        the same.
        
        The series can be read from a CSV file, where each line
        is a separate series. In each series, the first item in 
        the line is the name, and the rest are data numbers.
    """
    def __init__(self, filename=None):
        self.load_from_file(filename)
    
    def load_from_file(self, filename=None):
        if filename:
            # Get data from file: 
            data = np.transpose(nib.load(filename[0]).get_data(),
                                [1,2,3,4,5,0]).squeeze()

            # Use the water unsuppressed data to combine over coils:
            w_data, w_supp_data = ana.coil_combine(data.squeeze())
            self.timeseries = w_supp_data

    def series_names(self):
        """ Names of the data series
        """
        return self.names

        
    def series_count(self):
        return len(self.data)


    def series_len(self):
        """ Length of a data series
        """
        return self.datalen


    def get_series_data(self, name):
        return self.data[name]

        
    def get_spectra(self, lbr, cutoff):
        self.data = {}
        self.names = []
        d, ppm = self.calc_spectra(lbr, cutoff)
        for k,v in d.items():
            self.names.append(k)
            self.data[k] =  v
            self.datalen = len(v)
        self.ppm = ppm

    def calc_spectra(self, lbr, cutoff,
                     sampling_rate=5000, min_ppm=-0.7, max_ppm=4.3):
        # Once we've done that, we only care about the water-suppressed data
        self.f, self.spec = ana.get_spectra(nt.TimeSeries(self.timeseries,
                                            sampling_rate=sampling_rate),
                                            line_broadening=lbr,
                            filt_method = dict(lb=cutoff, filt_order=256))

        # The first echo (off-resonance) is in the first output 
        self.echo1 = self.spec[:,0]
        # The on-resonance is in the second:
        self.echo2 = self.spec[:,1]
        f_ppm = ut.freq_to_ppm(self.f)
        idx0 = np.argmin(np.abs(f_ppm - min_ppm))
        idx1 = np.argmin(np.abs(f_ppm - max_ppm))
        idx = slice(idx1, idx0)
        # Convert from Hz to ppm and extract the part you are interested in.
        self.f_ppm = f_ppm[idx]

        # Pack it into a dict:
        m_e1 = np.mean(self.echo1[:,idx], 0)
        m_e2 = np.mean(self.echo2[:,idx], 0)
        self.diff = m_e2 - m_e1
        
        return dict(echo1=m_e1, echo2=m_e2, diff=self.diff), f_ppm

if __name__ == "__main__":
    # In interactive work, there might already be an app there: 
    app=QtGui.QApplication.instance()
    # create QApplication if it doesnt exist 
    if not app: 
        app = QtGui.QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    app.exec_()
        
    
    
     
