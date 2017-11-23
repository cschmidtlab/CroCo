#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Windows for the CroCo cross-link converter GUI programme
"""

from PyQt5.QtWidgets import QWidget, QToolTip, QPushButton, QApplication,\
                            QMessageBox, QDesktopWidget, QGridLayout,\
                            QComboBox, QLabel, QFileDialog, QAction,\
                            QMainWindow, QHBoxLayout, QTabWidget, QLineEdit
from PyQt5.QtGui import QFont, QIcon

import PyQt5.QtCore as QtCore

import os
os.getcwd()

import sys
import CroCo_reader as xr
import CroCo_writer as xw

import pandas as pd

class CroCo_MainWindow(QMainWindow):

    # initialise properties
    input_format = 'pLink'
    output_format = 'xTable'
    output_dir = ''
    fnames = []

    def __init__(self):
        # initialise the parent classe
        super().__init__()
        # initial itself
        self.createComponents()
        self.createLayout()
        self.createMenu()

    def createComponents(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """

        self.setWindowTitle('The CroCo Crosslink Converter')
        self.setWindowIcon(QIcon('images/python-icon.png')) # icon in ubuntu is
                                                       # displayed on sidepanel
        QToolTip.setFont(QFont('SansSerif', 10))

        #############################################
        # Crosslink Converter Part
        #############################################

        # Define Dropdown Menu for Input
        self.input_lbl = QLabel('Input file format', self)
        self.input_dropdown = QComboBox(self)
        self.input_dropdown.addItem("pLink")
        self.input_dropdown.addItem("Kojak")
        self.input_dropdown.addItem("xTable")
        self.input_dropdown.activated[str].connect(self.set_input_format)

        # Dropdown for output
        self.output_lbl = QLabel('Output file format', self)
        self.output_dropdown = QComboBox(self)
        self.output_dropdown.addItem("xTable")
        self.output_dropdown.addItem("xVis")
        self.output_dropdown.addItem("xiNet")
        self.output_dropdown.activated[str].connect(self.set_output_format)

        # Start button
        self.sbtn = QPushButton('Start', self)
        self.sbtn.setToolTip('Start the conversion')
        self.sbtn.clicked.connect(self.start_conversion)
        self.sbtn.resize(self.sbtn.sizeHint())

        # Close button
        self.qbtn = QPushButton('Quit', self)
        # call a close Event on clicking the button
        self.qbtn.clicked.connect(self.close)
        self.qbtn.resize(self.qbtn.sizeHint())

        # Open file dialog
        self.fbtn = QPushButton('Load file(s)', self)
        self.fbtn.clicked.connect(self.set_input_files)

        # Output files button
        self.obtn = QPushButton('Output dir', self)
        self.obtn.clicked.connect(self.set_output_dir)

        # Label for inpath
        self.inpath_title = QLabel('Reading from:')
        self.inpath_lbl = QLabel('')

        # Label for outpath
        self.outpath_title = QLabel('Writing to:')
        self.outpath_lbl = QLabel(self.output_dir)

        #############################################
        # Spectrum annotation
        #############################################

        # Open xlink file dialog
        self.assign_lxlink_lbl = QLabel('Load xlinks:', self)
        self.assign_lxlink = QPushButton('Load xlink file', self)

        # Open mgf file dialog
        self.assign_lmgf_lbl = QLabel('Load mgf', self)
        self.assign_lmgf = QPushButton('Load file', self)

        # Output files button
        self.assign_obtn_lbl = QLabel('Output dir', self)
        self.assign_obtn = QPushButton('Set output dir', self)

        # Start button
        self.assign_start = QPushButton('Start assignment', self)

        # PSM annotation
        self.assign_accept = QPushButton('Accept', self)
        self.assign_unsure = QPushButton('Unsure', self)
        self.assign_decline = QPushButton('Decline', self)

        # Navigation
        self.assign_spectrum_lbl = QLabel('Current Spectrum:', self)
        self.assign_xlink_lbl = QLabel('Current xlink:', self)
        self.assign_xlink_input = QLineEdit()
        self.assign_goto = QPushButton('Go to', self)

        self.assign_next = QPushButton('Next', self)
        self.assign_prev = QPushButton('Previous', self)


    #####################################
    # Set Variables on self
    #####################################

        # assign to variables
    def set_input_format(self, text):
        self.input_format = text

    def set_output_format(self, text):
        self.output_format = text

        # open folder dialog for pLink (only one folder at once) or file
        # dialog for others (multiple files at once)
        # decorator explicitely defines method as slot!
    @QtCore.pyqtSlot()
    def set_input_files(self):
        if self.input_format == 'pLink':
            self.fnames = [QFileDialog.getExistingDirectory(self,\
                                                'Select directory')]
        else:
            self.fnames = QFileDialog.getOpenFileNames(self, 'Open file',\
                                                       '/home')[:-1][0]
        # update input label
        self.inpath_lbl.setText(os.path.basename(self.fnames[0]))
    @QtCore.pyqtSlot()
    def set_output_dir(self):
        self.output_dir = QFileDialog.getExistingDirectory(self,\
                                                'Select directory for output')
        # update output label
        self.outpath_lbl.setText(os.path.basename(self.output_dir))

    def createMenu(self):
        """
        Creates the application Menu
        """
        self.actionClose = QAction(self.tr('Quit'), self)
        # Only relevant for macOS
        self.actionClose.setMenuRole(QAction.QuitRole)
        self.actionClose.setShortcut('Ctrl+Q')
        self.actionClose.setStatusTip('Exit application')
        self.actionClose.triggered.connect(self.close)


        self.actionAbout = QAction(self.tr('About'), self)
        self.actionAbout.setMenuRole(QAction.AboutRole) # only for macOS
        self.actionClose.setStatusTip('About CroCo')
        self.actionAbout.triggered.connect(self.show_about)


        menuFile = self.menuBar().addMenu(self.tr('File'))
        menuFile.addAction(self.actionClose)


        menuAbout = self.menuBar().addMenu(self.tr('About'))
        menuAbout.addAction(self.actionAbout)


    def createLayout(self):
        """
        Takes the logical in- and outputs and links them to optical
        widgets such as buttons
        """

        ###############################################
        # Define Layout for converter page as grid
        ###############################################

        converter_layout = QGridLayout()
        converter_layout.setSpacing(10)
        converter_layout.setRowStretch(0, 1)
        converter_layout.setRowStretch(1,1)
        converter_layout.setRowStretch(2,1)
        converter_layout.setRowStretch(3,2)
        converter_layout.setRowStretch(4,1)

        # add the 1st line of elements to the grid
        converter_layout.addWidget(self.input_lbl, 0, 0)
        converter_layout.addWidget(self.input_dropdown, 0, 1)
        converter_layout.addWidget(self.output_lbl, 0, 2)
        converter_layout.addWidget(self.output_dropdown, 0, 3)
        # 2nd line
        converter_layout.addWidget(self.fbtn, 1, 1)
        converter_layout.addWidget(self.obtn, 1, 3)
        #3rd and 4th line
        converter_layout.addWidget(self.inpath_title, 2, 0)
        converter_layout.addWidget(self.inpath_lbl, 3, 0)
        converter_layout.addWidget(self.outpath_title, 2, 2)
        converter_layout.addWidget(self.outpath_lbl, 3, 2)
        # 5th line
        converter_layout.addWidget(self.sbtn, 4, 0)
        converter_layout.addWidget(self.qbtn, 4, 3)

        # generate widget for converter window
        converter_widget = QWidget()
        # apply layout
        converter_widget.setLayout(converter_layout)

        ###############################################
        # Define Layout for spectrum assignment page
        ###############################################

        assign_layout = QGridLayout()

        assign_layout.addWidget(self.assign_lxlink_lbl, 0, 0)
        assign_layout.addWidget(self.assign_lmgf_lbl, 0, 2)
        assign_layout.addWidget(self.assign_obtn_lbl, 0, 4)

        assign_layout.addWidget(self.assign_lxlink, 1, 0)
        assign_layout.addWidget(self.assign_lmgf, 1, 2)
        assign_layout.addWidget(self.assign_obtn, 1, 4)

        assign_layout.addWidget(self.assign_start, 2, 2)
        assign_layout.addWidget(self.assign_xlink_lbl, 2, 0)
        assign_layout.addWidget(self.assign_spectrum_lbl, 2, 3)

        assign_layout.addWidget(self.assign_accept, 3, 2)

        assign_layout.addWidget(self.assign_xlink_input, 4, 0)
        assign_layout.addWidget(self.assign_goto, 4, 1)
        assign_layout.addWidget(self.assign_unsure, 4, 2)

        assign_layout.addWidget(self.assign_prev, 5, 0)
        assign_layout.addWidget(self.assign_decline, 5, 2)
        assign_layout.addWidget(self.assign_next, 5, 4)
        # generate widget for assignment window
        assign_widget = QWidget()
        # apply layout
        assign_widget.setLayout(assign_layout)

        #############################################
        # Tabbing
        #############################################

        self.tabWidget = QTabWidget()
        self.tabWidget.addTab(converter_widget, 'Conversion')
        self.tabWidget.addTab(assign_widget, 'Assignment')

        ##############################
        # Define main layout
        ##############################
        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.tabWidget)
        # generate widget for main window
        mainWidget = QWidget()
        # apply layout
        mainWidget.setLayout(mainLayout)
        # assign widget to main window (i.e. self)
        self.setCentralWidget(mainWidget)

        # set window size
        self.resize(250, 200)
        # center
        self.center()
        # show the window
        self.show()

    @QtCore.pyqtSlot()
    def start_conversion(self):
        """
        Collect all necessary information from self and start the
        conversion by calling the actual conversion script
        """

        print('Going to convert {} from {} '.format(', '.join(self.fnames),
                                                   self.input_format) +
               'format to {} format'.format(self.output_format))

        in_dict = {'pLink': xr.ReadpLink,
                   'Kojak': xr.ReadKojak,
                   'xTable': pd.read_csv}

        out_dict = {'xTable': xw.WriteXtable,
                    'xVis': xw.WritexVis,
                    'xiNet': xw.WritexiNET}

        was_error = False

        for f in self.fnames:
            try:
                xtable = in_dict[self.input_format](f)
                print('{}: Table succesfully read!'.format(f))
            except Exception as e:
                self.print_warning(e)

            # if no user-defined output dir use current
            if self.output_dir == '':
                self.output_dir = os.path.dirname(f)

            # set filename for output file
            fname = os.path.splitext(os.path.split(f)[1])[0] + '_' + self.input_format +\
                    '_to_' + self.output_format

            # generate output path w/o extension
            outpath = os.path.join(self.output_dir, fname)

            try:
                out_dict[self.output_format](xtable, outpath)
                print('{}: Table successfully written '.format(f) +
                      'to {}!'.format(outpath))
            except Exception as e:
                self.print_warning('Conversion of {} was '.format(f) +
                                   'not successfull:{}'.format(str(e)))
                was_error = True
                break

        if not was_error:
            QMessageBox.information(self,
                            'Success!',
                            'File(s) successfully written ' +
                            'to {}!'.format(outpath))
    ##################################
    # Dialogs
    ##################################

    @QtCore.pyqtSlot()
    def show_about(self):
        QMessageBox.about(self,
                          self.tr('About CroCo'),
                          'Version 0.1 (Oct 2017) <br><br>Written by <a href="mailto:jub@halomem.de">Julian Bender</a> at Martin Luther University Halle Wittenberg, Germany')

    def print_warning(self, error):
        QMessageBox.warning(self, "Error!",
                                str(error))

    ##################################
    # Other
    ##################################

    def center(self):
        """
        Center the window by moving it into the middle of the screen
        """

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    # redefinition of the internal closing event
    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = CroCo()
    sys.exit(app.exec_())