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

from ui_croco_qt import Ui_MainWindow

import os
os.getcwd()

import sys
import CroCo_reader as xr
import CroCo_writer as xw

import pandas as pd

class CroCo_MainWindow(QMainWindow, Ui_MainWindow):

    # initialise properties
    input_format = 'pLink'
    output_format = 'xTable'
    output_dir = ''
    fnames = []

    def __init__(self):
        # initialise the parent classe
        super().__init__()
        # initial itself
        self.setupUi(self)
        self.createConnects()

    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """

        #############################################
        # Crosslink Converter Part
        #############################################

        # Define Dropdown Menu for Input
        self.convert_input_dropdown.activated[str].connect(self.set_input_format)

        # Dropdown for output
        self.convert_output_dropdown.activated[str].connect(self.set_output_format)

        # Start button
        self.convert_start.clicked.connect(self.start_conversion)

        # call a close Event on clicking the button
        self.convert_quit.clicked.connect(self.close)

        # Open file dialog
        self.convert_load_btn.clicked.connect(self.set_input_files)

        # Output files button
        self.convert_output_btn.clicked.connect(self.set_output_dir)

        #############################################
        # Spectrum annotation
        #############################################

        # Open xlink file dialog

        # Open mgf file dialog

        # Output files button
        
        # Start button

        # PSM annotation

        # Navigation

        #############################################
        # Menu
        #############################################
        
        # Quit Button
        self.actionQuit.triggered.connect(self.close)
        # About PopUp
        self.actionAbout.triggered.connect(self.show_about)
        
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
        self.convert_input_lbl.setText(os.path.basename(self.fnames[0]))
        
    @QtCore.pyqtSlot()
    def set_output_dir(self):
        self.output_dir = QFileDialog.getExistingDirectory(self,\
                                                'Select directory for output')
        # update output label
        self.convert_output_lbl.setText(os.path.basename(self.output_dir))

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

    # redefinition of the internal closing event
    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()