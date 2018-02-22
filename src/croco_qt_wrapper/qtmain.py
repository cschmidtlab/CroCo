#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Main Programme for for the CroCo cross-link converter GUI programme
"""

from PyQt5.QtWidgets import QWidget, QToolTip, QPushButton, QApplication,\
                            QMessageBox, QDesktopWidget, QGridLayout,\
                            QComboBox, QLabel, QFileDialog, QAction,\
                            QMainWindow, QHBoxLayout, QTabWidget, QLineEdit
import PyQt5.QtCore as QtCore

import os
import pandas as pd
import sys

croco_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

sys.path.append(os.path.abspath(os.path.join('..', croco_dir)))

from croco_qt_wrapper.ui.ui_mainwindow import Ui_MainWindow

# load croco
import croco

"""
Dict that contains all variables, options and settings that should be
passed within the programme
"""
Options = {'cv_input_format': 'pLink1',
           'cv_output_format': 'xTable',
           'cv_output_dir': '',
           'cv_fnames': []}

class CroCo_MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):
        # initialise the parent class
        super().__init__()
        # initiate itself
        self.setupUi(self)
        self.createConnects()
        # read Options
        global Options

        # set starting positions of dropdown menus
        self.cv_input_format = Options['cv_input_format']
        self.cv_output_format = Options['cv_output_format']
        self.cv_output_dir = Options['cv_output_dir']
        self.cv_fnames = Options['cv_fnames']

    #####################################
    # Connects
    #####################################

    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """

        #-----------------------------------------
        # Crosslink Converter
        #-----------------------------------------

        # Dropdown Menu for Input
        self.convert_input_dropdown.activated[str].\
            connect(self.set_conv_input_format)

        # Dropdown for output
        self.convert_output_dropdown.activated[str].\
            connect(self.set_conv_output_format)

        # Start button
        self.convert_start.clicked.connect(self.start_conversion)

        # call a close Event on clicking the button
        self.convert_quit.clicked.connect(self.close)

        # Open file dialog
        self.convert_load_btn.clicked.connect(self.set_conv_input_files)

        # Output files button
        self.convert_output_btn.clicked.connect(self.set_conv_outdir)

        #-------------------------------------------
        # Menu
        #-------------------------------------------

        # Quit Button
        self.actionQuit.triggered.connect(self.close)
        # About PopUp
        self.actionAbout.triggered.connect(self.show_about)

    #####################################
    # Definitions for Conversion
    #####################################

        # assign to variables
    def set_conv_input_format(self, text):
        print(text)
        self.cv_input_format = text

    def set_conv_output_format(self, text):
        self.cv_output_format = text

    # open folder dialog for pLink (only one folder at once) or file
    # dialog for others (multiple files at once)
    # decorator explicitely defines method as slot!
    @QtCore.pyqtSlot()
    def set_conv_input_files(self):
        if 'pLink' in self.cv_input_format:
            self.cv_fnames = [QFileDialog.getExistingDirectory(self,\
                                                'Select directory')]
        else:
            self.cv_fnames = QFileDialog.getOpenFileNames(self, 'Open file',\
                                                       '/home')[:-1][0]

        # avoid incorrect indexing if cv_fnames is empty (i.e. if user pressed cancel)
        if self.cv_fnames != '':
            # update input label
            self.convert_input_lbl.setText(os.path.basename(self.cv_fnames[0]))

    @QtCore.pyqtSlot()
    def set_conv_outdir(self):
        self.cv_output_dir = QFileDialog.getExistingDirectory(self,\
                                                'Select directory for output')
        # update output label
        self.convert_output_lbl.setText(os.path.basename(self.cv_output_dir))


    @QtCore.pyqtSlot()
    def start_conversion(self):
        """
        Collect all necessary information from self and start the
        conversion by calling the actual conversion script
        """

        print('Going to convert {} from {} '.format(', '.join(self.cv_fnames),
                                                   self.cv_input_format) +
               'format to {} format'.format(self.cv_output_format))

        in_dict = {'pLink1': croco.pLink1.Read,
                   'pLink2': croco.pLink2.Read,
                   'Kojak': croco.Kojak.Read,
                   'xQuest': croco.xQuest.Read,
                   'xTable': pd.read_csv}

        out_dict = {'xTable': croco.xTable.Write,
                    'xVis': croco.xVis.Write,
                    'xiNet': croco.xiNET.Write,
                    'DynamXL': croco.DynamXL.Write,
                    'xWalk': croco.xWalk.Write}

        was_error = False

        for f in self.cv_fnames:
            try:
                xtable = in_dict[self.cv_input_format](f)
                print('{}: Table succesfully read!'.format(f))
            except Exception as e:
                print_warning(self, e)

            # if no user-defined output dir use current
            if self.cv_output_dir == '':
                self.cv_output_dir = os.path.dirname(f)

            # set filename for output file
            fname = os.path.splitext(os.path.split(f)[1])[0] + '_' + self.cv_input_format +\
                    '_to_' + self.cv_output_format

            # generate output path w/o extension
            outpath = os.path.join(self.cv_output_dir, fname)

            try:
                out_dict[self.cv_output_format](xtable, outpath)
                print('{}: Table successfully written '.format(f) +
                      'to {}!'.format(outpath))
            except Exception as e:
                print_warning(self, 'Conversion of {} was '.format(f) +
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
                          'Version 0.1 (Dec 2017) <br><br>Written by <a href="mailto:jub@halomem.de">Julian Bender</a> at Martin Luther University Halle Wittenberg, Germany')

    # redefinition of the internal closing event
    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

def print_warning(self, error):
    QMessageBox.warning(self, "Error!",
                            str(error))
