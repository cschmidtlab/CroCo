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

import sys, os
import pandas as pd
import itertools as it

from croco.ui.ui_mainwindow import Ui_MainWindow
from croco.ui.ui_spectrum import Ui_SpectrumAssignment
from croco.ui.ui_spectrumoptions import Ui_SpectrumAssignmentOptions
import croco.reader as xr
import croco.writer as x

class CroCo_MainWindow(QMainWindow, Ui_MainWindow):

    #####################################
    # Options preset
    #####################################
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

    #####################################
    # Connects
    #####################################

    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """

        #-----------------------------------------
        # Crosslink Converter Part
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

        #-----------------------------------------
        # Spectrum annotation
        #-----------------------------------------

        # Open xlink file dialog
        self.assign_load_xlink_btn.clicked.connect(self.set_assign_xtable)

        # Open mgf file dialog
        self.assign_load_mgf_btn.clicked.connect(self.set_assign_mgf)

        # Output files button
        self.assign_output_btn.clicked.connect(self.set_assign_outdir)
        
        # Start button
        self.assign_start_btn.clicked.connect(self.start_assignment)

        # PSM annotation

        # Navigation

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
        self.input_format = text

    def set_conv_output_format(self, text):
        self.output_format = text

    # open folder dialog for pLink (only one folder at once) or file
    # dialog for others (multiple files at once)
    # decorator explicitely defines method as slot!
    @QtCore.pyqtSlot()
    def set_conv_input_files(self):
        if self.input_format == 'pLink':
            self.fnames = [QFileDialog.getExistingDirectory(self,\
                                                'Select directory')]
        else:
            self.fnames = QFileDialog.getOpenFileNames(self, 'Open file',\
                                                       '/home')[:-1][0]

        # avoid incorrect indexing if fnames is empty (i.e. if user pressed cancel)
        if self.fnames != '':
            # update input label
            self.convert_input_lbl.setText(os.path.basename(self.fnames[0]))
        
    @QtCore.pyqtSlot()
    def set_conv_outdir(self):
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

    #####################################
    # Definitions for Assignment
    #####################################

    @QtCore.pyqtSlot()
    def set_assign_xtable(self):
        self.assign_xtable = QFileDialog.getOpenFileName(self,
                                                            'Open xTable file')[0]
                                                            
        self.assign_xlink_lbl.setText(os.path.basename(self.assign_xtable))
        
    @QtCore.pyqtSlot()        
    def set_assign_mgf(self):
        self.assign_mgf = QFileDialog.getOpenFileName(self,
                                                            'Open MGF file')[0]
        self.assign_mgf_lbl.setText(os.path.basename(self.assign_mgf))
                                
    @QtCore.pyqtSlot()
    def set_assign_outdir(self):
        self.assign_outdir = QFileDialog.getExistingDirectory(self,\
                                                'Select directory for output')
        self.assign_outdir_lbl.setText(os.path.basename(self.assign_outdir))

    @QtCore.pyqtSlot()
    def start_assignment(self):
        self.assignwin = AssignmentWindow() # calls init
        self.assignwin.show() # shows windows
        # TODO Load first spectrum


    ##################################
    # Dialogs
    ##################################

    @QtCore.pyqtSlot()
    def show_about(self):
        QMessageBox.about(self,
                          self.tr('About CroCo'),
                          'Version 0.1 (Dec 2017) <br><br>Written by <a href="mailto:jub@halomem.de">Julian Bender</a> at Martin Luther University Halle Wittenberg, Germany')

class AssignmentWindow(QMainWindow, Ui_SpectrumAssignment):
    
    ##################################
    # Options presets
    ################################## 
    ionmass = 'monoisotopic'
    iontypes = ['a', 'b', 'y']
    maxmass = 2000
    
    def __init__(self):
        # initialise the parent classe
        super().__init__()
        # initial itself
        self.setupUi(self)
        self.createConnects()

    #####################################
    # Connects
    #####################################

    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """
        self.assign_options.triggered.connect(self.showAssignmentOptions)

    #####################################
    # Definitions for Assignment
    #####################################
    
    def iterFromxTable(self):
        """
        Reads xTable input file and converts data to an itertools cycle
        with all relevant information for spectrum assignment
        """
        
        try:
            self.xtable = pd.read_excel(self.assign_xtable)
        except:
            print_warning('Please provide a valid xTable input file')
        
        xlinks = data['pepseq1',
                      'xlink1',
                      'pepseq2',
                      'xlink2',
                      'xtype',
                      'scanno',
                      'prec_ch'].tolist()

        self.xlink_iter = it.cycle(xlinks)
    
    @QtCore.pyqtSlot()
    def showAssignmentOptions(self):
        """
        Opens the window with the options for spectrum assignment
        """
        self.assignoptions = AssignmentOptionsWindow() # calls init
        self.assignoptions.show() # shows window

class AssignmentOptionsWindow(QMainWindow, Ui_SpectrumAssignmentOptions):
    def __init__(self):
        # initialise the parent class
        super().__init__()
        # initial itself
        self.setupUi(self)
        self.createConnects()

    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """
        
        # Set the properties in the AssignmentWindow class explicitely
        AssignmentWindow.ionmass = self.options_ionmass.checkedButton().text()

        def set_options_ions(iontype):
            """
            Modifies the iontype list containing all ions that should be searched
            """
            if iontype not in AssignmentWindow.iontypes:
                AssignmentWindow.iontypes.append(iontype)
            else:
                AssignmentWindow.iontypes.remove(iontype)

        # set the connections for all ion-types
        # use the lambda operator as connects expects a callable funtion
        # i.e. no nested input
        self.options_aion.stateChanged.connect(lambda: set_options_ions('a'))
        self.options_bion.stateChanged.connect(lambda: set_options_ions('b'))
        self.options_cion.stateChanged.connect(lambda: set_options_ions('c'))
        self.options_xion.stateChanged.connect(lambda: set_options_ions('x'))
        self.options_yion.stateChanged.connect(lambda: set_options_ions('y'))
        self.options_zion.stateChanged.connect(lambda: set_options_ions('z'))


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