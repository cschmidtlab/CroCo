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

# pyQt for matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import sys, os
import pandas as pd
import numpy as np
import itertools as it

from croco.ui.ui_mainwindow import Ui_MainWindow
from croco.ui.ui_spectrum import Ui_SpectrumAssignment
from croco.ui.ui_spectrumoptions import Ui_SpectrumAssignmentOptions

import croco.reader as xr
import croco.writer as xw

import croco.SpectrumAssignment as assign
import croco.SpectrumReader as assignr


"""
Dict that contains all variables, options and settings that should be
passed within the programme
"""
Options = {'cv_input_format': 'pLink',
           'cv_output_format': 'xTable',
           'cv_output_dir': '',
           'cv_fnames': [],
           
           'as_input_format': 'pLink',
           'as_output_format': 'xTable',
           'as_output_dir': '',
           'as_fnames': [],
           
           'as_xtable_path': None,
           'as_mgf_path': None,
           'as_outdir': None,
           
           'as_filehandle': None,
           'as_xtable': None,
       
           'xlink_mass': 138.068,
       
           'min_charge': 1,
           'max_charge': 2,
           'ionmass': 'monoisotopic',
           'iontypes': ['a', 'b', 'y'],
           'max_mz': 2000,
           'min_ppm': -50,
           'max_ppm': 50,
           
           'score_threshold': None,
           'score_threshold_direction': 'Less than',
           
           'modifications': None
            }

class CroCo_MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):
        # initialise the parent classe
        super().__init__()
        # initial itself
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
        self.assign_load_xlink_btn.clicked.connect(self.set_as_xtable_path)

        # Open mgf file dialog
        self.assign_load_mgf_btn.clicked.connect(self.set_as_mgf_path)

        # Output files button
        self.assign_output_btn.clicked.connect(self.set_as_outdir)
        
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
        print(text)
        self.cv_input_format = text

    def set_conv_output_format(self, text):
        self.cv_output_format = text

    # open folder dialog for pLink (only one folder at once) or file
    # dialog for others (multiple files at once)
    # decorator explicitely defines method as slot!
    @QtCore.pyqtSlot()
    def set_conv_input_files(self):
        if self.cv_input_format == 'pLink':
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

        in_dict = {'pLink': xr.ReadpLink,
                   'Kojak': xr.ReadKojak,
                   'xTable': pd.read_csv}

        out_dict = {'xTable': xw.WriteXtable,
                    'xVis': xw.WritexVis,
                    'xiNet': xw.WritexiNET,
                    'DynamXL': xw.WriteDynamXL,
                    'xWalk': xw.WritexWalk}

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

    #####################################
    # Definitions for Assignment
    #####################################

    @QtCore.pyqtSlot()
    def set_as_xtable_path(self):
        self.as_xtable_path = QFileDialog.getOpenFileName(self,
                                                            'Open xTable file')[0]
                                                            
        self.assign_xlink_lbl.setText(os.path.basename(self.as_xtable_path))
        
        Options['as_xtable_path'] = self.as_xtable_path
        
    @QtCore.pyqtSlot()        
    def set_as_mgf_path(self):
        self.as_mgf_path = QFileDialog.getOpenFileName(self,
                                                            'Open MGF file')[0]
        self.assign_mgf_lbl.setText(os.path.basename(self.as_mgf_path))
        
        Options['as_mgf_path'] = self.as_mgf_path
                                
    @QtCore.pyqtSlot()
    def set_as_outdir(self):
        self.as_outdir = QFileDialog.getExistingDirectory(self,\
                                                'Select directory for output')
        self.assign_outdir_lbl.setText(os.path.basename(self.as_outdir))

        Options['as_outdir'] = self.as_outdir

    @QtCore.pyqtSlot()
    def start_assignment(self):
        
        ready = True
        warnings = []
        
        # check if all necessary files are present
        try:
            Options['as_filehandle'] =  open(self.as_mgf_path, 'r')
        except Exception as e:
            warnings.append('Please provide a valid MGF input file!')
            ready = False
        
        try:
            Options['as_xtable'] =  pd.read_excel(self.as_xtable_path)

        except Exception as e:
            warnings.append('Please provide a valid xTable input file!')
            ready = False
        
        if ready:
            self.assignwin = AssignmentWindow() # calls init
            self.assignwin.show() # shows windows
        else:
            print_warning(self, '\n\n'.join(warnings))

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

class AssignmentWindow(QMainWindow, Ui_SpectrumAssignment):
    
    def __init__(self):
        # initialise the parent classe
        super().__init__()
        # create Widgets that are then used in setupUI
        self.createFigure()
        # initial itself
        self.setupUi(self)
        # read options
        global Options
        
        self.createConnects()
        self.openMGF()
        self.iterFromxTable()
        self.loadSpectrum()

    #####################################
    # Connects
    #####################################

    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """
        self.assign_options.triggered.connect(self.showAssignmentOptions)
        
        self.assign_next_btn.clicked.connect(lambda: self.loadSpectrum())
        self.assign_previous_btn.clicked.connect(lambda: self.loadSpectrum(direction='prev'))
        
        self.assign_score_decline_btn.clicked.connect(lambda: self.rateSpectrum(score=0))
        self.assign_score_unsure_btn.clicked.connect(lambda: self.rateSpectrum(score=1))
        self.assign_score_accept_btn.clicked.connect(lambda: self.rateSpectrum(score=2))
        self.assign_score_good_btn.clicked.connect(lambda: self.rateSpectrum(score=3))

        self.assign_save_btn.clicked.connect(lambda: self.saveTable())

        self.assign_filter_btn.clicked.connect(self.filterByScore)

    def createFigure(self):
        """
        Define all widgets necessary to plot the spectrum
        """
        # a figure instance to plot on
        self.figure = Figure()
        
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

    #####################################
    # Definitions for Assignment
    #####################################
    @QtCore.pyqtSlot()
    def openMGF(self):
        """
        Opens MGF and calls the MGF-Indexer with an open as_filehandle
        """

        # call indexer
        self.spectrum2offset = assignr.IndexMGF(Options['as_filehandle'])

    def iterFromxTable(self):
        """
        Reads xTable input file and converts data to an itertools cycle
        with all relevant information for spectrum assignment
        """
        
        class xiter(list):
            """
            Iterator class to cycle over all xlinks read from an xTable
            document. Allows reverse iteration for going back to the
            previous spectrum
            """
            def __init__(self, list):
                self.list = list
                self.end = len(list) - 1
                self.index = -1

                global Options
                self.toShow()
                
            def next(self):
                """
                returns the next element of a list and the first if the last
                was already reached
                """
                if self.index == self.end:
                    # reset index in case the full list has passes
                    self.index = 0
                else:
                    self.index += 1
                
                # only return the next spectrum if it is not filtered out
                if self.toShow():
                    return self.list[self.index]
                else:
                    print('Skipping Spectrum:' + str(self.list[self.index][0]))
                    return self.next()
                    
                
            def prev(self):
                """
                inverse of next
                """
                if self.index == 0:
                    self.index = self.end
                else:
                    self.index -= 1

                if self.toShow():
                    return self.list[self.index]
                else:
                    print('Skipping Spectrum:' + str(self.list[self.index][0]))
                    return self.prev()

            def toShow(self):
                """
                Checks if a potential next spectrum fulufills all requirements
                i.e. score and not yet annotated
                """
                
                # only check for scores if max score is set
                if Options['score_threshold']:
                    if Options['score_threshold_direction'] == 'Less than':
                        if self.list[self.index][-2] < Options['score_threshold']:
                            return True
                        else:
                            return False
                    elif Options['score_threshold_direction'] == 'Greater than':
                        if self.list[self.index][-2] > Options['score_threshold']:
                            return True
                        else:
                            return False
                else:
                    return True

            def setCurrentScore(self, score):
                """
                Method to access the last element of the
                list of lists i.e. the manual score from outside
                the class
                """
                self.list[self.index][-1] = score
        
        self.xtable = Options['as_xtable']
        
        # create a column "Manual score" if not already present
        if not 'Manual score' in self.xtable.columns:
            self.xtable['Manual score'] = np.nan
         
        xlinks = self.xtable[['pepseq1',
                              'xlink1',
                              'pepseq2',
                              'xlink2',
                              'xtype',
                              'scanno',
                              'prec_ch',
                              'score',
                              'Manual score']].values

        indices = self.xtable.index.values
        
        # concatenate 1d array indices with 2d array xlinks
        xlinks = np.hstack((indices[:, None], xlinks))

        self.xlink_iter = xiter(xlinks)
        
    @QtCore.pyqtSlot()
    def loadSpectrum(self, direction='next'):
        """
        Generates a matplotlib figure object from a mgf file and a cross-link
        sequence
        """
        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        ax.clear() 
      
        if direction == 'next':
            [self.idx,
             pepseq1,
             xlink1,
             pepseq2,
             xlink2,
             xtype,
             scanno,
             prec_ch,
             score,
             man_score] = self.xlink_iter.next()
        else:
            [self.idx,
             pepseq1,
             xlink1,
             pepseq2,
             xlink2,
             xtype,
             scanno,
             prec_ch,
             score,
             man_score] = self.xlink_iter.prev()  

        mz2intens = assignr.ReadSpectrum(scanno, # spectrum
                                         Options['as_filehandle'], # file handle
                                         self.spectrum2offset) # spectrum dict

        ions2desc =  assign.IonsFromXlinkSequence(pepseq1,
                                                  xlink1,
                                                  pepseq2,
                                                  xlink2,
                                                  Options['xlink_mass'],
                                                  Options['iontypes'],
                                                  [Options['min_charge'],
                                                   Options['max_charge']],
                                                  Options['ionmass'],
                                                  Options['modifications'],
                                                  Options['max_mz'])
                
        assignment_error = assign.AssignAndPlotPSM(mz2intens,
                                                   ions2desc,
                                                   [Options['min_ppm'],
                                                    Options['max_ppm']],
                                                   ax=ax)

        print('Current score: {}'.format(score))

        # refresh canvas
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def filterByScore(self):
        """
        Read min or max score from GUI and store the respective string in
        the Options dict. Filtering is then performed by the iterator class
        """
        
        try:
            score = self.assign_filter_edit.text()
            Options['score_threshold'] = float(score)
            Options['score_threshold_direction'] = str(self.assign_filter_dropdown.currentText())
        except:
            print_warning('Please provde a proper number for score: {}'.format(score))

    @QtCore.pyqtSlot()
    def rateSpectrum(self, score=1):
        """
        Take the user input, assign the respective score value to a
        pd-DataFrame, and call the next spectrum
        """
        self.xtable.set_value(self.idx, 'Manual score', score)
        self.xlink_iter.setCurrentScore(score)
        
        print('You rated this spectrum with score {}'.format(score))
        
        self.loadSpectrum()

    @QtCore.pyqtSlot()
    def saveTable(self):
        """
        Saves the xtable on self with its input filename
        """
        try:
            self.xtable.to_excel(Options['as_xtable_path'])
            print('Table saved')
        except Exception as e:
            print_warning(self, 'There was an error saving your xTable file: {}'.\
                format(e))

    @QtCore.pyqtSlot()
    def showAssignmentOptions(self):
        """
        Opens the window with the options for spectrum assignment
        """
        self.assign_optionswin = AssignmentOptionsWindow() # calls init
        self.assign_optionswin.show() # shows window

    # redefinition of the internal closing event
    def closeEvent(self, event):
        try:
            # close the file when closing the window
            Options['as_filehandle'].close()
            event.accept()
        except:
            event.ignore()

            
class AssignmentOptionsWindow(QMainWindow, Ui_SpectrumAssignmentOptions):
    def __init__(self):
        # initialise the parent class
        super().__init__()
        # read Options
        global Options
        
        # initiate itself
        self.setupUi(self)
        self.createConnects()


    def createConnects(self):
        """
        Generates the logical wiring of the programme.
        Every input gets defined and assigned to an output
        """
        
        # Set the properties in the AssignmentWindow class explicitely
        self.ionmass = self.options_ionmass.checkedButton().text()

        self.iontypes = Options['iontypes']

        def set_options_ions(iontype):
            """
            Modifies the iontype list containing all ions that should be searched
            """
            if iontype not in self.iontypes:
                self.iontypes.append(iontype)
            else:
                self.iontypes.remove(iontype)

        # set the connections for all ion-types
        # use the lambda operator as connects expects a callable funtion
        # i.e. no nested input
        self.options_aion.stateChanged.connect(lambda: set_options_ions('a'))
        self.options_bion.stateChanged.connect(lambda: set_options_ions('b'))
        self.options_cion.stateChanged.connect(lambda: set_options_ions('c'))
        self.options_xion.stateChanged.connect(lambda: set_options_ions('x'))
        self.options_yion.stateChanged.connect(lambda: set_options_ions('y'))
        self.options_zion.stateChanged.connect(lambda: set_options_ions('z'))

    def saveOptions(self):
        pass
        
    def cancelOptions(self):
        pass


def print_warning(self, error):
    QMessageBox.warning(self, "Error!",
                            str(error))