# -*- coding: utf-8 -*-
"""

"""

from PyQt5.QtWidgets import QWidget, QToolTip, QPushButton, QApplication,\
                            QMessageBox, QDesktopWidget, QGridLayout,\
                            QComboBox, QLabel, QFileDialog, QAction,\
                            QMainWindow
from PyQt5.QtGui import QFont, QIcon

import PyQt5.QtCore as QtCore

import os
os.getcwd()

import sys
sys.path.append('./lib')
import reader as xr
import writer as xw

import pandas as pd

class CroCo(QMainWindow):

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
        self.setWindowIcon(QIcon('./python-icon.png')) # icon in ubuntu is
                                                       # displayed on sidepanel
        QToolTip.setFont(QFont('SansSerif', 10))

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

    @QtCore.pyqtSlot()
    def show_about(self):
        QMessageBox.about(self,
                          self.tr('About CroCo'),
                          'Version 0.1 (Oct 2017) <br><br>Written by <a href="mailto:jub@halomem.de">Julian Bender</a> at Martin Luther University Halle Wittenberg, Germany')

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

        # Define Layout as grid
        grid = QGridLayout()
        grid.setSpacing(10)
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1,1)
        grid.setRowStretch(2,1)
        grid.setRowStretch(3,2)
        grid.setRowStretch(4,1)


        # add the 1st line of elements to the grid
        grid.addWidget(self.input_lbl, 0, 0)
        grid.addWidget(self.input_dropdown, 0, 1)
        grid.addWidget(self.output_lbl, 0, 2)
        grid.addWidget(self.output_dropdown, 0, 3)

        # 2nd line
        grid.addWidget(self.fbtn, 1, 1)
        grid.addWidget(self.obtn, 1, 3)

        #3rd and 4th line
        grid.addWidget(self.inpath_title, 2, 0)
        grid.addWidget(self.inpath_lbl, 3, 0)
        grid.addWidget(self.outpath_title, 2, 2)
        grid.addWidget(self.outpath_lbl, 3, 2)

        # 5th line
        grid.addWidget(self.sbtn, 4, 0)
        grid.addWidget(self.qbtn, 4, 3)

        # generate widget for main window
        mainWidget = QWidget()
        # apply layout
        mainWidget.setLayout(grid)
        # assign widget to main window (i.e. self)
        self.setCentralWidget(mainWidget)
        
        # set window size
        self.resize(250, 200)
        # center
        self.center()
        # show the window
        self.show()


    def center(self):
        """
        Center the window by moving it into the middle of the screen
        """

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

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

        out_dict = {'xTable': xw.WriteXtable}

        was_error = False

        for f in self.fnames:
            try:
                xtable = in_dict[self.input_format](f)
                print('{}: Table succesfully read!'.format(f))
            except Exception as e:
                self.print_warning(e)

            # if no user-defined output dir use current
            if self.output_dir == '':
                self.output_dir = os.path.basename(f)

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
                                   'not successfull:'.format(str(e)))
                was_error = True
                break

        if not was_error:
            QMessageBox.information(self,
                            'Success!',
                            'File(s) successfully written ' +
                            'to {}!'.format(outpath))

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

if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = CroCo()
    sys.exit(app.exec_())