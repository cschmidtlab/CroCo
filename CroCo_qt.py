# -*- coding: utf-8 -*-
"""
Introduction to PyQt5

from http://zetcode.com/gui/pyqt5/firstprograms/
"""

from PyQt5.QtWidgets import QWidget, QToolTip, QPushButton, QApplication,\
                            QMessageBox, QDesktopWidget, QGridLayout,\
                            QComboBox, QLabel, QFileDialog
from PyQt5.QtGui import QFont, QIcon

import os
os.getcwd()

import sys
sys.path.append('./lib')
import xlink_parser as xp

import pandas as pd

class CroCo(QWidget):
    
    # initialise properties
    input_format = 'pLink'
    output_format = 'xTable'
    fnames = []
    
    def __init__(self):
        super().__init__()
        
        self.initUI()
        
    def initUI(self):

        self.setWindowTitle('Crosslink Converter')
        self.setWindowIcon(QIcon('./python-icon.png')) # icon in ubuntu is
                                                       # displayed on sidepanel        
        QToolTip.setFont(QFont('SansSerif', 10))        
        
        # Define Dropdown Menu for Input
        input_lbl = QLabel('Input file format', self)
        input_format = QComboBox(self)
        input_format.addItem("pLink")
        input_format.addItem("Kojak")
        input_format.addItem("xTable")   
        input_format.activated[str].connect(self.set_input_format)

        # Dropdown for output
        output_lbl = QLabel('Output file format', self)
        output_format = QComboBox(self)
        output_format.addItem("xTable")

        output_format.activated[str].connect(self.set_output_format)
        
        # Start button
        sbtn = QPushButton('Start', self)
        sbtn.setToolTip('Start the conversion')
        sbtn.clicked.connect(self.start_conversion)
        sbtn.resize(sbtn.sizeHint())

        # Close button
        qbtn = QPushButton('Quit', self)
        # call a close Event on clicking the button
        qbtn.clicked.connect(self.close)
        qbtn.resize(qbtn.sizeHint())

        # Open file dialog
        fbtn = QPushButton('Load file(s)', self)
        fbtn.clicked.connect(self.set_input_files)

        # Define Layout
        grid = QGridLayout()
        grid.setSpacing(10)
        
        grid.addWidget(input_lbl, 0, 0)
        grid.addWidget(input_format, 0, 1)
        grid.addWidget(output_lbl, 0, 2)
        grid.addWidget(output_format, 0, 3)
        
        grid.addWidget(fbtn, 1, 0)
        
        grid.addWidget(sbtn, 2, 0)
        grid.addWidget(qbtn, 2, 3)

        self.resize(250, 150)
        self.center()
        self.setLayout(grid) 
        self.show()

    def set_input_format(self, text):
        self.input_format = text

    def set_output_format(self, text):
        self.output_format = text

    def set_input_files(self):
        if self.input_format == 'pLink':
            self.fnames = [QFileDialog.getExistingDirectory(self, 'Select directory')[:-1]]
        else:
            self.fnames = QFileDialog.getOpenFileNames(self, 'Open file', '/home')[:-1][0]

    def center(self):
                
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

    def print_warning(self, error):
        QMessageBox.warning(self, "Error!",
                                error)
    def start_conversion(self):
        QMessageBox.information(self, "Information on variables",
                                "Going to convert {} from {} format to {} format".format(', '.join(self.fnames),
                                                                                         self.input_format,
                                                                                         self.output_format))
        
        in_dict = {'pLink': xp.ReadpLink,
                   'Kojak': xp.ReadKojak,
                   'xTable': pd.read_csv}
        
        out_dict = {'xTable': xp.WriteXtable}
        
        for f in self.fnames:
            try:
                xtable = in_dict[self.input_format](f)
            except Exception as e:
                self.print_warning(e)
            out_dict[self.output_format](xtable)

if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = CroCo()
    sys.exit(app.exec_())