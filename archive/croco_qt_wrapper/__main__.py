#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
The CroCo Cross-Link Converter GUI

Graphical interface to convert results from data analysis of chemical cross-linking /
mass-spectrometry experiments.

This script checks the filesystem and calls the main programme
"""

import sys, os
from PyQt5.QtWidgets import QApplication
import PyQt5.QtCore as QtCore

def main(argv):
    
    croco_script = sys.argv[0]
    # check if programme was called via symlink
    if os.path.islink(croco_script):
        croco_script = os.readlink(croco_script)
    # dir is the directory above the bin-dir
    croco_dir = os.path.abspath(os.path.join(os.path.dirname(croco_script), '..'))
    
    sys.path.append(os.path.abspath(os.path.join('..', croco_dir)))
        
    app = QApplication(argv)
    
    # Translation belongs here
    
    print('Running at {}'.format(os.getcwd()))
    from croco_qt_wrapper import qtmain
    
    mainwindow = qtmain.CroCo_MainWindow()
    mainwindow.show()
    
    sys.exit(app.exec_())

# print(__name__)

if __name__.endswith('__main__'):
    main(sys.argv)