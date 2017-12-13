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
    croco_dir = os.path.join(os.path.dirname(croco_script), '..')
    # resolve relative path
    prefix = os.path.abspath(os.path.normpath(croco_dir))
    
    src_dir = os.path.join(prefix, 'src')
    data_dir = os.path.join(prefix, 'data')
    share_dir = os.path.join(prefix, 'share', 'croco')
    
    if os.path.exists(src_dir):
        # started from local directory, not installed
        sys.stderr.write('Using modules from ' + src_dir + '\n')
        sys.path.append(src_dir)
    else:
        if os.name == 'nt':
            # installed on windows
            sys.stderr.write('Installed on Win \n')
        else:
            # started under *nix or Mac
            sys.stderr.write('Installed on *nix or Mac \n')

    # use relative paths from within data to simplify programme
    if os.path.exists(data_dir):
        os.chdir(data_dir)
        
    app = QApplication(argv)
    
    # Translation belongs here
    
    import croco.ui.mainwindow as mainwin
    
    mainwindow = mainwin.CroCo_MainWindow()
    mainwindow.show()
    
    sys.exit(app.exec_())

print(__name__)

if __name__.endswith('__main__'):
    main(sys.argv)