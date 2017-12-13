# -*- coding:utf-8 -*-
"""
Generate UI script for teh CroCo cross-link converter

This script bundles call to the PyQt5 compilers pyuic and pyrcc and
automatically generates usable Python files from QDesigner .ui and .qrc
documents.

It additionally corrects for an import path error generated by QDesigner
"""

import sys, os
import re
from subprocess import call

# needed to adjustt the import path in the ui-derived python-file
package_name = 'croco.ui'

# paths to bins
if sys.platform == "win32":
    bindir = r"C:\ProgramData\Anaconda3\Library\bin"
else:
    bindir = "/usr/bin"
    
uic_path = os.path.join(bindir, 'pyuic5.bat')
rcc_path = os.path.join(bindir, 'pyrcc5.exe')

# path of source files within the project
ui_path = 'qt'
rc_path = ""
# paths to write the compiled files to
ui_out_path = os.path.join('src', 'croco', 'ui')
# dict to set names for the conversion
ui_files = { "CroCo_qt.ui": "ui_mainwindow.py",
             "spectrum.ui": "ui_spectrum.py",
             "spectrumOptions.ui": "ui_spectrumoptions.py"}
rc_files = {"croco.qrc": "croco_rc.py"}

for file in ui_files:
    # generate file paths for subprocess call
    file_path = os.path.join(ui_path, file)
    out_file_path = os.path.join(ui_out_path, ui_files[file])
    # call the puic compiler
    call([uic_path, "-x", file_path, "-o", out_file_path])
    
    # Correcting path from QDesigner
    buffer = []
    with open(out_file_path, 'r') as i:
        # QDesigner generates an import command FILE_rc that has to be
        # relative to the package i.e. croco.FILE_rc
        pattern = re.compile(r'import (\w+)_rc$')
        for line in i:
            if pattern.match(line):
                group = pattern.match(line)
                line = 'import {}.{}_rc'.format(package_name, group[1])
            buffer.append(line)
    with open (out_file_path, 'w') as o:
        o.writelines(buffer)

# do the same for the resources compiler
for file in rc_files:
   file_path = os.path.join(rc_path, file)
   out_file_path = os.path.join(ui_out_path, rc_files[file])
   call([rcc_path, file_path, '-o', out_file_path])