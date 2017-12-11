# -*- coding:utf-8 -*-

import sys, os
from subprocess import call

if sys.platform == "win32":
    bindir = r"C:\ProgramData\Anaconda3\Library\bin"
else:
    bindir = "/usr/bin"

uic_path = os.path.join(bindir, 'pyuic5.bat')
rcc_path = os.path.join(bindir, 'pyrcc5.exe')

ui_path = 'ui'
rc_path = ""

out_path = 'src'

ui_files = { "CroCo_qt.ui": "ui_croco_qt.py" }
rc_files = {"croco.qrc": "croco_rc.py"}

for file in ui_files:
    file_path = os.path.join(ui_path, file)
    out_file_path = os.path.join(out_path, ui_files[file])
    call([uic_path, "-x", file_path, "-o", out_file_path])

for file in rc_files:
   file_path = os.path.join(rc_path, file)
   out_file_path = os.path.join(out_path, rc_files[file])
   call([rcc_path, file_path, '-o', out_file_path])