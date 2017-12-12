# -*- coding:utf-8 -*-

import sys, os
from cx_Freeze import setup, Executable

paths = sys.path
paths.append('src')  # path to look for modules (plus sys.path)

buildOptions = {'build_exe':
                    {'includes' : [],
                     'path' : paths,
                     'includes' : ['croco_qt', 'codecs'],
                     'excludes': [],
                     'include_files': []
                    }
                }

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
#if sys.platform == "win32":
#    base = "Win32GUI"

# set proper paths to tcl and tk libraries
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

setup(
    name = "croco",
    version = "0.1",
    description = "Cross-Link files conversion software",
    author = "Julian Bender",
    author_email = "jub@halomem.de",
    url = "www.halomem.de",
    long_description = """ The CroCo cross-link conversion engine aims to 
    simplify the integration of XL-MS data from different cross-link annotation
    programmes.""",
    options = buildOptions,
    executables = [Executable("bin/croco.py", base=base)]
    )  