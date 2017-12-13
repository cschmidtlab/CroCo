# -*- coding:utf-8 -*-

from distutils import setup


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
    package_dir = { "" : 'src'},
    packages = ['croco_qt']
    scripts = ['bin/croco']
    )