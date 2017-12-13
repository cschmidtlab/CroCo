# -*- coding:utf-8 -*-

import sys, os
from setuptools import setup, find_packages

# set proper paths to tcl and tk libraries
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

# see https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
setup(
    name = "croco_qt",
    version = "0.1",
    packages=find_packages('src'),  # include all packages under src
    package_dir={'':'src'},   # tell distutils packages are under src
    scripts = ['bin/croco_qt.py'],
    package_data = {'' : 'data'},
    #
    # can contain requirements such as 'docutils>=0.3'
    install_requires = [],
    #
    # metadata for upload to PyPI
    description = "Cross-Link files conversion software",
    author = "Julian Bender",
    author_email = "jub@halomem.de",
    url = "www.halomem.de",
    long_description = """ The CroCo cross-link conversion engine aims to 
    simplify the integration of XL-MS data from different cross-link annotation
    programmes."""
    )  