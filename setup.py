# -*- coding:utf-8 -*-

import sys, os
from setuptools import setup, find_packages

# set proper paths to tcl and tk libraries
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# see https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
setup(
    name = "croco_qt",
    version = "0.1",
    packages=find_packages('croco'),  # include all packages under src
    #package_dir={'':'src'},   # tell distutils packages are under src
    scripts = ['bin/croco_qt.py',
               'bin/croco_cl.py'],
    data_files=[('data/images', ['data/images/logo.png']),
                  ('data/config', ['data/config/modification.ini'])],
    #
    # can contain requirements such as 'docutils>=0.3'
    install_requires = [],
    #
    # metadata for upload to PyPI
    description = "Cross-Link files conversion software",
    author = "Julian Bender",
    author_email = "jub@halomem.de",
    url = "www.halomem.de",
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ]
    )  
