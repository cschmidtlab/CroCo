# -*- coding:utf-8 -*-

import sys, os
from setuptools import setup

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
    version = "0.3",
    packages=['croco', 'croco_wx_wrapper'],  # find all packages under src
    package_dir={'croco': 'src/croco',
                 'croco_wx_wrapper': 'src/croco_wx_wrapper'},
    package_data={'croco': ['data/*', 'lib/*']},
    entry_points={
        'gui_scripts': [
            'croco_wx = croco_wx_wrapper.__main__:main',
        ]
    },
 
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
