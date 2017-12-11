# -*- coding:utf-8 -*-

from distutils.core import setup

setup(
    name = "croco",
    version = "0.1",
    description = "Cross-Link files conversion software",
    author = "Julian Bender",
    author_email = "jub@halomem.de",
    url = "www.halomem.de",
    package_dir = {"": "src"},
    packages = ["croco"],
    scripts = ["bin/croco"],
    long_description = """ The CroCo cross-link conversion engine aims to 
    simplify the integration of XL-MS data from different cross-link annotation
    programmes."""
    )