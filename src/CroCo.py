#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
The CroCo Cross-Link Converter console

Convert results from data analysis of chemical cross-linking /
mass-spectrometry experiments.
"""

import argparse

import sys
sys.path.append('./lib')
import reader as xr
import writer as xw

import pandas as pd
import os

description = """The CroCo cross-link converter:
-------------------------------
Convert results from data analysis of chemical cross-linking mass-spectrometry experiments."""

epilog = """EXAMPLES:
python CroCo.py pLink testdata\plink\2.report\sample1 xTable:
Convert a pLink results file (as a folder) into an xTable file. The output will
be written to the folder containing the input file.

python CroCo.py -o testdata pLink testdata\plink\2.report\sample* xTable:
Convert all samples residing in the report folder to xTable and write
the results directly into testdata"""

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=description,
                                 epilog=epilog)

parser.add_argument("in_format", help="Format of the input (pLink, kojak, xTable)")

# TODO improve to directly read list
parser.add_argument("infiles",
                    help="Comma separated list of input file(s) to convert",
                    type=str)

parser.add_argument("out_format", help="Format to convert the file to (xTable, xVis, xiNet)")

parser.add_argument("-o",
                    "--outdir",
                    help="(Optional) Directory for the output file",
                    type=str)

args = parser.parse_args()

def print_warning(error):
    print("An error occurred: %s"%error)

in_dict = {'pLink': xr.ReadpLink,
           'Kojak': xr.ReadKojak,
           'xTable': pd.read_csv,
           'xQuest': xr.ReadxQuest}

out_dict = {'xTable': xw.WriteXtable,
            'xVis': xw.WritexVis,
            'xiNet': xw.WritexiNET,
            'dynamXL': xw.WriteDynamXL}
                        
infiles = list(args.infiles.split(','))
            
for f in infiles:
    try:
        xtable = in_dict[args.in_format](f)
        print('{}: Table succesfully read!'.format(f))
    except Exception as e:
        print_warning(e)
        break

    # if no user-defined output dir use current
    if not args.outdir:
        outdir = os.path.dirname(f)
    else:
        outdir = args.outdir

    # set filename for output file
    fname = os.path.splitext(os.path.split(f)[1])[0] + '_' + args.in_format +\
            '_to_' + args.out_format

    # generate output path w/o extension
    outpath = os.path.join(outdir, fname)

    try:
        out_dict[args.out_format](xtable, outpath)
        print('{}: Table successfully written '.format(f) +
                'to {}!'.format(outpath))
    except Exception as e:
        print_warning('Conversion of {} was '.format(f) +
                      'not successfull: {}'.format(str(e)))
        break