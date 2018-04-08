#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
The CroCo Cross-Link Converter console

Convert results from data analysis of chemical cross-linking /
mass-spectrometry experiments.
"""

import sys, os

#######################
# Initialise Paths
######################

croco_script = sys.argv[0]
# check if programme was called via symlink
if os.path.islink(croco_script):
    croco_script = os.readlink(croco_script)
# dir is the directory above the bin-dir
croco_dir = os.path.abspath(os.path.join(os.path.dirname(croco_script), '..'))

sys.path.append(os.path.abspath(os.path.join('..', croco_dir)))

import argparse
import pandas as pd

import croco

description = """The CroCo cross-link converter:
-------------------------------
Convert results from data analysis of chemical cross-linking mass-spectrometry experiments."""

epilog = r"""EXAMPLES:
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

in_dict = {'pLink1': croco.pLink1.Read,
           'pLink2': croco.pLink2.Read,
           'Kojak': croco.Kojak.Read,
           'xTable': pd.read_csv}

out_dict = {'xTable': croco.xTable.Write,
            'xVis': croco.xVis.Write,
            'xiNet': croco.xiNET.Write,
            'DynamXL': croco.DynamXL.Write,
            'xWalk': croco.xWalk.Write}
                        
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