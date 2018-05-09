# -*- coding: utf-8 -*-

"""
Functions to write xWalk data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd

def Write(xtable, outpath):
    """
    Converts xTable into a list format that can be used as
    input for the xWalk web-server at http://www.xwalk.org/cgi-bin/home.cgi

    Format is:

    residue #  \t  residue #

    As xWalk can only validate one protein at a time, the function
    generates several oouput files for all intra-protein cross-links

    :params: xtable: data table structure
    :params: outpath: path to write file
    """
    # filter inter-protein cross-links
    xtable = xtable[xtable['prot1'] == xtable['prot2']]

    for prot in set(xtable['prot1']):
        thisXtable = xtable[xtable['prot1'] == prot]
        thisXtable[['xpos1', 'xpos2']].drop_duplicates()\
                                      .to_csv('{}_{}.tsv'.format(outpath, alphanum_string(prot)),
                                              index = False,
                                              sep='\t')

def alphanum_string(s):
    """
    Method to clean strings from incorrect characters for file output
    """
    import re
    # new compiler that finds non-alphanumeric characters
    rex = re.compile(r'\W')
    # actually replace the strings
    result = rex.sub('_', s)

    # remove double occurences of _ in the string
    # initialize
    old_char = ''
    new_result = ''
    for char in result:
        if old_char != char:
            new_result += char
        else:
            # prevent removal of double occurences of other strings than _
            if old_char != '_':
                new_result += char
        old_char = char

    return new_result