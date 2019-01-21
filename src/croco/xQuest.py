# -*- coding: utf-8 -*-

"""
Functions to read Kojak data.

This script is part of the CroCo cross-link converter project
"""


import numpy as np
import pandas as pd

import re

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def process_xQuest_spectrum(spec_string):
    """
    Extract rawfile name, precursor charge and scan no from xQuest sequence
    string

    :returns: list of rawfile, scanno, prec_ch
    """
    spectrum_pattern = re.compile('(.+)\.(\d+)\.\d+\..+\.\d+\.\d+\.(\d+)')
    if spectrum_pattern.match(spec_string):
        match = spectrum_pattern.match(spec_string)
        # rawfile, scanno, prec_ch
        return match.groups()
    else:
        return np.nan

def process_xQuest_Id(Id_string):
    """
    Extract peptide sequence of the alpha (longer) and the beta (shorter)
    peptide as well as the relative positions of the cross-links within
    these sequences from an xQuest Id-string
    """
    Id_pattern = re.compile('(\w+)-(\w+)-a(\d+)-b(\d+)')
    if Id_pattern.match(Id_string):
        match = Id_pattern.match(Id_string)
        # pepseq1, pepseq2, xlink1, xlink2
        return match.groups()
    else:
        return np.nan

def Read(xquest_file, compact=False):
    """
    Read xQuest results file and return file in xTable format.

    :params: kojakpath: Kojak results file
    :params: rawfile: name of the corresponding rawfile

    :returns: xtable data table
    """

    ### Collect data and convert to pandas format

    print('Reading xQuest-file: ' + xquest_file)

    # only called if inter_file is not None
    try:
        xquest = pd.read_csv(xquest_file,
                             delimiter='\t')
    except:
        raise(FileNotFoundError('xQuest results file not found'))

    rename_dict = {'z':'prec_ch',
                   'Protein1':'prot1',
                   'Protein2': 'prot2',
                   'AbsPos1': 'xpos1',
                   'AbsPos2': 'xpos2',
                   'Type': 'type',
                   'ld-Score': 'score'}

    # Copy and rename selected columns to new xquest df
    try:
        xtable = xquest.loc[:, list(rename_dict.keys())]
        xtable.rename(index=str,
                    columns=rename_dict,
                    inplace=True)
    except Exception as e:
        raise Exception('Error during xQuest header renaming: %s' % e)

    # Extract rawfile, scanno and precursor charge from the mgf header string
    # used as Spectrum by xQuest
    xtable['rawfile'], xtable['scanno'], xtable['prec_ch'] =\
        zip(*xquest['Spectrum'].apply(process_xQuest_spectrum))

    # Extract peptide sequences and relative cross-link positions form the
    # xQuest ID-string
    xtable['pepseq1'], xtable['pepseq2'], xtable['xlink1'], xtable['xlink2'] =\
        zip(*xquest['Id'].apply(process_xQuest_Id))

    # Modifications are not defined in xQuest
    xtable['mod1'], xtable['mod2'] = "", ""

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        np.vectorize(hf.generateID)(xtable['type'], xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'])

    # Reassign the type for inter xlink to inter/intra/homomultimeric
    xtable.loc[xtable['type'] == 'inter', 'type'] =\
        np.vectorize(hf.categorizeInterPeptides)(xtable[xtable['type'] == 'inter']['prot1'],
                                                 xtable[xtable['type'] == 'inter']['pos1'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'],
                                                 xtable[xtable['type'] == 'inter']['prot2'],
                                                 xtable[xtable['type'] == 'inter']['pos2'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'])

    # xQuest does not incorporate decoy entries in the resutls table
    # but protein names can contain identifiers as reverse or decoy
    xtable['decoy'] = xtable['ID'].str.contains('reverse') |\
        xtable['ID'].str.contains('decoy')

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')
    
    xtable = hf.applyColOrder(xtable, col_order, compact)

    ### Return df
    return xtable