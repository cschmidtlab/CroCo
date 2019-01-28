# -*- coding: utf-8 -*-

"""
Functions to read and process data generated with the Kojak cross-link
search engine.
"""

import numpy as np
import pandas as pd

import os

if __name__ == '__main__':
    import HelperFunctions as hf
    import KojakFunctions as kj
else:
    from . import HelperFunctions as hf
    from . import KojakFunctions as kj

def init(this_order):
    """
    Initialises the column order when called from the GUI.
    No function if calling directly.
    """
    global col_order
    col_order = this_order

def Read(kojak_files, rawfile=None, compact=False, decoy_string='REVERSE'):
    """
    Read Kojak results file, calculate and process missing values required
    for xTable and return the xTable.

    Args:
        kojak_file (str): path or paths to Kojak results file(s)
        rawfile (str): name of the corresponding rawfile
        decoy_string (optional): string used in kojak to label decoys

    Returns:
        pandas.DataFrame: xtable data table
    """

    # convert to list if the input is only a single path
    if not isinstance(kojak_files, list):
        kojak_files = [kojak_files]

    allData = list()

    for file in kojak_files:


        print('Reading Kojak-file: ' + kojak_file)

        # only called if kojak_file is not None
        try:
            s = pd.read_csv(hf.FSCompatiblePath(kojak_file),
                                 skiprows = 1, # skip the Kojak version
                                 delimiter='\t')
            allData.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData)

    # remove lines containing non-identified PSMs (marked with '-' in both
    # Link columns
    xtable = xtable[(xtable['Link #1'] != '-') & (xtable['Link #2'] != '-')]

    # dropping lines causes fragmented index --> regenate the index
    xtable.reset_index(drop=True, inplace=True)

    # if split into mulitple rows if multiple candidate proteins were found to
    # match an experimental spectrum
    xtable = hf.split_concatenated_lists(xtable, where=['Protein #1', 'Protein #2'])

    xtable[['scanno', 'prec_ch', 'xlink1', 'xlink2', 'score']] =\
        xtable[['Scan Number', 'Charge', 'Link #1', 'Link #2', 'Score']].astype(int)

    # Extract peptide sequence, modification mass and position from the
    # Peptide #1 and Peptide #2 entries
    xtable = kj.extract_peptide(xtable)

    # transform unset xlinks to np.nan
    xtable[['xlink1', 'xlink2']] = xtable[['xlink1', 'xlink2']].replace(-1, np.nan)

    # extract protein name and relative cross-link position from the Protein #
    # entries
    xtable = kj.extract_protein(xtable)

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    # Must be calculated as float as NaN is not implemented in int
    xtable['pos1'] =\
        xtable['xpos1'].astype(float, errors='ignore') - \
        xtable['xlink1'].astype(float, errors='ignore') + 1
    xtable['pos2'] =\
        xtable['xpos2'].astype(float, errors='ignore') - \
        xtable['xlink2'].astype(float, errors='ignore') + 1

    # Calculate if a cross link is of inter or of loop type
    # Refine the inter type into inter/intra/homomultimeric
    # Generate ID for the xlinks
    xtable = kj.assign_ID_and_type(xtable)

    #sets the column decoy based on whether the decoy string is present in the
    # protein name or not
    xtable = kj.set_decoy(xtable, decoy_string)

    # set the rawfile name for xtable (None if not provided by call)
    xtable['rawfile'] = rawfile

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'Kojak'

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    xtable = hf.applyColOrder(xtable, col_order, compact)

    ### return xtable df
    return xtable

if __name__ == '__main__':
    kojak_file = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\kojak\20180518_JB_jb05a_l50.kojak.txt'

    col_order = ['rawfile', 'scanno', 'prec_ch',
                 'pepseq1', 'xlink1',
                 'pepseq2', 'xlink2', 'xtype',
                 'modmass1', 'modpos1', 'mod1',
                 'modmass2', 'modpos2', 'mod2',
                 'prot1', 'xpos1', 'prot2',
                 'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    init(col_order)

    xtable = Read(kojak_file, rawfile='Test')
