# -*- coding: utf-8 -*-

"""
Functions to read Percolator prcoesse Kojak data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

if __name__ == '__main__':
    import HelperFunctions as hf
    import KojakFunctions as kj
else:
    from . import HelperFunctions as hf
    from . import KojakFunctions as kj

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def calc_pos_from_xpos(xpos, xlink):
    """
    Calculates the absolute position of the first AA of a peptide
    sequence from the absolute position of the cross-link AA (xpos) its
    position within the sequence (xlink)

    Returns: pos - Absolute position of AA in sequence
    """
    try:
        if np.isnan(float(xpos)):
            return np.nan

        xpos = int(xpos)
        xlink = int(xlink)

        return xpos - xlink + 1

    except Exception as e:
        print('{}: xpos was {} and xlink was {}'.format(e, xpos, xlink))
        return np.nan

def Read(perc_files, percolator_string='.validated', decoy_string='REVERSE', rawfile=None, compact=False):
    """
    Collects unprocessed and percolated results and returns an xtable data array.

    Args:
        perc_file: path or list of paths to percolated Kojak file(s)
        percolator_string: user-defined string appended to the percolated filenames
        rawfile: name of the corresponding rawfile

    Returns:
        xtable: xtable data table
    """
    # convert to list if the input is only a single path
    if not isinstance(perc_files, list):
        perc_files = [perc_files]
    
    allData = list()
    
    for p_file in perc_files:
        ### Collect data and convert to pandas format
    
        print('Reading Percolator-file: ' + p_file)
    
        # only called if inter_file is not None
    
        percolated = pd.read_csv(hf.FSCompatiblePath(p_file),
                                 delimiter='\t',
                                 usecols=range(5),
                                 index_col=False, # avoid taking the first col as index
                                 engine='python')
    
        percolated.rename(columns={'PSMId': 'SpecId'}, inplace=True)
    
        unperc_file = p_file.replace(percolator_string, '')
    
        print('Reading Percolator input: ' + unperc_file)
    
        try:
            unpercolated = pd.read_csv(hf.FSCompatiblePath(unperc_file),
                                      delimiter = '\t',
                                      usecols=range(10),
                                      engine='python',
                                      index_col=False)
        except:
            raise FileNotFoundError(unperc_file)
    
        # Merge with left join (only keys that are in tje percolated DF will be re-
        # tained)
        xtable = pd.merge(percolated, unpercolated, on='SpecId', how='left')
    
        # Reading the Kojak-file is required to get additional information on the
        # matches such as the corresponding protein names
        kojak_file = unperc_file[0:unperc_file.find('.perc')] + '.kojak.txt'
    
        print('Reading Kojak-file: ' + kojak_file)
    
        try:
            kojak = pd.read_csv(hf.FSCompatiblePath(kojak_file),
                                skiprows = 1, # skip the Kojak version
                                delimiter='\t')
        except:
            raise FileNotFoundError("Could not find the kojak_file %s. Please move it into the same directory as the percolator files!" % kojak_file)
    
        kojak.rename(columns={'Scan Number': 'scannr'}, inplace=True)
    
        s = pd.merge(xtable, kojak, on=['scannr', 'Charge', 'dScore', 'Score'], how='left')
        
        allData.append(s)

    xtable = pd.concat(allData)

    # split ambiguous concatenated protein names
    xtable = hf.split_concatenated_lists(xtable, where=['Protein #1', 'Protein #2'])

    xtable[['scanno', 'prec_ch', 'xlink1', 'xlink2', 'score']] =\
        xtable[['scannr', 'Charge', 'Link #1', 'Link #2', 'Score']].astype(int)

    # Extract peptide sequence, modification mass and position from the
    # Peptide #1 and Peptide #2 entries
    xtable = kj.extract_peptide(xtable)

    # transform unset xlinks to np.nan
    xtable[['xlink1', 'xlink2']] = xtable[['xlink1', 'xlink2']].replace('-1', np.nan)

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

    xtable['search_engine'] = 'Kojak and Percolator'

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    xtable = hf.applyColOrder(xtable, col_order, compact)

    return xtable

if __name__ == '__main__':
    import os

    # defines the column headers required for xtable output
    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    os.chdir(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Kojak')
    perc_file = r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Kojak\20180518_JB_jb05a_l50.perc.loop.validated.txt'

    perc = Read(perc_file)

#    perc.to_excel('test.xls',
#                  index=False)