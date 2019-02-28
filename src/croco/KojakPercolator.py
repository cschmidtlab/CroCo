# -*- coding: utf-8 -*-

"""
Functions to read Percolator processed Kojak data.
"""

import numpy as np
import pandas as pd

if __name__ == '__main__':
    import HelperFunctions as hf
    import KojakFunctions as kj
else:
    from . import HelperFunctions as hf
    from . import KojakFunctions as kj

def Read(perc_files, rawfile=None, percolator_string='.validated', decoy_string='REVERSE', compact=False, col_order=None):
    """
    Collects unprocessed and percolated results and returns an xtable data array.

    Args:
        perc_file (str): path or list of paths to percolated Kojak file(s)
        percolator_string (str): user-defined string appended to the percolated filenames
        decoy_string (optional): string used in kojak to label decoys
        rawfile (str): name of the corresponding rawfile
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to compact the xTable to only those columns listed in col_order

    Returns:
        pandas.DataFrame: xtable data table
    """
    # convert to list if the input is only a single path
    if not isinstance(perc_files, list):
        perc_files = [perc_files]
    
    allData = list()

    kojak_dtypes = {'Scan Number': pd.Int64Dtype(),
                    'Charge': pd.Int64Dtype(),
                    'Link #1': pd.Int64Dtype(),
                    'Link #2': pd.Int64Dtype(),
                    'Score': float
                    }

    for p_file in perc_files:
        ### Collect data and convert to pandas format
    
        print('Reading Percolator-file: ' + p_file)
    
        try:
            percolated = pd.read_csv(hf.FSCompatiblePath(p_file),
                                     delimiter='\t',
                                     usecols=range(5),
                                     index_col=False, # avoid taking the first col as index
                                     engine='python')
        except FileNotFoundError:
            raise Exception("Could not find the percolated file %s. Please move it into the same directory as the percolator files!" % p_file)
    
        percolated.rename(columns={'PSMId': 'SpecId'}, inplace=True)
    
        unperc_file = p_file.replace(percolator_string, '')
    
        print('Reading Percolator input: ' + unperc_file)
    
        try:
            unpercolated = pd.read_csv(hf.FSCompatiblePath(unperc_file),
                                      delimiter = '\t',
                                      usecols=range(10),
                                      engine='python',
                                      index_col=False)
        except FileNotFoundError:
            raise Exception("Could not find the unpercolated file %s. Please move it into the same directory as the percolator files!" % unperc_file)
    
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
                                dtype=kojak_dtypes,
                                na_values='-',
                                delimiter='\t')
        except FileNotFoundError:
            raise Exception("Could not find the kojak_file %s. Please move it into the same directory as the percolator files!" % kojak_file)
    
        kojak.rename(columns={'Scan Number': 'scannr'}, inplace=True)
    
        s = pd.merge(xtable, kojak, on=['scannr', 'Charge', 'dScore', 'Score'], how='left')
        
        allData.append(s)

    xtable = pd.concat(allData, sort=False)

    # split ambiguous concatenated protein names
    xtable = hf.split_concatenated_lists(xtable, where=['Protein #1', 'Protein #2'])

    ### Process the data to comply to xTable format
    xtable = xtable.rename(columns={'Scan Number': 'scanno',
                                    'Charge': 'prec_ch',
                                    'Link #1': 'xlink1',
                                    'Link #2': 'xlink2',
                                    'Score': 'score'
                                    })

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

    xtable['search_engine'] = 'Kojak and Percolator'

    xtable = hf.applyColOrder(xtable, col_order, compact)
    
    return xtable

if __name__ == '__main__':
    # defines the column headers required for xtable output
    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    perc_file = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\20190226_croco_testfiles_kojak\20180615_KS_CL_9_msconvert.perc.intra.validated.txt'

    xtable = Read(perc_file)