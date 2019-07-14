# -*- coding: utf-8 -*-

"""
Functions to read and process data generated with the Kojak cross-link
search engine.
"""

import numpy as np
import pandas as pd

if __name__ == '__main__':
    import HelperFunctions as hf
    import KojakFunctions as kj
else:
    from . import HelperFunctions as hf
    from . import KojakFunctions as kj

def Read(kojak_files, rawfile=None, decoy_string='decoy', col_order=None, compact=False):
    """
    Read Kojak results file, calculate and process missing values required
    for xTable and return the xTable.

    Args:
        kojak_files (list): path or paths to Kojak results file(s)
        rawfile (str): name of the corresponding rawfile
        decoy_string (optional): string used in kojak to label decoys
        col_order (list) – List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Compact the xTable to only the columns given in col_order or not

    Returns:
        pandas.DataFrame: xtable data table
    """

    # convert to list if the input is only a single path
    if not isinstance(kojak_files, list):
        kojak_files = [kojak_files]

    allData = list()

    kojak_dtypes = {'Scan Number': pd.Int64Dtype(),
                    'Charge': pd.Int64Dtype(),
                    'Link #1': pd.Int64Dtype(),
                    'Link #2': pd.Int64Dtype(),
                    'Score': float
                    }

    for file in kojak_files:


        print('Reading Kojak-file: ' + file)

        # only called if kojak_file is not None
        try:
            s = pd.read_csv(hf.compatible_path(file),
                            skiprows = 1, # skip the Kojak version
                            dtype=kojak_dtypes,
                            na_values = '-',
                            delimiter='\t')
            allData.append(s)
        except Exception as e:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData)

    # remove lines containing non-identified PSMs (marked with '-' in both
    # Link columns
    xtable.dropna(axis=0, how='all', subset=['Link #1', 'Link #2'], inplace=True)

    # dropping lines causes fragmented index --> regenate the index
    xtable.reset_index(drop=True, inplace=True)

    # if split into mulitple rows if multiple candidate proteins were found to
    # match an experimental spectrum
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

    xtable['search_engine'] = 'Kojak'

    xtable = hf.order_columns(xtable, col_order, compact)

    ### return xtable df
    return xtable

if __name__ == '__main__':
    kojak_file = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\final\input\kojak\20180615_KS_CL_9_msconvert.kojak.txt'

    col_order = ['rawfile', 'scanno', 'prec_ch',
                 'pepseq1', 'xlink1',
                 'pepseq2', 'xlink2', 'xtype',
                 'modmass1', 'modpos1', 'mod1',
                 'modmass2', 'modpos2', 'mod2',
                 'prot1', 'xpos1', 'prot2',
                 'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    xtable = Read(kojak_file, col_order=col_order, rawfile='20180615_KS_CL_9')
