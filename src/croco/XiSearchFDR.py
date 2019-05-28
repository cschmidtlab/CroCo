# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data filtered with xiFDR.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import os

if __name__ == '__main__':
    import HelperFunctions as hf
    import Xi as xi
else:
    from . import Xi as xi
    from . import HelperFunctions as hf

def _assign_type(row):
    """
    Assign mono, loop, inter and intra link
    based on prot1, prot2, xlink1 and xlink2 entries

    Args:
        row (Series): a series or list containing prot1, prot2, xlink1, xlink2
    Returns:
        str or np.nan: type of cross-link (inter, intra, loop, mono)
    """
    prot1, prot2, xlink1, xlink2 = row

    prot1 = str(prot1)
    prot2 = str(prot2)
    xlink1 = str(xlink1)
    xlink2 = str(xlink2)

    if prot2 != 'nan' and prot1 == prot2:
        t = 'intra'
    elif prot2 != 'nan':
        t = 'inter'
    elif prot2 == 'nan' and xlink2 != 'nan':
        t = 'loop'
    elif prot1 != 'nan' and prot2 == 'nan' and xlink1 != 'nan':
        t = 'mono'
    else:
        t = np.nan
    return t

def _modifications_from_sequence(sequence):
    """
    Extract a modification name and its position from a sequence containing
    the modifications as lowercase characters
    
    Args:
        sequence (str): a sequence to be parsed
    Returns:
        str: sequence without the modification characters
        list of str: modification names
        list of int: modification positions within the peptide
    """
    mods = []
    modposns = []
    
    cur_pos = 0
    this_mod = ''
    clean_sequence = ''
    
    for char in sequence:
        if char.isupper():
            if this_mod != '':
                mods.append(this_mod)
                modposns.append(cur_pos)
            this_mod = ''
            clean_sequence += char
            cur_pos += 1
        else:
            this_mod += char
            
    return clean_sequence, mods, modposns

def Read(xifdr_files, modstring=None, col_order=None, compact=False):
    """
    Collects data from Xi spectrum search filtered by xiFDR and returns an xtable data array.

    Args:
        xifdr_files: path or list of paths to xiFDR file(s)
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: xtable data table
    """
    # convert to list if the input is only a single path
    if not isinstance(xifdr_files, list):
        xifdr_files = [xifdr_files]

    allData = list()

    xifdr_dtypes = {'scan': pd.Int64Dtype(),
                   'exp charge': pd.Int64Dtype(),
                   'PepSeq1': str,
                   'PepSeq2': str,
                   'LinkPos1': pd.Int16Dtype(),
                   'LinkPos2': pd.Int16Dtype(),
                   'Protein1': str,
                   'Protein2': str,
                   'ProteinLinkPos1': pd.Int16Dtype(),
                   'ProteinLinkPos2': pd.Int16Dtype(),
                   'PepPos1': pd.Int16Dtype(),
                   'PepPos2': pd.Int16Dtype(),
                   }

    for file in xifdr_files:

        print('Reading xiFDR-file: {}'.format(file))
        try:
            s = pd.read_csv(hf.FSCompatiblePath(file), delimiter=',', dtype=xifdr_dtypes)
            allData.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData)

    ### Process the data to comply to xTable format
    xtable = xtable.rename(columns={'scan': 'scanno',
                                    'exp charge': 'prec_ch',
                                    'PepSeq1': 'pepseq1',
                                    'LinkPos1': 'xlink1',
                                    'PepSeq2': 'pepseq2',
                                    'LinkPos2': 'xlink2',
                                    # modmass1
                                    # modmass2
                                    'Protein1': 'prot1',
                                    'ProteinLinkPos1': 'xpos1',
                                    'Protein2': 'prot2',
                                    'ProteinLinkPos2': 'xpos2',
                                    'PepPos1': 'pos1',
                                    'PepPos2': 'pos2',
                                    })

    # Extract clean sequence and modificiations from the sequence string
    xtable[['pepseq1', 'mod1', 'modpos1']] =\
        pd.DataFrame(xtable['pepseq1'].apply(_modifications_from_sequence).tolist(), index=xtable.index)
    xtable[['pepseq2', 'mod2', 'modpos2']] =\
        pd.DataFrame(xtable['pepseq2'].apply(_modifications_from_sequence).tolist(), index=xtable.index)

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
          _assign_type, axis=1)
    
    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for inter xlink to inter/intra/homomultimeric
        onlyInter = xtable['type'] == 'inter'
        xtable.loc[onlyInter, 'type'] =\
            np.vectorize(hf.categorizeInterPeptides)(xtable[onlyInter]['prot1'],
                                                     xtable[onlyInter]['pos1'],
                                                     xtable[onlyInter]['pepseq1'],
                                                     xtable[onlyInter]['prot2'],
                                                     xtable[onlyInter]['pos2'],
                                                     xtable[onlyInter]['pepseq1'])
        print('[xiFDR Read] categorized inter peptides')
    else:
        print('[xiFDR Read] skipped inter peptide categorization')    
 
    
    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        pd.Series(np.vectorize(hf.generateID,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                index=xtable.index).replace('nan', np.nan)

    xtable['decoy'] = xtable['Decoy1'] | xtable['Decoy2']

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearchFDR'

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

    os.chdir(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\xi\pXtract_msconvertStyle')
    xifdr_files = r'XiFDR_5_FDR_PSM_PSM_xiFDR1.0.22.csv'
    xtable = Read(xifdr_files, compact=True)
