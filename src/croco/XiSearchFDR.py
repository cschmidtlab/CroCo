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

def Read(xifdr_file, xi_file, compact=False):
    """
    Collects data from Xi spectrum search filtered by xiFDR and returns an xtable data array.

    Args:
        xifdr_file: results file from xiFDR (contains PSM_xiFDR)
        xi_file: path to percolated Kojak file
        keep (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: xtable data table
    """
    
    if isinstance(xifdr_file, list):
        if len(xifdr_file) > 1:
            raise Exception('[xiFDR Read] Sorry! Only one xiFDR file per conversion is allowed to unambiguously relate it to a xi-file')
        xifdr_file = xifdr_file[0]

    if not 'PSM_xiFDR' in xifdr_file:
        raise Exception('[xiFDR Read] The string "PSM_xiFDR" is missing in your input file. Did you choose the right file?')

    print('Reading xiFDR-file: {}'.format(xifdr_file))

    try:
        xifdr = pd.read_csv(hf.FSCompatiblePath(xifdr_file), delimiter=',')
    
        xifdr.rename(columns={'run': 'Run',
                              'scan': 'Scan',
                              'Protein1': 'Protein1_FDR',
                              'Protein2': 'Protein2_FDR'}, inplace=True)
    except Exception as e:
        raise Exception('[xiFDR Read] Error while reading the xiFDR file: {}'.format(e))

    print('Reading xi-file: {}'.format(xi_file))

    xiraw = pd.read_csv(hf.FSCompatiblePath(xi_file), delimiter=',')

    # Merge with left join (only keys that are in tje percolated DF will be re-
    # tained)
    data = pd.merge(xifdr, xiraw, on=['Run', 'Scan'], how='left')

    ### Process the data to comply to xTable format
    xtable = data.rename(columns={'Scan': 'scanno',
                                  'PrecoursorCharge': 'prec_ch',
                                  'BasePeptide1': 'pepseq1',
                                  'ProteinLink1': 'xpos1',
                                  'BasePeptide2': 'pepseq2',
                                  'ProteinLink2': 'xpos2',
                                  'ModificationMasses1': 'modmass1',
                                  'ModificationMasses2': 'modmass2',
                                  'Modifications1': 'mod1',
                                  'Modifications2': 'mod2',
                                  'Protein1': 'prot1',
                                  'Protein2': 'prot2',
                                  'Start1': 'pos1',
                                  'Start2': 'pos2',
                                  'Link1': 'xlink1',
                                  'Link2': 'xlink2',
                                  'ModificationPositions1': 'modpos1',
                                  'ModificationPositions2': 'modpos2',
                                  'match score': 'score'
                                  })

    xtable['rawfile'] = xtable['Source'].apply(xi.rawfile_from_source)

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        xi.assign_type, axis=1)

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

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearch and XiFDR'

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')
    
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

    os.chdir(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Xi')
    xi_file = r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Xi\20180612_croco_testfiles_XiVersion1.6.739.csv'
    xifdr_file = r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Xi\20180612_croco_testfiles_5_FDR_PSM_xiFDR1.0.22.csv'
    xi_df = Read(xi_file, xifdr_file, compact=True)

    xi_df.to_excel('test.xls',
                   index=False)

