# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data filtered with xiFDR.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import os

from . import Xi as xi

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def Read(xifdr_file, xi_file, keep=False):
    """
    Collects data from Xi spectrum search filtered by xiFDR and returns an xtable data array.

    Args:
        xi_file: path to percolated Kojak file
        xifdr_file: results file from xiFDR (e.g. PSM_xiFDR)
        keep (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: xtable data table
    """

    print('Reading xiFDR-file: {}'.format(xifdr_file))

    xifdr = pd.read_csv(xifdr_file, delimiter=',')

    xifdr.rename(columns={'run': 'Run',
                          'scan': 'Scan',
                          'Protein1': 'Protein1_FDR',
                          'Protein2': 'Protein2_FDR'}, inplace=True)

    print('Reading xi-file: {}'.format(xi_file))

    xiraw = pd.read_csv(xi_file, delimiter=',')

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

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].apply(\
        lambda row: '-'.join(str(element) for element in row), axis=1)

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        xi.assign_type, axis=1)

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearch and XiFDR'

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    if keep is True:
        # reorder columns to start with the xtable columns
        all_cols = list(xtable.columns.values)
        remaining_cols = [x for x in all_cols if x not in col_order]
        new_order = col_order + remaining_cols

        xtable = xtable[new_order]
    elif keep is False:
        xtable = xtable[col_order]

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
    xi_df = Read(xi_file, xifdr_file, keep=True)

    xi_df.to_excel('test.xls',
                   index=False)

