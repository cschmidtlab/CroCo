# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import os

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def assign_type(row):
    """
    Assign mono, loop, inter and intra link
    based on prot1, prot2, xlink1 and xlink2 entries
    """
    prot1, prot2, xlink1, xlink2 = row
    
    prot1 = str(prot1)
    prot2 = str(prot2)
    xlink1 = str(xlink1)
    xlink2 = str(xlink2)

    if prot2 != 'nan' and prot1 == prot2:
        type = 'intra'
    elif prot2 != 'nan':
        type = 'inter'
    elif prot2 == 'nan' and xlink2 != 'nan':
        type = 'loop'
    elif prot1 != 'nan' and prot2 == 'nan':
        type = 'mono'
    else:
        type = None
    return type

def rawfile_from_source(source_str):
    """
    Exctracts filename from string like
    E:\julian\20180612_croco_testfiles\mgf_msconvert\20180518_JB_jb05a_l100.mgf
    """
    try:
        return source_str.split('.')[-2].split('\\')[-1]
    except AttributeError as e:
        if np.isnan(float(source_str)):
            return np.nan
        else:
            raise Exception(e)
            
def Read(xi_file, keep=False):
    """
    Collects data from Xi spectrum search and returns an xtable data array.

    Args:
        xi_file: path to percolated Kojak file
        keep (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: xtable data table
    """

    print('Reading xi-file: {}'.format(xi_file))

    data = pd.read_csv(xi_file, delimiter=',')
    
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
    
    xtable['rawfile'] = xtable['Source'].apply(rawfile_from_source)

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].apply(\
        lambda row: '-'.join(str(element) for element in row), axis=1)

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        assign_type, axis=1)    
    
    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearch'

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
    
    xi = Read(xi_file)
    
    xi.to_excel('test.xls',
                 index=False)
    
    