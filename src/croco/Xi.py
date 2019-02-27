# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data.

"""

import numpy as np
import pandas as pd

import os

if __name__ in ['__main__', 'Xi']:
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf


def assign_type(row):
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
        t = None
    return t

def rawfile_from_source(source_str):
    """
    Exctracts filename from string like
    E:\julian\20180612_croco_testfiles\mgf_msconvert\20180518_JB_jb05a_l100.mgf

    Args:
        source_str (str): Path to a rawfile
    Returns:
        str: filename from path
    """
    try:
        return source_str.split('.')[-2].split('\\')[-1]
    except AttributeError as e:
        if np.isnan(float(source_str)):
            return np.nan
        else:
            raise Exception(e)

def Read(xi_files, col_order=None, compact=False):
    """
    Collects data from Xi spectrum search and returns an xtable data array.

    Args:
        xi_file: path or list of paths to xi file(s)
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to compact the xTable to only those columns listed in col_order
    Returns:
        pandas.DataFrame: xtable data table
    """

    # convert to list if the input is only a single path
    if not isinstance(xi_files, list):
        xi_files = [xi_files]

    allData = list()

    for file in xi_files:

        print('Reading xi-file: {}'.format(file))
        try:
            s = pd.read_csv(hf.FSCompatiblePath(file), delimiter=',')
            allData.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData)

    ### Process the data to comply to xTable format
    xtable = xtable.rename(columns={'Scan': 'scanno',
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

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        assign_type, axis=1)

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        np.vectorize(hf.generateID)(xtable['type'], xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'])

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for inter xlink to inter/intra/homomultimeric
        xtable.loc[xtable['type'] == 'inter', 'type'] =\
            np.vectorize(hf.categorizeInterPeptides)(xtable[xtable['type'] == 'inter']['prot1'],
                                                     xtable[xtable['type'] == 'inter']['pos1'],
                                                     xtable[xtable['type'] == 'inter']['pepseq1'],
                                                     xtable[xtable['type'] == 'inter']['prot2'],
                                                     xtable[xtable['type'] == 'inter']['pos2'],
                                                     xtable[xtable['type'] == 'inter']['pepseq2'])
        print('[xi Read] categorized inter peptides')
    else:
        print('[xi Read] skipped inter peptide categorization')

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearch'

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

    xi = Read(xi_file)

    xi.to_excel('test.xls',
                 index=False)

