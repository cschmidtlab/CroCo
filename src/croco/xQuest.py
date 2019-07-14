# -*- coding: utf-8 -*-

"""
Functions to read xQuest data.

"""


import numpy as np
import pandas as pd

import re

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf


def _process_xquest_spectrum(spec_string):
    """
    Extract rawfile name, precursor charge and scan no from xQuest spectrum
    string

    Args:
        spec_string: xQuest spectrum string
    Returns:
        str or np.nan: rawfile name
        int or np.nan: scan number
        int or np.nan: precursor charge
    """
    spectrum_pattern = re.compile('(.+)\.(\d+)\.\d+\..+\.\d+\.\d+\.(\d+)')
    if spectrum_pattern.match(spec_string):
        match = spectrum_pattern.match(spec_string)
        rawfile, scanno, prec_ch = match.groups()
        return str(rawfile), int(scanno), int(prec_ch)
    else:
        return np.nan

def _process_xquest_id(Id_string):
    """
    Extract peptide sequence of the alpha (longer) and the beta (shorter)
    peptide as well as the relative positions of the cross-links within
    these sequences from an xQuest Id-string

    Args:
        ID_string (str): an xQuest Id-String
        type (str): the xlink type from xQuest (monolink, inrtalink, xlink)
    Returns:
        str or np.nan: pepseq1
        str or np.nan: pepseq2
        int or np.nan: xlink1
        int or np.nan: xlink2
    """
    xlink_pattern = re.compile('^(\w+)-(\w+)-a(\d+)-b(\d+)')
    intralink_pattern = re.compile('^(\w+)-\D{1}(\d+)-\D{1}(\d+)')
    monolink_pattern = re.compile('^(\w+)-\D{1}(\d+)-\d+')

    if xlink_pattern.match(Id_string):
        match = xlink_pattern.match(Id_string)
        # pepseq1, pepseq2, xlink1, xlink2
        pepseq1, pepseq2, xlink1, xlink2 = match.groups()

#                return pepseq1, pepseq2, xlink1, xlink2
        return pepseq1, pepseq2, int(xlink1), int(xlink2)
    elif intralink_pattern.match(Id_string):
        match = intralink_pattern.match(Id_string)
        # pepseq1, pepseq2, xlink1, xlink2
        pepseq, xlink1, xlink2 = match.groups()

#                return pepseq, pepseq, xlink1, xlink2
        return pepseq, pepseq, int(xlink1), int(xlink2)
    elif monolink_pattern.match(Id_string):
        match = monolink_pattern.match(Id_string)
        # pepseq1, pepseq2, xlink1, xlink2
        pepseq, xlink = match.groups()

#                return pepseq, np.nan, xlink, np.nan
        return pepseq, np.nan, int(xlink), np.nan
    else:
        return np.nan, np.nan, np.nan, np.nan

def _categorize_xquest_type(XQType):
    """
    Extract protein name and absolute cross-link position from
    xQuest type string (xlink, loop, mono)

    Args:
        XQType (str): xquest type string
    Returns:
        str or np.nan: type of cross-link (inter, loop, mono)
    """

    if XQType == 'xlink':
        return 'inter'
    elif XQType == 'intralink':
        return 'loop'
    elif XQType == 'monolink':
        return 'mono'
    else:
        return np.nan


def Read(xQuest_files, col_order=None, compact=False):
    """
    Read xQuest results file and return file in xTable format.

    Args:
        xQuest_files (list): path to xQuest results file(s)
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to compact the xTable to only those columns listed in col_order

    Returns:
        pandas.DataFrame: xTable data table
    """

    # convert to list if the input is only a single path
    if not isinstance(xQuest_files, list):
        xQuest_files = [xQuest_files]

    allData = list()

    xQuest_dtypes = {'z': pd.Int64Dtype(),
                     'Protein1': str,
                     'Protein2': str,
                     'AbsPos1': pd.Int64Dtype(),
                     'AbsPos2': pd.Int64Dtype(),
                     'ld-Score': float}

    for file in xQuest_files:

        ### Collect data and convert to pandas format
        print('Reading xQuest-file: ' + file)

        # only called if inter_file is not None
#        try:
        s = pd.read_csv(hf.compatible_path(file),
                        delimiter='\t',
                        na_values='-',
                        dtype=xQuest_dtypes)
        allData.append(s)
#        except:
#            raise Exception('[xQuest Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData)

    rename_dict = {'z':'prec_ch',
                   'Protein1':'prot1',
                   'Protein2': 'prot2',
                   'AbsPos1': 'xpos1',
                   'AbsPos2': 'xpos2',
                   'ld-Score': 'score'}

    # Copy and rename selected columns to new xquest df
    try:
        xtable.rename(index=str,
                      columns=rename_dict,
                      inplace=True)
    except Exception as e:
        raise Exception('[xQuest Read] Error during xQuest header renaming: %s' % e)

    # Extract rawfile, scanno and precursor charge from the mgf header string
    # used as Spectrum by xQuest
    xtable[['rawfile', 'scanno', 'prec_ch']] =\
        pd.DataFrame(xtable['Spectrum'].apply(_process_xquest_spectrum).tolist(), index=xtable.index)

    print('[xQuest Read] Processed Spectrum entry')

    # Extract peptide sequences and relative cross-link positions form the
    # xQuest ID-string
    xtable[['pepseq1', 'pepseq2', 'xlink1', 'xlink2']] =\
        pd.DataFrame(xtable['Id'].apply(_process_xquest_id).tolist(), index=xtable.index)

    print('[xQuest Read] Processed xQuest ID' )

    # Modifications are not defined in xQuest
    xtable['mod1'], xtable['mod2'] = "", ""

    # calculate the absolute position of the first amino acide of the resp
    # peptides
    xtable['pos1'] = xtable['xpos1'] - xtable['xlink1'] + 1
    xtable['pos2'] = xtable['xpos2'] - xtable['xlink2'] + 1

    print('[xQuest Read] Calculated positions')

    # Assign mono
    xtable['type'] = xtable['Type'].apply(_categorize_xquest_type)

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for intra and inter xlink to inter/intra/homomultimeric
        intraAndInter = (xtable['type'] == 'inter') | (xtable['type'] == 'intra')
        xtable.loc[intraAndInter, 'type'] =\
            np.vectorize(hf.categorize_inter_peptides)(xtable[intraAndInter]['prot1'],
                                                     xtable[intraAndInter]['pos1'],
                                                     xtable[intraAndInter]['pepseq1'],
                                                     xtable[intraAndInter]['prot2'],
                                                     xtable[intraAndInter]['pos2'],
                                                     xtable[intraAndInter]['pepseq2'])
        print('[xQuest Read] categorized inter peptides')
    else:
        print('[xQuest Read] skipped inter peptide categorization')

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        pd.Series(np.vectorize(hf.generate_id,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                 index=xtable.index).replace('nan', np.nan)

    print('[xQuest Read] Generated ID')

    # xQuest does not incorporate decoy entries in the results table
    # but protein names can contain identifiers as reverse or decoy
    xtable['decoy'] = xtable['ID'].str.contains('reverse') |\
        xtable['ID'].str.contains('decoy')

    # the following properties cannot directly be inferred from the
    # xQuest results file
    # to avoid confusion with missing valued like np.nan, they are set to
    # UNKNOWN
    for header in ['xtype', 'modmass1', 'modpos1', 'modmass2', 'modpos2']:
        xtable[header] = np.nan

    xtable['search_engine'] = 'xQuest'

    xtable = hf.order_columns(xtable, col_order, compact)

    ### Return df
    return xtable

if __name__ == '__main__':
    """
    For testing purposes only
    """

    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    xtable = Read(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\xQuest\20190227_croco_PK_xquest_results_targetdecoy.xls', col_order=col_order)