# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data filtered with xiFDR.

This script is part of the CroCo cross-link converter project.
Last tested with Xi version 1.8.7 and xiFDR version 2.3.5.
"""

import numpy as np
import pandas as pd
import os, re

if __name__ == '__main__':
    import HelperFunctions as Help
    from Xi import XiConfig
else:
    from . import HelperFunctions as Help
    from .Xi import XiConfig


def _modifications_from_sequence(sequence, moddict):
    """
    Extract a modification name and its position from a sequence containing
    the modifications as lowercase characters
    
    Args:
        sequence (str): a sequence to be parsed
        moddict (dict): a dictionary mapping symbols for modified amino acids to tuples of corresponding amino acids and their masses
    Returns:
        str: sequence without the modification characters
        list of str: modification names
        list of int: modification positions within the peptide
    """
    mods = []
    modposns = []
    modmasses = []

    # list containign the symbal, start position, end position of the
    # symbol in the original sequence
    found = list()
    for symbol in moddict.keys():
        if symbol in sequence:
            for m in re.finditer(symbol, sequence):
                found.append((symbol, m.start(), m.end()))
            sequence = re.sub(symbol, moddict[symbol][0], sequence)

    # sort in place
    found.sort(key=lambda x: x[1])

    bias = 0
    for match in found:
        mods.append(match[0])
        modmasses.append(moddict[match[0]][1])
        modposns.append(1 + match[1] - bias)
        bias += len(match[0]) - len(moddict[match[0]][0])

    return sequence, mods, modposns, modmasses


def _mods_from_xi_config(xi_config):
    """
    Extract a modifications dictionary from a xi config file
    
    Args:
        xi_config (str): path to xi_config file
    
    Returns:
        dict: dictionary mapping modification symbols to a list of their unmodified counterpart and the modified mass
    """

    config_instance = XiConfig(xi_config)
    conf = config_instance.config

    moddict = dict()

    for k, v in conf.items():
        if 'modification' in k:
            for x in v:
                for aa in x['MODIFIED'].split(','):
                    symbol = aa + x['SYMBOLEXT']
                    if symbol not in moddict:
                        moddict[symbol] = aa, float(x['DELTAMASS'])
        elif 'loss' in k:
            for x in v:
                if {'NAME', 'AMINOACIDS', 'MASS'}.issubset(x):
                    for aa in x['AMINOACIDS'].split(','):
                        symbol = aa + x['NAME']
                        if symbol not in moddict:
                            moddict[symbol] = aa, -float(x['MASS'])
        elif 'crosslinker' in k:
            for x in v:
                if {'Name', 'LINKEDAMINOACIDS', 'MASS'}.issubset(x):
                    for aa in x['LINKEDAMINOACIDS'].split(','):
                        aa = aa[0]  # the string also contains the link probability
                        # crosslinker names are not unique, so we need to add the name to the symbol
                        # to make it unique
                        symbol = aa + x['Name']
                        if symbol not in moddict:
                            moddict[symbol] = aa, float(x['MASS'])

    return moddict


def Read(xifdr_files, xi_config, col_order=Help.default_col_order, compact=False):
    """
    Collects data from Xi spectrum search filtered by xiFDR and returns a xtable data array.

    Args:
        xifdr_files: path or list of paths to xiFDR file(s)
        xi_config: path to corresponding xi_config file
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: xtable data table
    """
    # convert to list if the input is only a single path
    if not isinstance(xifdr_files, list):
        xifdr_files = [xifdr_files]

    all_data = list()

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
            s = pd.read_csv(Help.compatible_path(file), delimiter=',', dtype=xifdr_dtypes)
            all_data.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(all_data)

    ### Process the data to comply to xTable format

    xtable = xtable.rename(columns={'exp charge': 'prec_ch',
                                    'LinkPos1': 'xlink1',
                                    'LinkPos2': 'xlink2',
                                    'Protein1': 'prot1',
                                    'ProteinLinkPos1': 'xpos1',
                                    'Protein2': 'prot2',
                                    'ProteinLinkPos2': 'xpos2',
                                    'PepPos1': 'pos1',
                                    'PepPos2': 'pos2',
                                    'Score': 'score',
                                    'run': 'rawfile',
                                    'scan': 'scanno',
                                    })

    moddict = _mods_from_xi_config(xi_config)

    # Extract clean sequence and modifications from the sequence string
    xtable[['pepseq1', 'mod1', 'modpos1', 'modmass1']] = \
        pd.DataFrame(xtable['PepSeq1'].apply(lambda x: _modifications_from_sequence(x, moddict)).tolist(),
                     index=xtable.index)
    xtable[['pepseq2', 'mod2', 'modpos2', 'modmass2']] = \
        pd.DataFrame(xtable['PepSeq2'].apply(lambda x: _modifications_from_sequence(x, moddict)).tolist(),
                     index=xtable.index)

    # assign categories of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(
        Help.assign_type, axis=1)

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for inter xlink to inter/intra/homomultimeric
        only_inter = xtable['type'] == 'inter'
        xtable.loc[only_inter, 'type'] = \
            np.vectorize(Help.categorize_inter_peptides)(xtable[only_inter]['prot1'],
                                                         xtable[only_inter]['pos1'],
                                                         xtable[only_inter]['pepseq1'],
                                                         xtable[only_inter]['prot2'],
                                                         xtable[only_inter]['pos2'],
                                                         xtable[only_inter]['pepseq1'])
        print('[xiFDR Read] categorized inter peptides')
    else:
        print('[xiFDR Read] skipped inter peptide categorization')

        # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] = \
        pd.Series(np.vectorize(Help.generate_id,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                  index=xtable.index).replace('nan', np.nan)

    xtable['decoy'] = xtable['Decoy1'] | xtable['Decoy2']

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearchFDR'

    xtable = Help.order_columns(xtable, col_order, compact)

    return xtable


if __name__ == '__main__':
    # defines the column headers required for xtable output
    col_order = ['rawfile', 'scanno', 'prec_ch',
                 'pepseq1', 'xlink1',
                 'pepseq2', 'xlink2', 'xtype',
                 'modmass1', 'modpos1', 'mod1',
                 'modmass2', 'modpos2', 'mod2',
                 'prot1', 'xpos1', 'prot2',
                 'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    os.chdir(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\final\input\xi_xiFDR')
    xifdr_files = r'XiFDR_5_FDR_PSM_PSM_xiFDR1.0.22.csv'
    xi_config = r'xi_config.conf'
    xtable = Read(xifdr_files, xi_config, compact=False)
