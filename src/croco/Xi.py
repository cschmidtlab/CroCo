# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data.

"""

import numpy as np
import pandas as pd
import re

if __name__ in ['__main__', 'Xi']:
    import HelperFunctions as Help
else:
    from . import HelperFunctions as Help


def _rawfile_from_source(source_str):
    r"""
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

    xi_dtypes = {'Scan': pd.Int64Dtype(),
                 'PrecoursorCharge': pd.Int64Dtype(),
                 'BasePeptide1': str,
                 'ProteinLink1': pd.Int16Dtype(),
                 'BasePeptide2': str,
                 'ProteinLink2': pd.Int16Dtype(),
                 'Protein1': str,
                 'Protein2': str,
                 'Start1': pd.Int32Dtype(),
                 'Start2': pd.Int32Dtype(),
                 'Link1': pd.Int16Dtype(),
                 'Link2': pd.Int16Dtype(),
                 'match score': float
                 }

    for file in xi_files:

        print('Reading xi-file: {}'.format(file))
        try:
            s = pd.read_csv(Help.compatible_path(file), delimiter=',', dtype=xi_dtypes)
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

    xtable['rawfile'] = xtable['Source'].apply(_rawfile_from_source)

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        Help.assign_type, axis=1)

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        pd.Series(np.vectorize(Help.generate_id,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                 index=xtable.index).replace('nan', np.nan)

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for inter xlink to inter/intra/homomultimeric
        onlyInter = xtable['type'] == 'inter'
        xtable.loc[onlyInter, 'type'] =\
            np.vectorize(Help.categorize_inter_peptides)(xtable[onlyInter]['prot1'],
                                                         xtable[onlyInter]['pos1'],
                                                         xtable[onlyInter]['pepseq1'],
                                                         xtable[onlyInter]['prot2'],
                                                         xtable[onlyInter]['pos2'],
                                                         xtable[onlyInter]['pepseq1'])
        print('[Xi Read] categorized inter peptides')
    else:
        print('[Xi Read] skipped inter peptide categorization')

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearch'

    xtable = Help.order_columns(xtable, col_order, compact)

    return xtable

class XiConfig:
    def __init__(self, config_file):
        self.filepath = config_file
        self.config = {}
        self.load_config()

    def load_config(self):
        with open(self.filepath, 'r') as f:
            lines = f.readlines()

        p_col = re.compile(":")  # matches ':' not preceded or followed by ':'
        p_semicol = re.compile(";")  # matches ';'

        for line in lines:
            if line.startswith('#'):
                continue
            if line.strip() == '':  # skip empty lines
                continue

            colon_pos = [m.start() for m in p_col.finditer(line)]
            semicolon_pos = [m.start() for m in p_semicol.finditer(line)]

            if len(semicolon_pos) == 0 and len(colon_pos) == 1:  # e.g. topmgchits:150
                # split on colon pos
                key = line[:colon_pos[0]].strip()
                value = line[colon_pos[0] + 1:].strip()
                self.config[key] = value

            elif len(semicolon_pos) == 0 and len(colon_pos) > 1:  # e.g. crosslinker:NonCovalentBound:Name:NonCovalent
                split_colon = colon_pos[1] if len(colon_pos) % 2 == 1 else colon_pos[0]
                key = line[:split_colon].strip()
                value = dict([list(map(str.strip, line[split_colon + 1:].split(':', 1)))])
                if key in self.config:
                    self.config[key].append(value)
                else:
                    self.config[key] = [value]

            elif len(semicolon_pos) > 0 and len(
                    colon_pos) > 1:  # e.g. modification:variable::SYMBOLEXT:bs3_hyd;MODIFIED:S,T,Y;DELTAMASS:156.0786347
                first_semicolon_pos = semicolon_pos[0]
                last_colon_before_semicolon_pos = [pos for pos in colon_pos if pos < first_semicolon_pos][-2]
                key = line[:last_colon_before_semicolon_pos].strip(' :')  # strip leading and trailing spaces and colons
                value = {x.split(':')[0].strip(): x.split(':')[1].strip() if ':' in x else True for x in
                         line[last_colon_before_semicolon_pos + 1:].split(';')}
                if key in self.config:
                    self.config[key].append(value)
                else:
                    self.config[key] = [value]

            else:
                raise ValueError(f"Invalid line format: {line.strip()}")


if __name__ == '__main__':

    # defines the column headers required for xtable output
    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    xi_file = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\Xi\XI_results_XiVersion1.6.739.csv'

    xtable = Read(xi_file)
