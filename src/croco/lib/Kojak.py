# -*- coding: utf-8 -*-

"""
Functions to read Kojak data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import re

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def calc_pos_from_xpos(xpos, xlink):
    """
    Calculates the absolute position of the first AA of a peptide
    sequence from the absolute position of the cross-link AA (xpos) its
    position within the sequence (xlink)

    Returns: pos - Absolute position of AA in sequence
    """
    xpos = int(xpos)
    xlink = int(xlink)

    return xpos - xlink + 1

def process_kojak_peptide(peptide_string):
    """
    Return Modifications and the peptide sequence
    from a Kojak sequence string such as M[15.99]TDSKYFTTNK
    """

    pattern = re.compile('\[(.*?)\]')
    mods = re.findall(pattern, peptide_string)

    pattern = re.compile('([A-Z]+)')
    sequence = ''.join(re.findall(pattern, peptide_string))

    if mods == []:
        mods = np.nan

    return mods, sequence

def process_kojak_protein(protein_string):
    """
    Return protein name and absolute cross-link position from
    a kojak string such as
    sp|P07340|AT1B1_RAT Sodium/potassium-transporting ATPase subunit beta-1 OS=Rattus norvegicus GN=Atp1(13);
    """
    pattern = re.compile('([^ ]+) .+?(?:\((\d+)\))?;')
    if pattern.match(protein_string):
        match = pattern.match(protein_string)
        prot, xpos = match.groups()
        if xpos == None: # re.match returns None (not NaN) if a substring doesnt match
            return prot, np.nan
        else:
            return prot, xpos
    else:
        return np.nan, np.nan


def Read(kojak_file, rawfile=None, keep=False):
    """
    reads pLink results file and returns an xtable data array.

    :params: koajk_file: path to Kojak results file
    :params: rawfile: name of the corresponding rawfile

    :returns: xtable data table
    """

    ### Collect data and convert to pandas format

    print('Reading Kojak-file: ' + kojak_file)

    # only called if inter_file is not None
    if kojak_file:
        data = pd.read_csv(kojak_file,
                           skiprows = 1, # skip the Kojak version
                           delimiter='\t')
    else:
        return FileNotFoundError('Kojak txt file not found')

    ### Convert data inside pandas df

    # remove lines containing non-identified PSMs (marked with '-' in both
    # Link columns
    data = data[(data['Link #1'] != '-') | (data['Link #2'] != '-')]

    # transfer scan no from read data to xtable
    xtable = pd.DataFrame(data['Scan Number'])

    # rename column
    xtable.columns = ['scanno']

    xtable['prec_ch'] = data['Charge']

    # apply the function and assign the result to multiple new columns
    xtable['mod1'], xtable['pepseq1'] =\
           zip(*data['Peptide #1'].apply(process_kojak_peptide))

    xtable['mod2'], xtable['pepseq2'] =\
           zip(*data['Peptide #2'].apply(process_kojak_peptide))

    xtable['xlink1'] = data['Link #1']
    xtable['xlink2'] = data['Link #2']

    # transform unset xlinks to np.nan
    xtable.replace('-1', np.nan, inplace=True)

    xtable['prot1'], xtable['xpos1'] =\
            zip(*data['Protein #1'].apply(process_kojak_protein))

    xtable['prot2'], xtable['xpos2'] =\
            zip(*data['Protein #2'].apply(process_kojak_protein))

    xtable['score'] = data['Score']

    # set a decoy indicator where at least one protein is reversed
    xtable['decoy'] = np.where(xtable['prot1'].str.contains('REVERSE') |\
                                 xtable['prot2'].str.contains('REVERSE'),
                                 True, False)

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable.loc[xtable['prot2'].notnull(), 'type'] = 'inter'
    xtable.loc[(xtable['prot2'].notnull())\
        & (xtable['prot1'] == xtable['prot2']), 'type'] = 'intra'
    xtable.loc[xtable['prot2'].isnull() & xtable['xlink2'].notnull(), 'type'] = 'loop'
    xtable.loc[xtable['xlink1'] == '-1', 'type'] = 'mono'

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    # Must be calculated as float as NaN is not implemented in int
    xtable['pos1'] = xtable['xpos1'].astype(float, errors='ignore') - \
                     xtable['xlink1'].astype(float, errors='ignore') + 1
    xtable['pos2'] = xtable['xpos2'].astype(float, errors='ignore') - \
                     xtable['xlink2'].astype(float, errors='ignore') + 1

    # set the rawfile name for xtable (None if not provided by call)
    xtable['rawfile'] = rawfile

    xtable['xtype'] = np.nan
    
    xtable['search_engine'] = 'Kojak'
    
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

    ### return xtable df

    return xtable