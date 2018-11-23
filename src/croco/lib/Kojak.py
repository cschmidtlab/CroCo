# -*- coding: utf-8 -*-

"""
Functions to read Kojak data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import re

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

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

    mods = []
    sequence = ''
    is_mod = False

    for char in peptide_string:
        if char == '[':
            is_mod = True
            theMod = ''
        elif char == ']':
            is_mod = False
            mods.append(theMod)
        elif is_mod == False:
            if char.isalpha():
                sequence += char
        else:
            theMod += char

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


def Read(kojak_file, rawfile=None, compact=False):
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
        return FileNotFoundError('Kojak txt file not found: {}'.format(kojak_file))

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
        np.vectorize(hf.generateID)(xtable['type'], xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'])

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    # Must be calculated as float as NaN is not implemented in int
    xtable['pos1'] = xtable['xpos1'].astype(float, errors='ignore') - \
                     xtable['xlink1'].astype(float, errors='ignore') + 1
    xtable['pos2'] = xtable['xpos2'].astype(float, errors='ignore') - \
                     xtable['xlink2'].astype(float, errors='ignore') + 1

    # Reassign the type for inter xlink to inter/intra/homomultimeric
    xtable.loc[xtable['type'] == 'inter', 'type'] =\
        np.vectorize(hf.categorizeInterPeptides)(xtable[xtable['type'] == 'inter']['prot1'],
                                                 xtable[xtable['type'] == 'inter']['pos1'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'],
                                                 xtable[xtable['type'] == 'inter']['prot2'],
                                                 xtable[xtable['type'] == 'inter']['pos2'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'])
    # set the rawfile name for xtable (None if not provided by call)
    xtable['rawfile'] = rawfile

    xtable['xtype'] = np.nan
    
    xtable['search_engine'] = 'Kojak'
    
    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')
    
    xtable = hf.applyColOrder(xtable, col_order, compact)

    ### return xtable df
    return xtable

if __name__ == '__main__':
    infile = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\kojak\2017_08_04_SVs_BS3_13.kojak.txt'
    xtable = Read(infile, rawfile='Test')
    