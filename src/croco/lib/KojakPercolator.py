# -*- coding: utf-8 -*-

"""
Functions to read Percolator prcoesse Kojak data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import re

from . import HelperFunctions as HeFn

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
    try:
        if np.isnan(float(xpos)):
            return np.nan

        xpos = int(xpos)
        xlink = int(xlink)
            
        return xpos - xlink + 1
    
    except Exception as e:
        print('{}: xpos was {} and xlink was {}'.format(e, xpos, xlink))
        return np.nan

def process_kojak_peptide(peptide_string):
    """
    Return Modifications and the peptide sequence
    from a Kojak sequence string such as M[15.99]TDSKYFTTNK
    """

    pattern = re.compile('\[(.*?)\]')
    modmasses = re.findall(pattern, peptide_string)

    # define positions for modifications
    count = 0
    modpositions = []
    for e in peptide_string:
        if e in 'RHKDESTNQCUGPAVILMFYW':
            count += 1
        elif e == '[': # start of modification
            modpositions.append(count)

    pattern = re.compile('([A-Z]+)')
    sequence = ''.join(re.findall(pattern, peptide_string))

    mods = np.nan
    if modpositions != []:
        mods = 'Unnamed'

    return mods, modmasses, modpositions, sequence

def process_kojak_protein(protein_string):
    """
    Return protein name and absolute cross-link position from
    a kojak string such as
    sp|P07340|AT1B1_RAT Sodium/potassium-transporting ATPase subunit beta-1 OS=Rattus norvegicus GN=Atp1(13);
    SPA_STAAU(260);
    """
    
    pattern = re.compile('([^ ]*)(?: ?.*)(?:\((\d+)\));?$')
 
    if protein_string != '':
        if pattern.match(protein_string):
            match = pattern.match(protein_string)
            prot, xpos = match.groups()
            if xpos == None: # re.match returns None (not NaN) if a substring doesnt match
                return prot, np.nan
            else:
                return prot, xpos
        else:
            return np.nan, np.nan

def generate_ID(type, prot1, xpos1, prot2, xpos2):
    """
    Return a link ID based on the type of the xlink
    """

    if type in ['mono', 'loop']:
        return '-'.join([str(prot1), str(xpos1)])
    else:
        return '-'.join([str(prot1), str(xpos1), str(prot2), str(xpos2)])

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

def split_concaten(dataframe, where, delimiter=';'):
    """
    Splits each row of a dataframe that contains a delimiter-separated
    string into two columns with each element of the string in each row.
    
    Args:
        dataframe: dataframe to operate on
        where: column-name in which to find the strings
        delimiter: (optional) the delimiter-character to look for
    
    Returns:
        dataframe: Modified dataframe
    """
    

def Read(perc_file, percolator_string='.validated', rawfile=None, keep=False):
    """
    Collects unprocessed and percolated results and returns an xtable data array.

    Args:
        perc_file: path to percolated Kojak file
        percolator_string: user-defined string appended to the percolated filenames
        rawfile: name of the corresponding rawfile

    Returns:
        xtable: xtable data table
    """
    
    ### Collect data and convert to pandas format

    print('Reading Percolated Percolator-file: ' + perc_file)

    # only called if inter_file is not None

    percolated = pd.read_csv(perc_file,
                             delimiter='\t',
                             usecols=range(5),
                             index_col=False, # avoid taking the first col as index
                             engine='python')
                             
    percolated.rename(columns={'PSMId': 'SpecId'}, inplace=True)
        
    unperc_file = perc_file.replace(percolator_string, '')

    print('Reading Unpercolated Percolator-file: ' + unperc_file)

    try:
        unpercolated = pd.read_csv(unperc_file,
                                  delimiter = '\t',
                                  usecols=range(10),
                                  engine='python',
                                  index_col=False)
    except:
        raise FileNotFoundError(unperc_file)
    
    # Merge with left join (only keys that are in tje percolated DF will be re-
    # tained)
    data = pd.merge(percolated, unpercolated, on='SpecId', how='left')
        
    # Reading the Kojak-file is required to get additional information on the
    # matches such as the corresponding protein names
    kojak_file = unperc_file[0:unperc_file.find('.perc')] + '.kojak.txt'
    
    print('Reading Kojak-file: ' + kojak_file)

    try:
        kojak = pd.read_csv(kojak_file,
                            skiprows = 1, # skip the Kojak version
                            delimiter='\t')
    except:
        raise FileNotFoundError(kojak_file)
        
    kojak.rename(columns={'Scan Number': 'scannr'}, inplace=True)

    data = pd.merge(data, kojak, on=['scannr', 'Charge', 'dScore', 'Score'], how='left')

    # split ambiguous concatenated protein names
    data = HeFn.split_concatenated_lists(data, where=['Protein #1', 'Protein #2'])

    print(sum(data['split_entry'] == False))

    ### Process the data to comply to xTable format
    xtable = data.rename(columns={'scannr': 'scanno',
                                  'Charge': 'prec_ch',
                                  'Link #1': 'xlink1',
                                  'Link #2': 'xlink2',
                                  })
    
    # change the identifier for empty entries
    xtable['xlink1'].replace(to_replace='-1', value=np.nan, inplace=True)
    xtable['xlink2'].replace(to_replace='-1', value=np.nan, inplace=True)

    # apply the function and assign the result to multiple new columns
    xtable['mod1'], xtable['modmass1'], xtable['modpos1'], xtable['pepseq1'] =\
           zip(*xtable['Peptide #1'].apply(process_kojak_peptide))

    xtable['mod2'], xtable['modmass2'], xtable['modpos2'], xtable['pepseq2'] =\
           zip(*xtable['Peptide #2'].apply(process_kojak_peptide))

    xtable['prot1'], xtable['xpos1'] =\
            zip(*xtable['Protein #1'].apply(process_kojak_protein))

    xtable['prot2'], xtable['xpos2'] =\
            zip(*xtable['Protein #2'].apply(process_kojak_protein))

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        assign_type, axis=1)
    
    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        np.vectorize(generate_ID)(xtable['type'], xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'])


    xtable['decoy'] = 'T-' in xtable['SpecId']

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    # Must be calculated as float as NaN is not implemented in int
    xtable['pos1'] = xtable.apply(lambda row: calc_pos_from_xpos(row['xpos1'], row['xlink1']), axis=1)
    xtable['pos2'] = xtable.apply(lambda row: calc_pos_from_xpos(row['xpos2'], row['xlink2']), axis=1)

    # set the rawfile name for xtable (None if not provided by call)
    xtable['rawfile'] = rawfile

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'Kojak and Percolator'

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
    import os
    
    # defines the column headers required for xtable output
    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']
    
    os.chdir(r'Z:\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Kojak')
    perc_file = r'Z:\02_experiments\05_croco_dataset\002_20180425\crosslink_search\Kojak\20180518_JB_jb05a_l50.perc.inter.validated.txt'
    
    perc = Read(perc_file)
    
    perc.to_excel('test.xls',
                  index=False)