# -*- coding: utf-8 -*-
"""
Functions that are collectively used by croco.Kojak and croco.KojakPercolator.
"""

import re
import numpy as np
import pandas as pd


if __name__ == '__main__' or __name__ == 'KojakFunctions':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def extract_peptide(xtable):
    """
    Extract peptide sequence, modification mass and position from the
    Peptide #1 and Peptide #2 entries

    Args:
        xtable (pandas.DataFrame): xTable data structure with "Peptide #1" and "Peptide #2" columns
    Returns:
        pandas.DataFrame: xTable with modmass, modpos, pepseq and mod
    """
    pep1notNull = xtable['Peptide #1'].notnull()
    pep2notNull = xtable['Peptide #2'].notnull()

    # the index corresponds to the index of the slice of the dataframe
    # as original row numbers are retained during conversion, values can directly
    # be inserted at the right row
    xtable[['modmass1', 'modpos1', 'pepseq1']] =\
        pd.DataFrame(xtable.loc[pep1notNull, 'Peptide #1'].apply(process_kojak_peptide).tolist(),
                     index=xtable.loc[pep1notNull, 'Peptide #1'].index)

    if sum(pep2notNull) > 0:
        xtable[['modmass2', 'modpos2', 'pepseq2']] =\
            pd.DataFrame(xtable.loc[pep2notNull, 'Peptide #2'].apply(process_kojak_peptide).tolist(),
                         index=xtable.loc[pep2notNull, 'Peptide #2'].index)
    else:
        xtable['modmass2'] = np.nan
        xtable['modpos2'] = np.nan
        xtable['pepseq2'] = np.nan

    # use the modification masses as labels
    xtable['mod1'] = xtable['modmass1'].apply(lambda x: x if hf.isnan(x) else [str(y) for y in x])
    xtable['mod2'] = xtable['modmass2'].apply(lambda x: x if hf.isnan(x) else [str(y) for y in x])

    return xtable

def extract_protein(xtable):
    """
    Extract protein name and relative cross-link position from the Protein #
    entries

    Args:
        xtable (pandas.DataFrame): xTable data structure with "Protein #1", "Protein #2", xpos1, xlink1, and xlink2 columns
    Returns:
        pandas.DataFrame: xTable with prot and xpos
    """
    xtable[['prot1', 'xpos1']] =\
        xtable['Protein #1'].str.extract(r'^(\w+)(?:\((.*?)\))?;$')

    xtable.loc[xtable['xpos1'].notnull(), 'xpos1'] =\
        xtable.loc[xtable['xpos1'].notnull(), 'xpos1'].astype(int)
    xtable['xpos1'] = xtable['xpos1'].astype(pd.Int64Dtype())

    xtable[['prot2', 'xpos2']] =\
        xtable['Protein #2'].str.extract(r'^(\w+)(?:\((.*?)\))?;$')

    xtable.loc[xtable['xpos2'].notnull(), 'xpos2'] =\
        xtable.loc[xtable['xpos2'].notnull(), 'xpos2'].astype(int)
    xtable['xpos2'] = xtable['xpos2'].astype(pd.Int64Dtype())

    # xpos2 for loop links is not directly stored but can be inferred
    # from xpos1 and xlink2
    isLoopLink = xtable['prot1'].notnull() & xtable['prot2'].isnull()
    xtable.loc[isLoopLink, 'xpos2'] = xtable.loc[isLoopLink, ['xpos1', 'xlink2']].sum(axis=1).astype(pd.Int64Dtype()) - xtable.loc[isLoopLink, 'xlink1']
    
    return xtable

def assign_ID_and_type(xtable):
    """
    Calculate if a cross link is of inter or of loop type
    Refine the inter type into inter/intra/homomultimeric
    Generate ID for the xlinks

    Args:
        xtable (pandas.DataFrame): Table data structure with "prot", "pos", "pepseq"
    Returns:
        pandas.DataFrame: xTable with type and ID
    """

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable.loc[xtable['prot2'].notnull(), 'type'] = 'inter'
    xtable.loc[xtable['prot2'].isnull() & xtable['xlink2'].notnull(), 'type'] = 'loop'
    # Kojak does not generate monolinked peptides but peptides modified
    # with the hydrolysed xlinker mass
    xtable.loc[xtable['xlink1'].isnull() & xtable['xlink2'].isnull(), 'type'] = 'linear or mono'

    # Reassign the type for inter xlink to inter/intra/homomultimeric
    isInterLink = xtable['type'] == 'inter'
    # only perform if the selection is not all false
    if sum(isInterLink) > 0:
        xtable.loc[isInterLink, 'type'] =\
            np.vectorize(hf.categorize_inter_peptides)(xtable[isInterLink]['prot1'],
                                                       xtable[isInterLink]['pos1'],
                                                       xtable[isInterLink]['pepseq1'],
                                                       xtable[isInterLink]['prot2'],
                                                       xtable[isInterLink]['pos2'],
                                                       xtable[isInterLink]['pepseq2'])

    # only apply the operation requiring at least prot1 and xpos1 to those
    # lines that are loop, intra or interlinks
    type_identified = xtable['type'].notna()
    # generate an ID for every crosslink position within the protein(s)
    xtable.loc[type_identified, 'ID'] =\
        pd.Series(np.vectorize(hf.generate_id,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                 index=xtable.index).replace('nan', np.nan)

    return xtable

def set_decoy(xtable, decoy_string):
    """
    sets the column decoy based on whether the decoy string is present in the
    protein name or not

    Args:
        xtable (pandas.DataFrame): xTable with "prot" columns titles
        decoy_string (str): Kojak decoy string
    Returns:
        pandas.DataFrame: xTable with decoy column
    """
    # Check if all prot2 are null (may be in only loop dfs)
    if xtable['prot2'].isnull().all():
        xtable['decoy'] = np.where(xtable['prot1'].str.contains(decoy_string), True, False)
    else:
        # set a decoy indicator where at least one protein is reversed
        xtable['decoy'] = np.where(xtable['prot1'].str.contains(decoy_string) |\
                                     xtable['prot2'].str.contains(decoy_string),
                                     True, False)

    return xtable

def process_kojak_peptide(peptide_string):
    """
    Return Modifications, their localisation and the peptide sequence
    from a Kojak sequence string such as M[15.99]TDSKYFTTNK.

    If modifications are found, two lists with modification masses, positions
    and the raw peptide sequence are returned.
    If no modififications are found within a peptide string, the function
    returns np.nan, np.nan and the sequence.

    Args:
        peptide_string (str): a Kojak peptide string
    Returns:
        list of float or np.nan: list of modification masses
        list of int or np.nan: list of modification positions within the peptide
        str: peptide sequence without modifications
    """

    modmasses = []
    sequence = ''
    modposns = []
    is_mod = False

    posInStr = 0
    for char in peptide_string:
        if char == '[':
            is_mod = True
            theMod = ''
        elif char == ']':
            is_mod = False
            modmasses.append(float(theMod))
            modposns.append(int(posInStr))
        elif is_mod == False:
            if char.isalpha():
                sequence += char
                posInStr += 1
        else:
            theMod += char

    if modmasses == []:
        modmasses = np.nan

    return modmasses, modposns, sequence

#def process_kojak_protein(protein_string):
#    """
#    Return protein name and absolute cross-link position from
#    a kojak string such as
#    sp|P07340|AT1B1_RAT Sodium/potassium-transporting ATPase subunit beta-1 OS=Rattus norvegicus GN=Atp1(13);
#
#    Args:
#        protein_string(str): a kojak protein string
#
#    Returns:
#        str or np.nan: protein name
#        int or np.nan: position
#    """
#    # RE: group1: everything until the first (lazy) brackets
#    # group2 (optional) everything inside the brackets
#    pattern = re.compile('^([^\(]+?)(?:\((\d*)\))?;')
#    if pattern.match(protein_string):
#        match = pattern.match(protein_string)
#        prot, xpos = match.groups()
#        if xpos == None: # re.match returns None (not NaN) if a substring doesnt match
#            return prot, np.nan
#        else:
#            return prot, int(xpos)
#    else:
#        return np.nan, np.nan
