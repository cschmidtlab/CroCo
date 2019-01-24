# -*- coding: utf-8 -*-
"""
Kojak Helper Functions

@author: User
"""

import re
import numpy as np
import pandas as pd


if __name__ == '__main__' or __name__ == 'KojakFunctions':
    import HelperFunctions as hf
    import KojakFunctions as kj
else:
    from . import HelperFunctions as hf
    from . import KojakFunctions as kj

def extract_peptide(xtable):
    """
    Extract peptide sequence, modification mass and position from the
    Peptide #1 and Peptide #2 entries
    """
    xtable['modmass1'], xtable['modpos1'], xtable['pepseq1'] =\
        [pd.Series(x) for x in zip(*xtable['Peptide #1'].apply(kj.process_kojak_peptide))]

    xtable['modmass2'], xtable['modpos2'], xtable['pepseq2'] =\
           [pd.Series(x) for x in zip(*xtable['Peptide #2'].apply(kj.process_kojak_peptide))]

    # use the modification masses as labels
    xtable['mod1'] = xtable['modmass1'].apply(lambda x: x if hf.isNaN(x) else [str(y) for y in x])
    xtable['mod2'] = xtable['modmass2'].apply(lambda x: x if hf.isNaN(x) else [str(y) for y in x])

    return xtable

def extract_protein(xtable):
    """
    extract protein name and relative cross-link position from the Protein #
    entries
    """
    xtable['prot1'], xtable['xpos1'] =\
           [pd.Series(x) for x in zip(*xtable['Protein #1'].apply(kj.process_kojak_protein))]

    # calculate xpos2 as if both peptides are identical
    # if the peptides are not, xpos2 is replaces below
    xtable['xpos2'] = xtable['xpos1'] + xtable['xlink2'] - xtable['xlink1']
    xtable['prot2'] = np.nan

    # prot2 and xpos2 for inter and intra cross-links are derived from Protein #2
    # select the lines that have a Protein #2 entry
    hasSecProtein = xtable['Protein #2'] != '-'
    if sum(hasSecProtein) > 0:
        prot2, xpos2 =\
            [pd.Series(x) for x in zip(*xtable.loc[hasSecProtein, 'Protein #2'].apply(kj.process_kojak_protein))]
        xtable.loc[hasSecProtein, 'prot2'] = prot2
        xtable.loc[hasSecProtein, 'xpos2'] = xpos2

    return xtable

def assign_ID_and_type(xtable):
    """
    Calculate if a cross link is of inter or of loop type
    Refine the inter type into inter/intra/homomultimeric
    Generate ID for the xlinks
    """
    
    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable.loc[xtable['prot2'].notnull(), 'type'] = 'inter'
    xtable.loc[xtable['prot2'].isnull() & xtable['xlink2'].notnull(), 'type'] = 'loop'
    # TODO: Kojak does not generate monolinked peptides but peptides modified
    # with the hydrolysed xlinker mass

    # Reassign the type for inter xlink to inter/intra/homomultimeric
    isInterLink = xtable['type'] == 'inter'
    # only perform if the selection is not all false
    if sum(isInterLink) > 0:
        xtable.loc[isInterLink, 'type'] =\
            np.vectorize(hf.categorizeInterPeptides)(xtable[isInterLink]['prot1'],
                                                     xtable[isInterLink]['pos1'],
                                                     xtable[isInterLink]['pepseq1'],
                                                     xtable[isInterLink]['prot2'],
                                                     xtable[isInterLink]['pos2'],
                                                     xtable[isInterLink]['pepseq1'])

    # only apply the operation requiring at least prot1 and xpos1 to those
    # lines that are loop, intra or interlinks
    type_identified = xtable['type'].notna()
    # generate an ID for every crosslink position within the protein(s)
    xtable.loc[type_identified, 'ID'] =\
        np.vectorize(hf.generateID)(xtable.loc[type_identified, 'type'],
                                    xtable.loc[type_identified, 'prot1'],
                                    xtable.loc[type_identified, 'xpos1'],
                                    xtable.loc[type_identified, 'prot2'],
                                    xtable.loc[type_identified, 'xpos2'])

    return xtable

def set_decoy(xtable, decoy_string):
    """
    sets the column decoy based on whether the decoy string is present in the
    protein name or not
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
    from a Kojak sequence string such as M[15.99]TDSKYFTTNK
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

def process_kojak_protein(protein_string):
    """
    Return protein name and absolute cross-link position from
    a kojak string such as
    sp|P07340|AT1B1_RAT Sodium/potassium-transporting ATPase subunit beta-1 OS=Rattus norvegicus GN=Atp1(13);
    """
    # RE: group1: everything until the first (lazy) brackets
    # group2 (optional) everything inside the brackets
    pattern = re.compile('^([^\(]+?)(?:\((\d*)\))?;')
    if pattern.match(protein_string):
        match = pattern.match(protein_string)
        prot, xpos = match.groups()
        if xpos == None: # re.match returns None (not NaN) if a substring doesnt match
            return prot, np.nan
        else:
            return prot, int(xpos)
    else:
        return np.nan, np.nan
