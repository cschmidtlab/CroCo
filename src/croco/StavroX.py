# -*- coding: utf-8 -*-

"""
Functions to read StavroX processed crosslink data.
"""

import numpy as np
import pandas as pd

import os, re
if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf


def _type_from_proteins(protein1_string, protein2_string):
    """
    Decide whether its a mono, loop, intra or inter-protein link
    based on the fields Protein 1 and Protein 2 of a StavroX results
    table

    Args:
        protein1_string (str): A Protein1 string from StavroX
        protein2_string (str): A Protein2 string from StavroX
    Returns:
        str: type of cross-link (mono, loop, intra, inter)
    """

    if protein2_string == 'H2O':
        return 'mono'
    elif protein2_string == 'intrapeptidal':
        return 'loop'
    elif protein1_string == protein2_string:
        return 'intra'
    else:
        return 'inter'

def _clear_protname(protname_string):
    """
    Clear fasta-header and return the cleaned string. Also remove StavroX
    specific protein2-terms like H20 and intrapeptidal

    Args:
        protname_string (str): A StavroX-Protein string

    Returns:
        str: cleaned protein string (my be empty string)
    """
    if protname_string.startswith('>'):
        protname_string = protname_string.split(' ')[0]
        return protname_string[1:]
    elif protname_string == 'H2O':
        return ''
    elif protname_string == 'intrapeptidal':
        return ''
    else:
        return protname_string

def _clear_xlink(xlink_string):
    """
    Removes all non-numerical characters from a string.

    Args:
        xlink_string (str): Contains the xlink position with the linked AA
    Returns:
        int: the xlink position without the amino acid label
    """
    xlink_string = re.sub('\D+', '', xlink_string)

    return int(xlink_string)

def _calc_xpos2(t, pos1, pos2, xlink2):
    """
    Calculate the absolute position of a beta-peptide cross-link
    based on what type the xlink is

    Args:
        t(str): type of the crosslink (mono, loop, inter/intra/homomultimeric)
        pos1(int): absolute position of the first AA of the alpha-peptide
        pos2(int): absolute position of the first AA of the beta-peptide
        xlink2(int): relative position of the cross-link in the beta-peptide
    Returns:
        int: absolute position of the cross-link in the beta-peptide
    """

    try:
        if t == 'mono':
            return np.nan
        elif t == 'loop':
            return int(pos1) + int(xlink2) - 1
        else:
            return int(pos2) + int(xlink2) - 1
    except:
        return np.nan

def _mods_and_sequences_from_peptides(peptide1, peptide2, mod_dict):
    """
    Extract modification position and name from two StavroX peptide strings.
    Set peptide2 to np.nan for monolinks and to peptide1 for loop-links

    Args:
        peptide1 (str): Entry of the StavroX Peptide 1 column
        peptide2 (str): Entry of the StavroX Peptide 2 column
        mod_dict (dict): dict mapping modification abbreviations to lists of [Modified AA, Modification name, Modification mass]
    Returns:
        list: name of the modification(s) of peptide1
        list: position of the modification(s) of peptide1
        list: mass(es) of the modification(s) of peptide1
        list: name of the modification(s) of peptide2
        list: position of the modification(s) of peptide2
        list: mass(es) of the modification(s) of peptide2
    """

    if str(peptide2) == '0': # mono-link
        peptide2 = np.nan
        (mod1, modpos1, modmass1, sequence1), (mod2, modpos2, modmass2, sequence2) =\
            _calculate_mod_modpos_modmass(peptide1), ([], [], [])
    elif str(peptide2) == '1': # loop-link
        (mod1, modpos1, modmass1, sequence1), (mod2, modpos2, modmass2, sequence2) =\
            _calculate_mod_modpos_modmass(peptide1, mod_dict), _calculate_mod_modpos_modmass(peptide1, mod_dict)
    else: # regular inter-peptide link
        (mod1, modpos1, modmass1, sequence1), (mod2, modpos2, modmass2, sequence2) =\
            _calculate_mod_modpos_modmass(peptide1, mod_dict), _calculate_mod_modpos_modmass(peptide2, mod_dict)

    return mod1, modpos1, modmass1, sequence1, mod2, modpos2, modmass2, sequence2

def _calculate_mod_modpos_modmass(peptide_string, mod_dict):
    """
    Extraction modification, name of the modification and modification mass
    from a StavroX peptide string

    Args:
        peptide_string (str): Entry of the StavroX Peptide 1/2 column
        mod_dict (dict): dict mapping modification abbreviations to lists of [Modified AA, Modification name, Modification mass]

    Returns:
        list: name of the modification(s)
        list: position of the modification(s)
        list: mass(es) of the modification(s)
        str: sequence without modification(s)
    """

    # if there is no peptide, there are no modifications
    try:
        a_float = float(peptide_string)
        if np.isnan(a_float):
            return [], [], []
        else:
            raise Exception('Unusual peptide string recognised: {}'.format(a_float))
    except:
        modpos = []
        modmass= []
        mod = []

        sequence = ''
        for idx, char in enumerate(peptide_string):
            if char in mod_dict.keys():
                # avoid setting N- and C-terminal ends as modification
                mod.append(mod_dict[char][1])
                modmass.append(mod_dict[char][2])
                modpos.append(idx)
                char = mod_dict[char][0]
            if char in '[]{}':
                continue
            sequence += char

        return mod, modpos, modmass, sequence


def _parse_ssf(ssf_file):
    """
    Parses a StavroX properties.ssf file to extract modification identifiers
    and masses

    Args:
        ssf_file (str): path to the properties.ssf file from StavroX

    Returns:
        dict: Dict mapping modification identifiers a list of [modified AA, name of modification, modmass]
    """
    ssf_dict = {}
    current_list = None

    with open(ssf_file, 'r') as f:
        line = f.readline()
        while line:
            if 'ELEMENTS' in line:
                current_list = 'elements'
                ssf_dict[current_list] = []
            elif 'AMINOACIDS' in line:
                current_list = 'AAs'
                ssf_dict[current_list] = []
            elif 'VARMODIFICATION' in line:
                current_list = 'varmods'
                ssf_dict[current_list] = []
            elif 'STATMODIFICATION' in line:
                current_list = 'fixedmods'
                ssf_dict[current_list] = []
            elif 'END' in line:
                current_list = None
            elif current_list != None:
                ssf_dict[current_list].append(line.strip().split(';'))
            line = f.readline()

    elements2mass = dict((x, y) for x, y in ssf_dict['elements'])

    AA2formula = dict((y,z) for x,y,z in ssf_dict['AAs'])
    AA2name = dict((y,x) for x,y,z in ssf_dict['AAs'])
    
    mod_dict = dict()
    for x,y,z in ssf_dict['varmods']:
        mod_dict[y] = [x,]
    for x,y in ssf_dict['fixedmods']:
        mod_dict[y] = [x,]

    for symbol in mod_dict.keys():
        modformula = AA2formula[symbol]
        elements, stoichiometries = _elemental_composition(modformula, elements2mass)
        modmass = 0
        for element, amount in zip(elements, stoichiometries):
            modmass += int(amount) * float(elements2mass[element])
            
        unmodformula = AA2formula[mod_dict[symbol][0]]
        elements, stoichiometries = _elemental_composition(unmodformula, elements2mass)
        unmodmass = 0
        for element, amount in zip(elements, stoichiometries):
            unmodmass += int(amount) * float(elements2mass[element])
        
        mod_dict[symbol].append(AA2name[symbol])
        mod_dict[symbol].append(modmass-unmodmass)

#
#    mod2mass = {}
#    mod2name = {}
#
#    for AA in AA2formula:
#        if AA not in 'RHKDESTNQCUGPAVILMFYW':
#            formula = AA2formula[AA]
#            elements, stoichiometries = _elemental_composition(formula, elements2mass)
#            mass = 0
#            for element, amount in zip(elements, stoichiometries):
#                mass += int(amount) * float(elements2mass[element])
#
#            mod2mass[AA] = mass
#            mod2name[AA] = AA2name[AA]

    return mod_dict


def _elemental_composition(formula, elements2mass):
    """
    Compute two lists with elements and amount of atoms from
    a chemical cormula

    Args:
        formula (string): a chemical formula like C6H12N2O
        elements2mass (dict): a dict mapping elements to mass

    Returns:
        elements: a list of the elements in the formula
        stoichiometries: a list of the amounts of these elements
    """

    elements = []
    this_stoichiometrie = None
    stoichiometries = []
    collected_stoichiometry = True # set initial to true to avoid
    # setting element 1 to 1

    # formula is e.g. C6H12N2O
    for char in formula:
        # check if character is an element
        if char in elements2mass.keys():
            elements.append(char)
            # only enter if this_stoichiometrie has been filled
            if this_stoichiometrie:
                # if stoichiometrie is not given as 1
                stoichiometries.append(this_stoichiometrie)
                this_stoichiometrie = ''
            else:
                this_stoichiometrie = ''
                if collected_stoichiometry is False:
                    stoichiometries.append(1)
            collected_stoichiometry = False
        else:
            this_stoichiometrie += char
            collected_stoichiometry = True
    if collected_stoichiometry is False:
        stoichiometries.append(1)
    else:
        stoichiometries.append(this_stoichiometrie)

    return elements, stoichiometries

def Read(stavrox_files, ssf_file, col_order=None, compact=False):
    """
    Collect data from StavroX spectrum search and return an xtable data array.

    Args:
        stavrox_files: path or list of paths to StavroX output file(s)
        ssf_file: properties.ssf to load modification IDs and masses
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to compact the xTable to only those columns listed in col_order
    Returns:
        pandas.DataFrame: xtable data table
    """

    # convert to list if the input is only a single path
    if not isinstance(stavrox_files, list):
        stavrox_files = [stavrox_files]

    allData = list()

    stavrox_dtypes = {'Scan number': str,
                      'Charge': 'int16',
                      'Protein 1 From': pd.Int64Dtype(),
                      'Protein 2 From': pd.Int64Dtype(),
                      'Protein 1': str,
                      'Protein 2': str,
                      'Peptide 1': str,
                      'Peptide 2': str,
                      'best linkage position peptide 1': str,
                      'best linkage position peptide 2': str,
                      'Score': float}

    for file in stavrox_files:

        print('Reading StavroX-file: {}'.format(file))

        try:
            with open(hf.compatible_path(file), 'r') as f:
                firstline = f.readline()

            headers = list()

            for element in firstline.split(';'):
                if element not in ['From', 'To']:
                    # there is a typo in the StavroX output files
                    if element == 'Peptide2':
                        headers.append('Peptide 2')
                        last_saved = 'Peptide 2'
                    else:
                        headers.append(element)
                        last_saved = element
                else:
                    headers.append(last_saved + ' ' + element)
            # the spectrum UUID column contains a semicolon delimiter! This messes up the whole csv reading
            headers.insert(-1, 'Spectrum UUID2')

            print(headers)

            # Reassign the column headers to avoid duplicate From and To fields
            s = pd.read_csv(hf.compatible_path(file),
                            delimiter=';',
                            header=0,
                            index_col = False,
                            dtype=stavrox_dtypes,
                            names = headers)

            allData.append(s)
        except:
            raise Exception('[StavroX Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData)

    ### Process the data to comply to xTable format
    xtable = xtable.rename(columns={'Protein 1 From': 'pos1',
                                    'Protein 2 From': 'pos2',
                                    'Score': 'score'
                                    })

    # the field Scan number contains the mgf file header. Use Regex to extract
    # scan no
    xtable[['rawfile', 'scanno', 'prec_ch']] = xtable['Scan number'].str.extract(hf.regexDict['mgfTITLE'])

    print('[StavroX Read] Parsed MGF title')

    # calculate the type of line (i.e. mono, loop, intra or inter)
    xtable['type'] = np.vectorize(_type_from_proteins)(xtable['Protein 1'], xtable['Protein 2'])

    print('[StavroX Read] inferred type')

    # Set pos1 to 1 if Nterminal cross-link
    xtable['pos1'] = xtable['pos1'].replace(0, 1)
    xtable['pos2'] = xtable['pos2'].replace(0, 1)

    print('[StavroX Read] changed xlink positions')

    # remove for example preceding > in UniProt headers
    xtable['prot1'] = xtable['Protein 1'].apply(_clear_protname)
    xtable['prot2'] = xtable['Protein 2'].apply(_clear_protname)

    print('[StavroX Read] Cleared Protein names')

    # Best linkage position also contains the linked AA (that is already given
    # by sequence and link position)
    # --> xtract only the numerical part
    # StavroX treats the N-terminus as 0th position --> replace by 1
    xtable['xlink1'] = xtable['best linkage position peptide 1'].apply(_clear_xlink).replace(0, 1)
    xtable['xlink2'] = xtable['best linkage position peptide 2'].apply(_clear_xlink).replace(0, 1)

    print('[StavroX Read] Found xlink position')

    # calculate absolute position of xlink as sum of start of peptide
    # and relative position of the xlink
    xtable['xpos1'] = xtable['xlink1'].astype(int) + xtable['pos1'].astype(int) - 1

    # xpos2 has to be calculated separately for inter/intra, loop and mono-peptides
    xtable['xpos2'] =\
        np.vectorize(_calc_xpos2)(xtable['type'], xtable['pos1'], xtable['pos2'], xtable['xlink2'])

    print('[StavroX Read] Generated xpos')

    mod_dict = _parse_ssf(hf.compatible_path(ssf_file))

    print('[StavroX Read] parsed SSF')

    # Extract the modification mass and position from the peptide string
    xtable[['mod1', 'modpos1', 'modmass1', 'pepseq1', 'mod2', 'modpos2', 'modmass2', 'pepseq2']] =\
        pd.DataFrame(xtable[['Peptide 1', 'Peptide 2']].apply(\
            lambda row: _mods_and_sequences_from_peptides(row['Peptide 1'],
                                           row['Peptide 2'],
                                           mod_dict),
                    axis=1).tolist(), index=xtable.index)

    print('[StavroX Read] Extracted modifications and sequences')

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        pd.Series(np.vectorize(hf.generate_id,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                 index=xtable.index).replace('nan', np.nan)

    print('[StavroX Read] Generated ID')

    # Stavrox does not run on isotope-labeled xlinkers
    xtable['xtype'] = np.nan

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for inter xlink to inter/intra/homomultimeric
        onlyInter = xtable['type'] == 'inter'
        xtable.loc[onlyInter, 'type'] =\
            np.vectorize(hf.categorize_inter_peptides)(xtable[onlyInter]['prot1'],
                                                     xtable[onlyInter]['pos1'],
                                                     xtable[onlyInter]['pepseq1'],
                                                     xtable[onlyInter]['prot2'],
                                                     xtable[onlyInter]['pos2'],
                                                     xtable[onlyInter]['pepseq1'])
        print('[StavroX Read] categorized inter peptides')
    else:
        print('[StavroX Read] skipped inter peptide categorization')

    # StavroX already filters decoys
    xtable['decoy'] = False

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    xtable['search_engine'] = 'StavroX'

    xtable = hf.order_columns(xtable, col_order, compact)

    return xtable

if __name__ == '__main__':

    # defines the column headers required for xtable output
    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2',
                  'type', 'score',
                  'ID',
                  'pos1', 'pos2', 'decoy']

    stavrox_files = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\stavrox\20180615_KS_CL_9_pXtract.csv'

    ssf_file = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\stavrox\StavroxSettings.ssf'
    xtable = Read(stavrox_files, ssf_file, col_order=col_order, compact=False)
