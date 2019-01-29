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

def rawfile_and_scanno_from_Scan(scan_number_string):
    """
    Extract rawfile name and scan number from a string like
    20180518_JB_jb05a_l50.10069.10069.5.dta
    
    Args:
        scan_number_string (str): A spectrum ID string from StavroX
    
    Returns:
        str or np.nan: rawfile
        int or np.nan: scan number
    """

    pattern = re.compile(r'^([^\.]+)\.(\d+)\.\d+\.\d+\.\w+$')
    if pattern.match(scan_number_string):
        m = pattern.match(scan_number_string)
        rawfile, scanno = m.group(1,2)
        scanno = int(scanno)
        return rawfile, scanno
    else:
        return np.nan, np.nan

def type_from_proteins(protein1_string, protein2_string):
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

def clear_protname(protname_string):
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

def clear_xlink(xlink_string):
    """
    Removes all non-numerical characters from a string.
    
    Args:
        xlink_string (str): Contains the xlink position with the linked AA
    Returns:
        int: the xlink position without the amino acid label
    """
    xlink_string = re.sub('\D+', '', xlink_string)

    return int(xlink_string)

def calc_xpos2(t, pos1, pos2, xlink2):
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
            return 0
        elif t == 'loop':
            return int(pos1) + int(xlink2)
        else:
            return int(pos2) + int(xlink2)
    except:
        return 0

def pepseqs_from_peptides(peptide1, peptide2):
    """
    Extract the peptide sequences for peptide 1 and peptide 2 from the
    StavroX Peptide entry.
    Set peptide2 to np.nan for monolinks and to peptide1 for loop-links

    Args:
        peptide1 (str): Stavrox peptide string for peptide 1
        peptide2 (str): Stavrox peptide string for peptide 2

    Returns:
        pepseq1 (str): corrected sequence for peptide1
        pepseq2 (str): corrected sequence for peptide2
    """

    def clear_pepname(peptide_string):
        """
        Remove the preceding and trailing identifiers for free N- and C-termini
        []{} denote N- and C-termini in StavroX
        """

        if peptide_string[0] in '[]{}':
            peptide_string = peptide_string[1:]
        if peptide_string[-1] in '[]{}':
            peptide_string = peptide_string[:-1]
        # make all upper case to allow parsing in later programmes
        return peptide_string.upper()

    if str(peptide2) == '0': # mono-link
        peptide2 = np.nan
        return clear_pepname(peptide1), ''
    elif str(peptide2) == '1': # loop-link
        peptide2 = peptide1
        return clear_pepname(peptide1), clear_pepname(peptide2)
    else: # regular inter-peptide link
        return clear_pepname(peptide1), clear_pepname(peptide2)

def mods_from_peptides(peptide1, peptide2, mod2name, mod2mass):
    """
    Extract modification position and name from two StavroX peptide strings.
    Set peptide2 to np.nan for monolinks and to peptide1 for loop-links

    Args:
        peptide1 (str): Entry of the StavroX Peptide 1 column
        peptide2 (str): Entry of the StavroX Peptide 2 column
        mod2name (dict): dict mapping modification abbreviations to names
        mod2mass (dict): dict mapping modification abbreviations to masses
    Returns:
        list: name of the modification(s) of peptide1
        list: position of the modification(s) of peptide1
        list: mass(es) of the modification(s) of peptide1
        list: name of the modification(s) of peptide2
        list: position of the modification(s) of peptide2
        list: mass(es) of the modification(s) of peptide2
    """

    def calculate(peptide_string, mod2name=mod2name, mod2mass=mod2mass):
        # if there is no peptide, there are no modifications
        try:
            a_float = float(peptide_string)
            if np.isnan(a_float):
                return [], [], []
            else:
                raise Exception('Unusual peptide string recognised: {}'.format(a_float))
        except:
            # clean up peptide string
            if peptide_string[1] in '[]{}':
                peptide_string = peptide_string[1:]
            if peptide_string[-1] in '[]{}':
                peptide_string = peptide_string[:-1]

            modpos = []
            modmass= []
            mod = []

            position = 0
            for char in peptide_string:
                if char in mod2name.keys():
                    # avoid setting N- and C-terminal ends as modification
                    if char not in '[]{}':
                        position += 1
                        mod.append(mod2name[char])
                        modmass.append(mod2mass[char])
                        modpos.append(position)

            return mod, modpos, modmass

    if str(peptide2) == '0': # mono-link
        peptide2 = np.nan
        (mod1, modpos1, modmass1), (mod2, modpos2, modmass2) =\
            calculate(peptide1), ([], [], [])
    elif str(peptide2) == '1': # loop-link
        peptide2 = peptide1
        (mod1, modpos1, modmass1), (mod2, modpos2, modmass2) =\
            calculate(peptide1), calculate(peptide2)
    else: # regular inter-peptide link
        (mod1, modpos1, modmass1), (mod2, modpos2, modmass2) =\
            calculate(peptide1), calculate(peptide2)


    return mod1, modpos1, modmass1, mod2, modpos2, modmass2

def ParseSSF(ssf_file):
    """
    Parses a StavroX properties.ssf file to extract modification identifiers
    and masses

    Args:
        ssf_file (str): path to the properties.ssf file from StavroX

    Returns:
        dict: Dict mapping modification identifiers to their mass
        dict: Dict mapping modification identifiers to their name
    """

    def ElementalComposition(formula, elements2mass):
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
            elif 'END' in line:
                current_list = None
            elif current_list != None:
                ssf_dict[current_list].append(line.strip().split(';'))
            line = f.readline()

    elements2mass = dict((x, y) for x, y in ssf_dict['elements'])

    AA2formula = dict((y,z) for x,y,z in ssf_dict['AAs'])
    AA2name = dict((y,x) for x,y,z in ssf_dict['AAs'])

    mod2mass = {}
    mod2name = {}

    for AA in AA2formula:
        if AA not in 'RHKDESTNQCUGPAVILMFYW':
            formula = AA2formula[AA]
            elements, stoichiometries = ElementalComposition(formula, elements2mass)
            mass = 0
            for element, amount in zip(elements, stoichiometries):
                mass += int(amount) * float(elements2mass[element])

            mod2mass[AA] = mass
            mod2name[AA] = AA2name[AA]

    return mod2mass, mod2name

def Read(stavrox_files, ssf_file, col_order=None, compact=False):
    """
    Collect data from StavroX spectrum search and return an xtable data array.

    Args:
        stavrox_files: path to StavroX output file
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
    
    for file in stavrox_files:

        print('Reading StavroX-file: {}'.format(stavrox_files))

        try:
            # Reassign the column headers to avoid duplicate From and To fields
            s = pd.read_csv(hf.FSCompatiblePath(file),
                               delimiter=';',
                               header=0,
                               index_col = False,
                               names = ['Score',
                                        'm/z',
                                        'Charge',
                                        'M+H+',
                                        'Calculated MassDeviation in ppm',
                                        'Deviation in ppm',
                                        'Peptide 1',
                                        'Protein 1',
                                        'Protein 1 From',
                                        'Protein 1 To',
                                        'Peptide 2',
                                        'Protein 2',
                                        'Protein 2 From',
                                        'Protein 2 To',
                                        'Scan number',
                                        'is Selected in Table',
                                        'Candidate identifier',
                                        'Folder Number',
                                        'Retention time in sec',
                                        'miscellaneous',
                                        'best linkage position peptide 1',
                                        'best linkage position peptide 2',
                                        'All linkage positions',
                                        'Spectrum UUID',
                                        'local False discovery rate',
                                        'Unknown'])
    
            allData.append(s)
        except:
            raise Exception('[StavroX Read] Failed opening file: {}'.format(file))
    
    xtable = pd.concat(allData)
    
    ### Process the data to comply to xTable format
    xtable = xtable.rename(columns={'Scan': 'scanno',
                                    'Charge': 'prec_ch',
                                    'Protein 1 From': 'pos1',
                                    'Protein 2 From': 'pos2',
                                    'Score': 'score'
                                    })

    # calculate the type of line (i.e. mono, loop, intra or inter)
    xtable['type'] = np.vectorize(type_from_proteins)(xtable['Protein 1'], xtable['Protein 2'])

    # Set pos1 to 1 if Nterminal cross-link
    xtable.loc[(xtable['type'] != 'mono') & (xtable['pos1'] == 0), ['pos1', 'pos2']] = 1


    # Extract the Peptide sequence e.g. from [KPDT]AGT]
    #xtable[['pepseq1', 'pepseq2']] = data[['Peptide 1', 'Peptide 2']].apply(pepseqs_from_peptides, axis=1)
    xtable['pepseq1'], xtable['pepseq2'] =\
        np.vectorize(pepseqs_from_peptides)(xtable['Peptide 1'], xtable['Peptide 2'])

    # extract rawfile name and scan no from the Scan header e.g.
    # 20180518_JB_jb05a_m100.10636.10636.2.dta
    xtable['rawfile'], xtable['scanno'] =\
        zip(*xtable['Scan number'].apply(rawfile_and_scanno_from_Scan))

    # remove for example preceding > in UniProt headers
    xtable['prot1'] = xtable['Protein 1'].apply(clear_protname)
    xtable['prot2'] = xtable['Protein 2'].apply(clear_protname)

    # Best linkage position also contains the linked AA (that is already given
    # by sequence and link position)
    # --> xtract only the numerical part
    xtable['xlink1'] = xtable['best linkage position peptide 1'].apply(clear_xlink)
    xtable['xlink2'] = xtable['best linkage position peptide 2'].apply(clear_xlink)

    # calculate absolute position of xlink as sum of start of peptide
    # and relative position of the xlink
    xtable['xpos1'] = xtable['xlink1'].astype(int) + xtable['pos1'].astype(int)

    # xpos2 has to be calculated separately for inter/intra, loop and mono-peptides
    xtable['xpos2'] =\
        np.vectorize(calc_xpos2)(xtable['type'], xtable['pos1'], xtable['pos2'], xtable['xlink2'])

    mod2mass, mod2name = ParseSSF(hf.FSCompatiblePath(ssf_file))

    # Extract the modification mass and position from the peptide string
    xtable['mod1'], xtable['modpos1'], xtable['modmass1'],\
    xtable['mod2'], xtable['modpos2'], xtable['modmass2'] =\
        zip(*xtable[['Peptide 1', 'Peptide 2']].apply(\
            lambda row: mods_from_peptides(row['Peptide 1'],
                                           row['Peptide 2'],
                                           mod2name,
                                           mod2mass),
           axis=1))


    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        np.vectorize(hf.generateID)(xtable['type'], xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'])

    # Stavrox does not run on isotope-labeled xlinkers
    xtable['xtype'] = np.nan

    # Reassign the type for inter xlink to inter/intra/homomultimeric
    xtable.loc[xtable['type'] == 'inter', 'type'] =\
        np.vectorize(hf.categorizeInterPeptides)(xtable[xtable['type'] == 'inter']['prot1'],
                                                 xtable[xtable['type'] == 'inter']['pos1'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'],
                                                 xtable[xtable['type'] == 'inter']['prot2'],
                                                 xtable[xtable['type'] == 'inter']['pos2'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'])

    # StavroX already filters decoys
    xtable['decoy'] = False

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
                  'xpos2',
                  'type', 'score',
                  'ID',
                  'pos1', 'pos2', 'decoy']

    os.chdir(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\StavroX')
    stavrox_files = [r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\StavroX\l50.csv',
                     r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\StavroX\l100.csv']

    ssf_file = r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\StavroX\properties.ssf'
    stvrx = Read(stavrox_files, ssf_file, col_order=col_order, compact=False)

    stvrx.to_excel('test.xls',
                   index=False)