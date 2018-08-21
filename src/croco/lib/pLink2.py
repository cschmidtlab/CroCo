# -*- coding: utf-8 -*-
"""
Read pLink2 files
@author: aretaon
"""

import pandas as pd
import os, sys
import re
import numpy as np

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

def plink2peptide2pandas(filepath):
    """
    Read a pLink peptide results file and return a pandas dictionary

    :param filepath: Path to a pLink results file e.g. filtered_cross-linked_peptides.csv

    :returns: pandas dataframe
    """

    with open(filepath, 'r') as fh:

        # read the first header-line into list
        headers1 = fh.readline().strip().split(',')
        # read the second header-line and remove it from list
        # avoid the first entry as it is only the line-indicator
        headers2 = fh.readline().strip().split(',')[1:]

        data = {} # init of data dict for pandas
        for h in headers1 + headers2:
            data[h] = []

        # read the first line
        line = fh.readline()
        while line:

            # avoid reading in an empty line
            if line == '':
                continue

            if line[0].isdigit(): # indicates a headers1 line
                # save the line and use it when printing all following lines
                # corresponding to that title-line
                line1_data = line.strip().split(',')

            else:
                line2_data = line.strip().split(',')[1:]
                for i in range(len(line1_data)):
                    # iterate through the data and link them to their resp
                    # header by index
                    data[headers1[i]].append(line1_data[i])
                for i in range(len(line2_data)):
                    data[headers2[i]].append(line2_data[i])

            # read the next line
            line=fh.readline()

        # remove empty keys from dict
        if '' in data:
            del data['']

    return pd.DataFrame.from_dict(data)


def process_plink2_title(spec_string):
    """
    Extract rawfile name, precursor charge and scan no from pLink sequence
    string

    :returns: list of rawfile, scanno, prec_ch
    """
    # the pattern of the title string is 20171215_JB04_Sec06.10959.10959.2
    # in pLink 2.3 and 20171215_JB04_Sec06.10959.10959.2.0 in pLink 2.1
    pextract_pattern = re.compile('(.+?)\.\d+\.(\d+)\.(\d+)\.*\d*')
    if pextract_pattern.match(spec_string):
        match = pextract_pattern.match(spec_string)
        # rawfile, scanno, prec_ch
        return match.groups()
    else:
        return np.nan

def process_plink2_sequence(row):
    """
    Extract peptide sequences and cross-link positions from
    pLink sequence string e.g. YVPTAGKLTVVILEAK(7)-LTVVILEAK(2):1
    Can differentiate between mono, loop and cross-link information

    :param row: a row of a dataframe with headers 'Peptide' and 'type'

    :returns: list of pepseq1, pepseq2, xpos1, xpos2, xtype
    """

    #TODO: How to recognize heavy labels?

    xtype = 0 # all cross-links are light

    if row['type'] == 'inter':
        pattern = re.compile(r'(\w+)\((\d+)\)-(\w*)\((\d*)\)')
        match = pattern.match(row['Peptide'])
        pepseq1, xlink1, pepseq2, xlink2 = match.groups()

    elif row['type'] == 'loop':
        pattern = re.compile(r'(\w+)\((\d+)\)\((\d*)\)')
        match = pattern.match(row['Peptide'])
        pepseq1, xlink1, xlink2 = match.groups()
        pepseq2 = np.nan

    elif row['type'] == 'mono':
        pattern = re.compile(r'(\w+)\((\d*)\)')
        match = pattern.match(row['Peptide'])
        pepseq1, xlink1 = match.groups()
        xlink2, pepseq2 = [np.nan] * 2

    else:
        pepseq1, xlink1, pepseq2, xlink2 = [np.nan] * 4

    return pepseq1, xlink1, pepseq2, xlink2, xtype

def process_plink2_proteins(row):
    """
    Extract protein name and absolute cross-link position from
    pLink protein string e.g.
    sp|P63045|VAMP2_RAT(79)-sp|P63045|VAMP2_RAT(59)
    """

    if row['type'] == 'inter':
        pattern = re.compile('(.+?)\((\d+)\)-?([^\(]*)\(?(\d*)\)?')
        match = pattern.match(row['Proteins'])
        prot1, xpos1, prot2, xpos2 = match.groups()

        xpos1 = int(xpos1)
        xpos2 = int(xpos2)
        prot1 = str(prot1.strip())
        prot2 = str(prot2.strip())

    elif row['type'] == 'loop':
        pattern = re.compile(r'(.+?)\((\d+)\)\((\d*)\)')
        match = pattern.match(row['Proteins'])
        prot1, xpos1, xpos2 = match.groups()
        xpos1 = int(xpos1)
        xpos2 = int(xpos2)
        prot1 = str(prot1.strip())
        prot2 = prot1

    elif row['type'] == 'mono':
        pattern = re.compile(r'(.+?)\((\d+)\)')
        match = pattern.match(row['Proteins'])
        prot1, xpos1 = match.groups()
        xpos1 = int(xpos1)
        xpos2, prot2 = [np.nan] * 2
        prot1 = str(prot1.strip())

    else:
        prot1, xpos1, pepseq2, xpos2 = [np.nan] * 4

    return prot1, xpos1, prot2, xpos2

def calculate_abs_pos(row):
    """
    Return the absolute position of the first AA of both peptides.
    If only one peptide present (mono-link) return NaN for the second position
    """

    pos1 = int(row['xpos1']) - int(row['xlink1']) + 1

    if np.isnan(float(row['xpos2'])):
        pos2 = np.nan
    else:
        pos2 = int(row['xpos2']) - int(row['xlink2']) + 1

    return pos1, pos2

def read_plink_modifications(filepath):
    """
    Open a pLink modification.ini file and extract all modifications with
    their names as dict.

    :params: filepath: Path to modifications.ini

    :returns: mod_dict: Dict mapping names to masses
    """

    pattern = re.compile(r'^(.*)=\w+ \w+ (-?[0-9]\d*\.\d+)? -?[0-9]\d*\.\d+')
    mod_dict = {}

    with open(filepath, 'r') as f:
        for line in f:
            if pattern.match(line):
                match = pattern.match(line)
                name, mass = match.groups()
                mod_dict[name] = mass

    return mod_dict

def Read(plinkdir, compact=False):
    """
    reads pLink2 report dir and returns an xtabel data array.

    Args:
        plinkdir: plink2 reports subdir (reports)
        keep (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: data table
    """

    ### Collect data, convert to pandas format and merge

    # Initialise file names as None to use implicit booleaness
    inter_file = None
    loop_file = None
    mono_file = None

    for e in os.listdir(plinkdir):
        if 'filtered_cross-linked_peptides.csv' in e:
            inter_file = e

        if 'filtered_loop-linked_peptides' in e:
            loop_file = e

        if 'filtered_mono-linked_peptides' in e:
            mono_file = e

    frames = []
    # only called if inter_file is not None
    if inter_file:
        if os.path.exists(os.path.join(plinkdir, inter_file)):
            print('Reading pLink inter-file: ' + inter_file)
            inter_df = plink2peptide2pandas(os.path.join(plinkdir, inter_file))
            inter_df['type'] = 'inter'
            frames.append(inter_df)
    if loop_file:
        if os.path.exists(os.path.join(plinkdir, loop_file)):
            print('Reading pLink loop-file: ' + loop_file)
            loop_df = plink2peptide2pandas(os.path.join(plinkdir, loop_file))
            loop_df['type'] = 'loop'
            frames.append(loop_df)
    if mono_file:
        if os.path.exists(os.path.join(plinkdir, mono_file)):
            print('Reading pLink mono-file: ' + mono_file)
            mono_df =  plink2peptide2pandas(os.path.join(plinkdir, mono_file))
            mono_df['type'] = 'mono'
            frames.append(mono_df)

    xtable = pd.concat(frames)

    ### Convert data inside pandas df

    # split title column into three
    xtable['rawfile'], xtable['scanno'], xtable['prec_ch'] =\
        zip(*xtable['Title'].apply(process_plink2_title))

    # Directly assign the re group matches into new columns
    xtable['pepseq1'], xtable['xlink1'], xtable['pepseq2'],\
    xtable['xlink2'], xtable['xtype'] =\
        zip(*xtable.apply(process_plink2_sequence,
                          axis=1))

    xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'] =\
            zip(*xtable.apply(process_plink2_proteins,
                              axis=1))

    xtable['score'] = xtable['Score']

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        np.vectorize(hf.generateID)(xtable['type'], xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'])

    # calculate absolute position of first AA of peptide
    xtable['pos1'], xtable['pos2'] = zip(*xtable.apply(calculate_abs_pos,
                                                       axis=1))
    # add a label referring to the ordering in the pLink results table
    xtable['Order'] = xtable[['Peptide_Order', 'Spectrum_Order']].apply(lambda x: ','.join(x), axis=1)

    # set the sequence of loop links to be the same as the corresponding pepseq1
    xtable.loc[xtable['type'] == 'loop', 'pepseq2'] = xtable[xtable['type'] == 'loop']['pepseq1']

    # manually set decoy to reverse as pLink hat its own internal target-decoy
    # algorithm
    xtable['decoy'] = False

    # Reassign the type for inter xlink to inter/intra/homomultimeric
    xtable.loc[xtable['type'] == 'inter', 'type'] =\
        np.vectorize(hf.categorizeInterPeptides)(xtable[xtable['type'] == 'inter']['prot1'],
                                                 xtable[xtable['type'] == 'inter']['pos1'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'],
                                                 xtable[xtable['type'] == 'inter']['prot2'],
                                                 xtable[xtable['type'] == 'inter']['pos2'],
                                                 xtable[xtable['type'] == 'inter']['pepseq1'])

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    ## generate the mod_dict linking pLink modification names to masses
    # in case of calling croco from the source folder structure...
    file_dir, file_name = os.path.split(__file__)
    if os.path.exists(os.path.join(file_dir,
                                   '../data/modification.ini')):
        modifi_dir = os.path.abspath(os.path.join(file_dir,
                                                  '../data/modification.ini'))
    # ... or calling from a exe-file in a folder-setup with the data folder at top-level
    elif os.path.exists('./data/modification.ini'):
        modifi_dir = os.path.abspath('./data/modification.ini')
    # ... or calling from within a single bundled exe-file
    else:
        try:
            # PyInstaller creates a temp folder and stores its path in _MEIPASS
            base_path = sys._MEIPASS
            modifi_dir =  os.path.abspath(\
                os.path.join(base_path, './data/modification.ini'))
        # ... or something went wrong
        except:
            raise Exception('Modifications.ini not found')

    # load pLink modifications.ini from data-folder
    mod_dict = read_plink_modifications(os.path.abspath(modifi_dir))

    # extract modification information
    pattern = re.compile(r'(.*)\((\d+)\)')

    pepseq1 = xtable['pepseq1'].tolist()
    Modifications = xtable['Modifications'].tolist()

    if len(pepseq1) == len(Modifications):
        print('Len of pepseq1 and Modification match!')
    else:
        print('Len of pepseq1 and Modification dont match!')

    modmass1 = []
    mod1 = []
    modpos1 = []
    modmass2 = []
    mod2 = []
    modpos2 = []

    # iterate over all lines in the input file
    for idx, modstr in enumerate(Modifications):

        this_modmass1 = []
        this_mod1 = []
        this_modmass2 = []
        this_mod2 = []
        this_modpos1 = []
        this_modpos2 = []

        # Extract annotations from every item in the modstring
        for mod in modstr.split(';'):

            if pattern.match(mod):
                match = pattern.match(mod)
                mod, modpos = match.groups()

                # transform modification names to masses
                try:
                    mass = mod_dict[mod]
                except:
                    # use the input string if no subsitution found
                    mass = mod

                seqlen1 = len(pepseq1[idx])
                # pLink assigns additional modification position to the C-term
                # of the first peptide, the xlinker and the N-term of the
                # second peptide
                if int(modpos) > (seqlen1 + 3):
                    this_mod2.append(mod)
                    this_modpos2.append(int(modpos) - (seqlen1 + 3))
                    this_modmass2.append(mass)
                # C-term of first peptide
                elif int(modpos) == (seqlen1 + 1):
                    this_mod2.append(mod)
                    this_modpos2.append(int(modpos) - 1)
                    this_modmass2.append(mass)
                # cannot assign modifications to xlinker in xTable
                elif int(modpos) == (seqlen1 + 2):
                    pass
                # Modification on N-term of second peptide
                elif int(modpos) == (seqlen1 + 3):
                    this_mod2.append(mod)
                    this_modpos2.append(1)
                    this_modmass2.append(mass)
                else:
                    this_mod1.append(mod)
                    this_modpos1.append(modpos)
                    this_modmass1.append(mass)

        # multiple modifications of one peptide are stored as ;-delimited strings
        modmass1.append(this_modmass1)
        mod1.append(this_mod1)
        modpos1.append(this_modpos1)
        modmass2.append(this_modmass2)
        mod2.append(this_mod2)
        modpos2.append(this_modpos2)

    xtable['modmass1'] = modmass1
    xtable['mod1'] = mod1
    xtable['modpos1'] = modpos1
    xtable['modmass2'] = modmass2
    xtable['mod2'] = mod2
    xtable['modpos2'] = modpos2

    xtable['search_engine'] = 'pLink2'

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')
    
    xtable = hf.applyColOrder(xtable, col_order, compact)

    ### return xtable df
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
    
    init(col_order)
    
    xtable = Read(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\pLink2\Sec23_24_only\reports')