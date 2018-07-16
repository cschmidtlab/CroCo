# -*- coding: utf-8 -*-

"""
Functions to read pLink1 data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import os
import sys
import re

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def plinkprotein2pandas(filepath):
    """
    Read a pLink protein results file and return a pandas dictionary

    :params: filepath: Path to a pLink results file e.g. _inter_combine.protein.xls

    :returns: pandas dataframe
    """
    with open(filepath, 'r') as fh:

        data = {} # init of data dict for pandas

        # initialise the protein_header line
        header1 = ''
        # initialise the peptide header-line
        header2 = ''

        # protein level header starts with Order tag
        header1_exp = re.compile(r'^Order\s.*')
        # protein level entry starts with the order entry
        entry1_exp = re.compile(r'^\d+\s.*')
        # peptide level header starts with tab
        header2_exp = re.compile(r'^\s+Order.*')
        # peptide entry starts with whitespace foloowed by order entry
        entry2_exp = re.compile(r'^\s+\d+.*')

        for line in fh.readlines():

            # beginning of a protein-block
            if header1_exp.match(line):
                # only read the header-line on its first occurence
                if header1 == '':
                    # Get header names
                    header1 = line.strip().split('\t')
                    # set the column names for the dict
                    for h in header1:
                        data[h] = []

            # body of a protein-block
            elif entry1_exp.match(line):
                # Read and store header1 data to write with every case of header2
                entry1_data = line.strip().split('\t')

            # beginning of a peptide block
            elif header2_exp.match(line):
                if header2 == '':
                    header2 = line.strip().split('\t')
                    header2 = [x if x != 'Order' else 'Order2' for x in header2 ]
                    for h in header2:
                        data[h] = []

            # beginning of a peptide entry
            elif entry2_exp.match(line):
                entry2_data = line.strip().split('\t')
                # add the corresponding level1 entries
                for idx, d in enumerate(entry1_data):
                    data[header1[idx]].append(d)
                for idx, d in enumerate(entry2_data):
                    data[header2[idx]].append(d)

    return pd.DataFrame.from_dict(data)

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

def process_plink_sequence(seq_string):
    """
    Extract peptide sequences and cross-link positions from
    pLink sequence string e.g. YVPTAGKLTVVILEAK(7)-LTVVILEAK(2):1

    :returns: list of pepseq1, pepseq2, xpos1, xpos2, xtype
    """
    pattern = re.compile('(\w+)\((\d+)\)-(\w+)\((\d+)\):(\d+)')
    try:
        match = pattern.match(seq_string)
        # pepseq1, xpos1, pepseq2, xpos2, xtype
        return match.groups()

    except Exception as e:
        print(e)
        return np.nan

def process_plink_spectrum(spec_string):
    """
    Extract rawfile name, precursor charge and scan no from pLink sequence
    string

    :returns: list of rawfile, scanno, prec_ch
    """
    pextract_pattern = re.compile('(.+)\.(\d+)\.\d+\.(\d+)\.dta')
    if pextract_pattern.match(spec_string):
        match = pextract_pattern.match(spec_string)
        # rawfile, scanno, prec_ch
        return match.groups()
    else:
        return np.nan

def process_plink_proteins(prot_string):
    """
    Extract protein name and absolute cross-link position from
    pLink protein string e.g.
    sp|P63045|VAMP2_RAT(79)-sp|P63045|VAMP2_RAT(59)
    """
    pattern = re.compile('(.+?)\((\d+)\)-?([^\(]*)\(?(\d*)\)?')
    match = pattern.match(prot_string)
    return match.groups()

def Read(plinkdir):
    """
    reads pLink report dir and returns an xtabel data array.

    :params: plinkdir: plink report subdir (e.g. sample1)

    :returns: xtable data table
    :returns: xinfo meta-data object
    """

    ### Collect data, convert to pandas format and merge

    # Initialise file names as None to use implicit booleaness
    inter_file = None
    loop_file = None
    mono_file = None

    for e in os.listdir(plinkdir):
        if '_inter_combine.protein.xls' in e:
            inter_file = e

        if '_loop_combine.protein.xls' in e:
            loop_file = e

        if '_mono_combine.protein.xls' in e:
            mono_file = e

    frames = []
    # only called if inter_file is not None
    if inter_file:
        print('Reading pLink inter-file: ' + inter_file)
        inter_df = plinkprotein2pandas(os.path.join(plinkdir, inter_file))
        inter_df['type'] = 'inter'
        frames.append(inter_df)
    if loop_file:
        print('Reading pLink loop-file: ' + loop_file)
        loop_df = plinkprotein2pandas(os.path.join(plinkdir, loop_file))
        loop_df['type'] = 'loop'
        frames.append(loop_df)
    if mono_file:
        print('Reading pLink mono-file: ' + mono_file)
        mono_df =  plinkprotein2pandas(os.path.join(plinkdir, mono_file))
        mono_df['type'] = 'mono'
        frames.append(mono_df)

    data = pd.concat(frames)

    ### Convert data inside pandas df


    # init xtable with column containing lists of rawfile, scanno, prec_ch
    xtable = pd.DataFrame(data['Spectrum'].apply(process_plink_spectrum))
    # split column into three
    xtable['rawfile'], xtable['scanno'], xtable['prec_ch'] =\
        zip(*xtable['Spectrum'])
    # drop the original column
    xtable.drop('Spectrum',
                axis = 1,
                inplace=True)

    # Directly assign the re group matches into new columns
    xtable['pepseq1'], xtable['xlink1'], xtable['pepseq2'],\
    xtable['xlink2'], xtable['xtype'] =\
        zip(*data['Sequence'].apply(process_plink_sequence))

    xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'] =\
            zip(*data['Proteins'].apply(process_plink_proteins))

    xtable['type'] = data['type']
    xtable['score'] = data['Score']

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    xtable['pos1'] = xtable['xpos1'].astype(int, errors='ignore') - \
                     xtable['xlink1'].astype(int, errors='ignore') + 1
    xtable['pos2'] = xtable['xpos2'].astype(int, errors='ignore') - \
                     xtable['xlink2'].astype(int, errors='ignore') + 1

    # add a lobel referring to the ordering in the pLink results table
    xtable['Order'] = data[['Order', 'Order2']].apply(lambda x: ','.join(x), axis=1)

    # add a path to the plink PSM image
    xtable['PSM image'] = os.path.join(plinkdir, 'psm') + os.path.sep +\
                            data['Spectrum'].astype(str) + '.png'
    # clear path slashes and assign an excel suitable hyperlink
    xtable['PSM image'] = xtable['PSM image'].astype(str).\
                            apply(lambda x: '=HYPERLINK("' + os.path.abspath(x) +\
                                  '", "' + x.split(os.path.sep)[-1] + '")')

    # manually set decoy to reverse as pLink hat its own internal target-decoy
    # algorithm
    xtable['decoy'] = False

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    # save original score in columns
    xtable['plink score'] = xtable['score']

    # compute a minus log P score for better comparison with higher=better scores
    xtable['score'] = -np.log(xtable['score'])

    # generate the mod_dict linking pLink modification names to masses
    
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
    pattern = re.compile(r'(\d+),.*\((.*)\)')

    pepseq1 = xtable['pepseq1'].tolist()
    Modification = data['Modification'].tolist()

    if len(pepseq1) == len(Modification):
        print('Len of pepseq1 and Modification match!')
    else:
        print('Len of pepseq1 and Modification dont match!')

    mod1 = []
    modmass1 = []
    modpos1 = []
    modmass2 = []
    modpos2 = []
    mod2 = []

    # iterate over all lines in the input file
    for idx, modstr in enumerate(Modification):

        this_mod1 = []
        this_modmass1 = []
        this_modpos1 = []
        this_modmass2 = []
        this_modpos2 = []
        this_mod2 = []

        # Extract annotations from every item in the modstring
        for mod in modstr.split(';'):

            if pattern.match(mod):
                match = pattern.match(mod)
                modpos, mod = match.groups()

                # transform modification names to masses
                try:
                    mass = mod_dict[mod]
                except:
                    # use the input string if no subsitution found
                    mass = mod
                
                seqlen1 = len(pepseq1[idx]) + 1
                if int(modpos) > seqlen1:
                    this_modpos2.append(str(int(modpos) - seqlen1))
                    this_modmass2.append(mass)
                    this_mod2.append(mod)
                else:
                    this_modpos1.append(modpos)
                    this_modmass1.append(mass)
                    this_mod1.append(mod)

        # multiple modifications of one peptide are stored as ;-delimited strings
        modmass1.append(';'.join(this_modmass1))
        modpos1.append(';'.join(this_modpos1))
        mod1.append(';'.join(this_mod1))
        modmass2.append(';'.join(this_modmass2))
        modpos2.append(';'.join(this_modpos2))
        mod2.append(';'.join(this_mod2))

    xtable['mod1'] = mod1
    xtable['modmass1'] = modmass1
    xtable['modpos1'] = modpos1
    xtable['mod2'] = mod2
    xtable['modmass2'] = modmass2
    xtable['modpos2'] = modpos2

    # reorder columns
    # append pLink specific columns
    col_order.extend(['PSM image', 'plink score'])
    col_order[0] = 'Order' # append to front

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')
    
    # reorder columns to start with the xtable columns
    all_cols = list(xtable.columns.values)
    remaining_cols = [x for x in all_cols if x not in col_order]
    new_order = col_order + remaining_cols

    xtable = xtable[new_order]
    ### return xtable df

    return xtable