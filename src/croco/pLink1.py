# -*- coding: utf-8 -*-

"""
Functions to read pLink1 data.

"""

import numpy as np
import pandas as pd

import os
import sys
import re

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Read(plinkdirs, col_order=None, compact=False):
    """
    Read pLink report dir and return an xtabel data array.

    Args:
        plinkdirs (list): plink report subdir (e.g. sample1)
        col_order (list) – List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Compact the xTable to only the columns given in col_order or not

    Returns:
        pandas.DataFrame: xTable data table
    """

    def plinkprotein2pandas(filepath):
        """
        Read a pLink protein results file and return a pandas dictionary.
    
        Args:
            filepath (str): Path to a pLink results file e.g. _inter_combine.protein.xls
    
        Returns:
            pandas.DataFrame
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
                        
            try:
                data = pd.DataFrame.from_dict(data)
                if len(data) > 0:
                    return data
                else:
                    raise Exception('Generated xtable had a length of 0!')
            except:
                raise Exception('Could not generat xtable. Please check file at: {}'.format(filepath))
    
    def read_plink_modifications(filepath):
        """
        Open a pLink modification.ini file and extract all modifications with
        their names as dict.
        
        Args:
            filepath (str): Path to modifications.ini
        
        Returns
            dict: mod_dict mapping pLink modification names to masses
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
    
        Args:
            seq_string (str): pLink sequence string
        Returns:
            list or np.nan: [pepseq1, pepseq2, xpos1, xpos2, xtype]
        """
        pattern = re.compile('(\w+)\((\d+)\)-(\w+)\((\d+)\):(\d+)')
        try:
            match = pattern.match(seq_string)
            pepseq1, xpos1, pepseq2, xpos2, xtype = match.groups()
            return str(pepseq1), int(xpos1), str(pepseq2), int(xpos2), xtype
    
        except Exception as e:
            print(e)
            return np.nan

    def process_plink_spectrum(spec_string):
        """
        Extract rawfile name, precursor charge and scan no from pLink spectrum
        string such as 2017_08_04_SVs_BS3_16.17079.17079.4.dta
    
        Args:
            spec_string: pLink spectrum string
    
        Returns:
            list or np.nan: [rawfile, scanno, prec_ch]
        """
        # the pattern of the title string is 20171215_JB04_Sec06.10959.10959.2
        # in pLink 2.3 and 20171215_JB04_Sec06.10959.10959.2.0 in pLink 2.1
        pextract_pattern = re.compile('(.+?)\.\d+\.(\d+)\.(\d+)\.*\d*')
        if pextract_pattern.match(spec_string):
            match = pextract_pattern.match(spec_string)
            rawfile, scanno, prec_ch = match.groups()
            return str(rawfile), int(scanno), int(prec_ch)
        else:
            return np.nan
    
    def process_plink_proteins(prot_string):
        """
        Extract protein name and absolute cross-link position from
        pLink protein string e.g.
        sp|P63045|VAMP2_RAT(79)-sp|P63045|VAMP2_RAT(59)
        
        Args:
            prot_string: pLink protein string
        
        Returns:
            list or np.nan: [prot1, xpos1, prot2, xpos2]
        """
        pattern = re.compile('(.+?)\((\d+)\)-?([^\(]*)\(?(\d*)\)?')
        match = pattern.match(prot_string)
        prot1, xpos1, prot2, xpos2 = match.groups()
        return str(prot1), int(xpos1), str(prot2), int(xpos2)


    ### Collect data, convert to pandas format and merge

    print('[pLink1 Read] This is pLink1 Reader')

    # convert to list if the input is only a single path
    if not isinstance(plinkdirs, list):
        plinkdirs = [plinkdirs]
    
    allData = list()

    plink_dtypes = {'Spectrum': str,
                    'Sequence': str,
                    'Proteins': str,
                    'type': str,
                    'Score': float,
                    'Order': int,
                    'Order2': int}

    for file in plinkdirs:

        # Initialise file names as None to use implicit booleaness
        inter_file = None
        loop_file = None
        mono_file = None

        for e in os.listdir(hf.FSCompatiblePath(file)):
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
            inter_df = plinkprotein2pandas(hf.FSCompatiblePath(os.path.join(file, inter_file)))
            inter_df['type'] = 'inter'
            frames.append(inter_df)
        if loop_file:
            print('Reading pLink loop-file: ' + loop_file)
            loop_df = plinkprotein2pandas(hf.FSCompatiblePath(os.path.join(file, loop_file)))
            loop_df['type'] = 'loop'
            frames.append(loop_df)
        if mono_file:
            print('Reading pLink mono-file: ' + mono_file)
            mono_df =  plinkprotein2pandas(hf.FSCompatiblePath(os.path.join(file, mono_file)))
            mono_df['type'] = 'mono'
            frames.append(mono_df)
    
        s = pd.concat(frames)

        allData.append(s)

    xtable = pd.concat(allData).astype(dtype=plink_dtypes)
    ### Convert data inside pandas df

    # rawfile, scanno, prec_ch
    xtable[['rawfile', 'scanno', 'prec_ch']] =\
        pd.DataFrame(xtable['Spectrum'].apply(process_plink_spectrum).tolist(), index=xtable.index)

    # Directly assign the re group matches into new columns
    xtable[['pepseq1', 'xlink1', 'pepseq2', 'xlink2', 'xtype']] =\
        pd.DataFrame(xtable['Sequence'].apply(process_plink_sequence).tolist(), index=xtable.index)

    xtable[['prot1', 'xpos1', 'prot2', 'xpos2']] =\
            pd.DataFrame(xtable['Proteins'].apply(process_plink_proteins).tolist(), index=xtable.index)

    xtable['score'] = xtable['Score']

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        pd.Series(np.vectorize(hf.generateID,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                 index=xtable.index).replace('nan', np.nan)

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    xtable['pos1'] = xtable['xpos1'].astype(int, errors='ignore') - \
                     xtable['xlink1'].astype(int, errors='ignore') + 1
    xtable['pos2'] = xtable['xpos2'].astype(int, errors='ignore') - \
                     xtable['xlink2'].astype(int, errors='ignore') + 1

    # add a lobel referring to the ordering in the pLink results table
    xtable['Order'] = xtable[['Order', 'Order2']].apply(lambda x: ','.join(str(x)), axis=1)

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for intra and inter xlink to inter/intra/homomultimeric
        intraAndInter = (xtable['type'] == 'inter') | (xtable['type'] == 'intra')
        xtable.loc[intraAndInter, 'type'] =\
            np.vectorize(hf.categorizeInterPeptides)(xtable[intraAndInter]['prot1'],
                                                     xtable[intraAndInter]['pos1'],
                                                     xtable[intraAndInter]['pepseq1'],
                                                     xtable[intraAndInter]['prot2'],
                                                     xtable[intraAndInter]['pos2'],
                                                     xtable[intraAndInter]['pepseq1'])
        print('[xQuest Read] categorized inter peptides')
    else:
        print('[xQuest Read] skipped inter peptide categorization')

    # manually set decoy to reverse as pLink hat its own internal target-decoy
    # algorithm
    xtable['decoy'] = False

    # generate the mod_dict linking pLink modification names to masses
    
    # in case of calling croco from the source folder structure...
    file_dir, file_name = os.path.split(__file__)
    if os.path.exists(os.path.join(file_dir,
                                   '../data/modification.ini')):
        modifi_dir = os.path.abspath(os.path.join(file_dir,
                                                  '../data/modification.ini'))
    # ... or calling from a exe-file in a folder-setup with the data folder at top-level
    elif os.path.exists(os.path.join(file_dir,
                                     './data/modification.ini')):
        modifi_dir = os.path.abspath(os.path.join(file_dir,
                                                  './data/modification.ini'))
    # ... or calling from within a single bundled exe-file
    else:
        try:
            # PyInstaller creates a temp folder and stores its path in _MEIPASS
            base_path = sys._MEIPASS
            modifi_dir =  os.path.abspath(\
                os.path.join(base_path, './data/modification.ini'))
        # ... or something went wrong
        except:
            raise Exception('Modifications.ini not found. CWD is ' + file_dir)

    # load pLink modifications.ini from data-folder
    mod_dict = read_plink_modifications(os.path.abspath(modifi_dir))

    # extract modification information
    pattern = re.compile(r'(\d+),.*\((.*)\)')

    pepseq1 = xtable['pepseq1'].tolist()
    Modification = xtable['Modification'].tolist()

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
                
                seqlen1 = len(pepseq1[idx])
                # pLink assigns additional modification position to the C-term
                # of the first peptide, the xlinker and the N-term of the
                # second peptide
                if int(modpos) > seqlen1:
                    this_mod2.append(mod)
                    this_modpos2.append(int(modpos) - seqlen1)
                    this_modmass2.append(mass)
                else:
                    this_mod1.append(mod)
                    this_modpos1.append(modpos)
                    this_modmass1.append(mass)

        # multiple modifications of one peptide are stored as ;-delimited strings
        modmass1.append(this_modmass1)
        modpos1.append(this_modpos1)
        mod1.append(this_mod1)
        modmass2.append(this_modmass2)
        modpos2.append(this_modpos2)
        mod2.append(this_mod2)

    xtable['mod1'] = mod1
    xtable['modmass1'] = modmass1
    xtable['modpos1'] = modpos1
    xtable['mod2'] = mod2
    xtable['modmass2'] = modmass2
    xtable['modpos2'] = modpos2

    xtable['search_engine'] = 'pLink1'

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
        
    xtable = Read(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\pLink1_results\2.report\sample1', col_order=col_order)