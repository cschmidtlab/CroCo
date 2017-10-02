
"""
Functions to read-in and write-out cross-link information from different sources.

=== Definition of the xtable format ===

rawfile - name of the corresponding rawfile
scanno - Scan the cross-link was identified from
prec_ch - charge of the precursor ion

prot1 - fasta-header of the protein corresponding to the first peptide
prot2 - same for the second peptide

pepseq1 - Sequence of the longer (alpha) peptide
pepseq2 - Sequence of the shorter (beta) peptide

xlink1 - Position of the cross-linker within the longer peptide
xlink2 - Position of the cross-linker within the shorter peptide

xpos1 - Absolute position of the cross-linker of the longer peptide
xpos2 - Absolute position of the cross-linker of the shorter peptide (only if interlink)

# modpos1 - relative position of a modification within peptide 1
# modpos2 - same for peptide 2
mod1 - type of modification 1
mod2 - type of modification 2

ID - Identifier for the position of a cross-link between two proteins
decoy - true or false

score - not normalised score
type - inter, intra, loop or mono
xtype - heavy or light label

=== Definition of the xinfo format ===

xlinker - Name of the used reagent
link_dist - nominally bridged distance in Angstrom

"""

import numpy as np
import pandas as pd

import os
import re

def plinkspectra2pandas(filepath):
    """
    Read a pLink spectra results file and return a pandas dictionary
    
    :params: filepath: Path to a pLink results file e.g. _inter_combine.spectrum.xls
    
    :returns: pandas dataframe
    """
    
    with open(filepath, 'r') as fh:

        # read the first header-line into list
        headers1 = fh.readline().strip().split('\t')
        # read the second header-line and remove it from list
        # avoid the first entry as it is only the line-indicator
        headers2 = fh.readline().strip().split('\t')[1:]
        headers2 = [x if x != 'Score' else 'Score2' for x in headers2 ]
        
        data = {} # init of data dict for pandas
        for h in headers1 + headers2:
            data[h] = []

        # read the first line
        line = fh.readline()
        while line:
                    
            if line.startswith('*'): # indicates a headers2 line
                line_data = line.strip().split('\t')[1:] # remove the first element
                for i in range(len(line_data)):
                    # iterate through the data and link them to their resp
                    # header by index
                    data[headers2[i]].append(line_data[i])
            else:
                line_data = line.strip().split('\t')
                for i in range(len(line_data)):
                    data[headers1[i]].append(line_data[i])
    
            # read the next line
            line=fh.readline()
        
    return pd.DataFrame.from_dict(data)

def plinkpeptide2pandas(filepath):
    """
    Read a pLink peptide results file and return a pandas dictionary
    
    :params: filepath: Path to a pLink results file e.g. _inter_combine.peptide.xls
    
    :returns: pandas dataframe
    """
    
    with open(filepath, 'r') as fh:

        # read the first header-line into list
        headers1 = fh.readline().strip().split('\t')
        # read the second header-line and remove it from list
        # avoid the first entry as it is only the line-indicator
        headers2 = fh.readline().strip().split('\t')
        headers2 = [x if x != 'Order' else 'Order2' for x in headers2 ]
        
        data = {} # init of data dict for pandas
        for h in headers1 + headers2:
            data[h] = []

        # read the first line
        line = fh.readline()
        while line:
                    
            if line[0].isdigit(): # indicates a headers1 line
                # save the line and use it when printing all following lines
                # corresponding to that title-line
                line1_data = line.strip().split('\t')

            else:
                line2_data = line.strip().split('\t')
                for i in range(len(line1_data)):
                    # iterate through the data and link them to their resp
                    # header by index
                    data[headers1[i]].append(line1_data[i])
                for i in range(len(line2_data)):
                    data[headers2[i]].append(line2_data[i])
    
            # read the next line
            line=fh.readline()
        
    return pd.DataFrame.from_dict(data)


def process_sequence(seq_string):
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

def process_spectrum(spec_string):
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

def process_proteins(prot_string):
    """
    Extract protein name and absolute cross-link position from
    pLink protein string e.g.
    sp|P63045|VAMP2_RAT(79)-sp|P63045|VAMP2_RAT(59)
    """
    pattern = re.compile('(.+?)\((\d+)\)-?([^\(]*)\(?(\d*)\)?')
    match = pattern.match(prot_string)
    return match.groups()

def process_kojak_peptide(peptide_string):
    """
    Return Modifications and the peptide sequence
    from a Kojak sequence string such as M[15.99]TDSKYFTTNK
    """
    
    pattern = re.compile('\[(.*?)\]')
    mods = re.findall(pattern, peptide_string)
    # mod_pos = []
    # for n,i in enumerate(mods):
    #     # calculate mod position by finding index in string
    #     # and subtracting the amount of brackets up to the
    #     # nth modification within the peptide
    #     mod_pos.append(peptide_string.index(i) - len(mod_pos[-1:]) - (1 + 2*n))

    pattern = re.compile('([A-Z]+)')
    sequence = ''.join(re.findall(pattern, peptide_string))

    if mods == []:
        mods = np.nan
    #     mod_pos = np.nan
    
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
        return match.groups()
    else:
        return np.nan, np.nan

def ReadpLink(plinkdir):
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
        if '_inter_combine.peptide' in e:
            inter_file = e
            
        if '_loop_combine.peptide' in e:
            loop_file = e
        
        if '_mono_combine.peptide' in e:
            mono_file = e
    
    frames = []
    # only called if inter_file is not None
    if inter_file:
        print('Reading pLink inter-file: ' + inter_file)
        inter_df = plinkpeptide2pandas(os.path.join(plinkdir, inter_file))
        inter_df['type'] = 'inter'
        frames.append(inter_df)
    if loop_file:
        print('Reading pLink loop-file: ' + loop_file)
        loop_df = plinkpeptide2pandas(os.path.join(plinkdir, loop_file))
        loop_df['type'] = 'loop'
        frames.append(loop_df)
    if mono_file:
        print('Reading pLink mono-file: ' + mono_file)  
        mono_df =  plinkpeptide2pandas(os.path.join(plinkdir, mono_file))
        mono_df['type'] = 'mono'
        frames.append(mono_df)
    
    data = pd.concat(frames)
        
    ### Convert data inside pandas df
    

    # init xtable with column containing lists of rawfile, scanno, prec_ch
    xtable = pd.DataFrame(data['Spectrum'].apply(process_spectrum))
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
        zip(*data['Sequence'].apply(process_sequence))
    
    xtable['prot1'], xtable['xpos1'], xtable['prot2'], xtable['xpos2'] =\
            zip(*data['Proteins'].apply(process_proteins))

    xtable['type'] = data['type']
    xtable['score'] = data['Score']
    
    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # manually set decoy to reverse as pLink hat its own internal target-decoy
    # algorithm
    xtable['decoy'] = False

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')
    
    # compute a minus log P score for better comparison with higher=better scores
    xtable['-log score'] = -np.log(xtable['score'])


    ### return xtable df
    
    return xtable


def ReadKojak(kojakpath, rawfile):
    """
    reads pLink results file and returns an xtable data array.
    
    :params: kojakpath: Kojak results file
    :params: rawfile: name of the corresponding rawfile

    :returns: xtable data table
    :returns: xinfo meta-data object

    """
    
    ### Collect data and convert to pandas format

    for e in os.listdir(kojakpath):
        if '.kojak.txt' in e:
            kojak_file = e
    
    print('Reading Kojak-file: ' + kojak_file)

    # only called if inter_file is not None
    if kojak_file:
        data = pd.read_csv(os.path.join(kojakpath, kojak_file),
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

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    ### return xtable df
    
    return xtable

def ReadManual(filepath):
    """
    Read Excel file with hand-curated results
    """
    
    cols = ['order1',
            'order2',
            'Spectrum',
            'Sequence',
            'is_xlink',
            'Score',
            'Calc_M',
            'Delta_M',
            'ppm',
            'Modification',
            'Sample',
            'Engine',
            'MatchedIons',
            'MissCleaveNum',
            'Rank',
            'Proteins']
            
    
    true_hits = pd.read_excel(filepath,
                            sheetname='yes',
                            header=None,
                            names=cols)
                            
    true_hits.dropna(axis=0,
                how='all',
                inplace=True)
                
    true_hits['score'] = 1 # set true hits to a score of 1 and maybes to 0.5
                            
    maybe_hits = pd.read_excel(filepath,
                            sheetname='naja',
                            header=None,
                            names=cols)
    
    maybe_hits.dropna(axis=0,
                how='all',
                inplace=True)
                
    maybe_hits['score'] = 0.5
    
    no_hits = pd.read_excel(filepath,
                            sheetname='no',
                            header=None,
                            names=cols)
    
    no_hits.dropna(axis=0,
                how='all',
                inplace=True)
                
    no_hits['score'] = 0.0
    
    
    ann_data = pd.concat((true_hits, maybe_hits, no_hits))   
    
    ann_data.reset_index(drop=True, inplace=True)
    
    ann_data['pepseq1'], ann_data['xlink1'], ann_data['pepseq2'],\
    ann_data['xlink2'], ann_data['xtype'] =\
        zip(*ann_data['Sequence'].apply(process_sequence))
    
    ann_data['prot1'], ann_data['xpos1'], ann_data['prot2'], ann_data['xpos2'] =\
            zip(*ann_data['Proteins'].apply(process_proteins))
                
    ann_data['rawfile'], ann_data['scanno'], ann_data['prec_ch'] =\
        zip(*ann_data['Spectrum'].apply(process_spectrum))

    # generate an ID for every crosslink position within the protein(s)
    ann_data['ID'] =\
        ann_data[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # manually set reverse status to false
    ann_data['decoy'] = False

    ann_data = ann_data.apply(pd.to_numeric, errors = 'ignore')

    return ann_data

def WritePercolator(xtable, outdir):
    """
    writes an xtable data structure to file (in Percolator format)
    
    :params: xtable: data table structure
    :params: outdir: filepath to write Kojak-formatted file to
    
    :returns: True if successfull
    """

def WriteXtable(xtable, usecols, outdir):
    """
    writes an xtable data structure to file (in xtable format)
    
    :params: xtable: data table structure
    :params: usecols: columns to use for export, leave empty for all
    :params: outdir: filepath to write file to
    
    :returns: True if successfull
    """
    
    