
"""
Functions to read-in cross-link information from different sources.

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

pos1 - Absolute postiton of the first AA of the first peptide
pos2 - Absolute postiton of the first AA of the second peptide

xpos1 - Absolute position of the cross-linker of the longer peptide
xpos2 - Absolute position of the cross-linker of the shorter peptide (only if interlink)

mod1 - relative position of a modification within peptide 1 (;-delimited string)
mod2 - same for peptide 2
modmass1 - mass of modification 1 (;-delimited string)
modmass2 - mass of modification 2

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

# defines the column headers required for xtable output
col_order = [ 'rawfile', 'scanno', 'prec_ch',
              'pepseq1', 'xlink1',
              'pepseq2', 'xlink2', 'xtype', 'modmass1', 'mod1', 'modmass2', 'mod2',
              'prot1', 'xpos1', 'prot2',
              'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

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
        prot, xpos = match.groups()
        if xpos == None: # re.match returns None (not NaN) if a substring doesnt match
            return prot, np.nan
        else:
            return prot, xpos
    else:
        return np.nan, np.nan


def process_xQuest_spectrum(spec_string):
    """
    Extract rawfile name, precursor charge and scan no from xQuest sequence
    string

    :returns: list of rawfile, scanno, prec_ch
    """
    spectrum_pattern = re.compile('(.+)\.(\d+)\.\d+\..+\.\d+\.\d+\.(\d+)')
    if spectrum_pattern.match(spec_string):
        match = spectrum_pattern.match(spec_string)
        # rawfile, scanno, prec_ch
        return match.groups()
    else:
        return np.nan

def process_xQuest_Id(Id_string):
    """
    Extract peptide sequence of the alpha (longer) and the beta (shorter)
    peptide as well as the relative positions of the cross-links within
    these sequences from an xQuest Id-string
    """
    Id_pattern = re.compile('(\w+)-(\w+)-a(\d+)-b(\d+)')
    if Id_pattern.match(Id_string):
        match = Id_pattern.match(Id_string)
        # pepseq1, pepseq2, xlink1, xlink2
        return match.groups()
    else:
        return np.nan

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
    mod_dict = read_plink_modifications('../data/config/pLink/modification.ini')

    # extract modification information
    pattern = re.compile(r'(\d+),.*\((.*)\)')

    pepseq1 = xtable['pepseq1'].tolist()
    Modification = data['Modification'].tolist()

    if len(pepseq1) == len(Modification):
        print('Len of pepseq1 and Modification match!')
    else:
        print('Len of pepseq1 and Modification dont match!')

    modmass1 = []
    mod1 = []
    modmass2 = []
    mod2 = []

    # iterate over all lines in the input file
    for idx, modstr in enumerate(Modification):

        this_modmass1 = []
        this_mod1 = []
        this_modmass2 = []
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
                    this_mod2.append(str(int(modpos) - seqlen1))
                    this_modmass2.append(mass)
                else:
                    this_mod1.append(modpos)
                    this_modmass1.append(mass)

        # multiple modifications of one peptide are stored as ;-delimited strings
        modmass1.append(';'.join(this_modmass1))
        mod1.append(';'.join(this_mod1))
        modmass2.append(';'.join(this_modmass2))
        mod2.append(';'.join(this_mod2))

    xtable['modmass1'] = modmass1
    xtable['mod1'] = mod1
    xtable['modmass2'] = modmass2
    xtable['mod2'] = mod2

    # reorder columns
    # append pLink specific columns
    ['Order'] + col_order # append to front
    col_order.extend(['PSM image', 'plink score'])
    xtable = xtable[col_order]

    ### return xtable df

    return xtable


def ReadKojak(kojak_file, rawfile=None):
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
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # calculate absolute position of first AA of peptide
    # ignoring errors avoids raising error in case on NaN -> returns NaN
    # as pos
    # Must be calculated as float as NaN is not implemented in int
    xtable['pos1'] = xtable['xpos1'].astype(float, errors='ignore') - \
                     xtable['xlink1'].astype(float, errors='ignore') + 1
    xtable['pos2'] = xtable['xpos2'].astype(float, errors='ignore') - \
                     xtable['xlink2'].astype(float, errors='ignore') + 1

    # set the rawfile name for xtable (None if not provided by call)
    xtable['rawfile'] = rawfile

    xtable['xtype'] = np.nan

    # reassign dtypes for every element in the df
    # errors ignore leaves the dtype as object for every
    # non-numeric element
    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    # reorder columns
    xtable = xtable[col_order]

    ### return xtable df

    return xtable

def ReadxQuest(xquest_file):
    """
    Read xQuest results file and return file in xTable format.

    :params: kojakpath: Kojak results file
    :params: rawfile: name of the corresponding rawfile

    :returns: xtable data table
    """

    ### Collect data and convert to pandas format

    print('Reading xQuest-file: ' + xquest_file)

    # only called if inter_file is not None
    try:
        xquest = pd.read_csv(xquest_file,
                             delimiter='\t')
    except:
        raise(FileNotFoundError('xQuest results file not found'))

    rename_dict = {'z':'prec_ch',
                   'Protein1':'prot1',
                   'Protein2': 'prot2',
                   'AbsPos1': 'xpos1',
                   'AbsPos2': 'xpos2',
                   'Type': 'type',
                   'ld-Score': 'score'}

    # Copy and rename selected columns to new xquest df
    try:
        xtable = xquest.loc[:, list(rename_dict.keys())]
        xtable.rename(index=str,
                    columns=rename_dict,
                    inplace=True)
    except Exception as e:
        raise Exception('Error during xQuest header renaming: %s' % e)

    # Extract rawfile, scanno and precursor charge from the mgf header string
    # used as Spectrum by xQuest
    xtable['rawfile'], xtable['scanno'], xtable['prec_ch'] =\
        zip(*xquest['Spectrum'].apply(process_xQuest_spectrum))

    # Extract peptide sequences and relative cross-link positions form the
    # xQuest ID-string
    xtable['pepseq1'], xtable['pepseq2'], xtable['xlink1'], xtable['xlink2'] =\
        zip(*xquest['Id'].apply(process_xQuest_Id))

    # Modifications are not defined in xQuest
    xtable['mod1'], xtable['mod2'] = "", ""

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        xtable[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # xQuest does not incorporate decoy entries in the resutls table
    # but protein names can contain identifiers as reverse or decoy
    xtable['decoy'] = xtable['ID'].str.contains('reverse') |\
        xtable['ID'].str.contains('decoy')

    ### Return df
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
        zip(*ann_data['Sequence'].apply(process_plink_sequence))

    ann_data['prot1'], ann_data['xpos1'], ann_data['prot2'], ann_data['xpos2'] =\
            zip(*ann_data['Proteins'].apply(process_plink_proteins))

    ann_data['rawfile'], ann_data['scanno'], ann_data['prec_ch'] =\
        zip(*ann_data['Spectrum'].apply(process_plink_spectrum))

    # generate an ID for every crosslink position within the protein(s)
    ann_data['ID'] =\
        ann_data[['prot1', 'xpos1', 'prot2', 'xpos2']].astype(str).apply(\
            lambda x: '-'.join(x), axis=1)

    # manually set reverse status to false
    ann_data['decoy'] = False

    ann_data = ann_data.apply(pd.to_numeric, errors = 'ignore')

    return ann_data