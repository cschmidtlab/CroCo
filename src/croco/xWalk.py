# -*- coding: utf-8 -*-

"""
Functions to write data as input for xWalk.

"""

import pandas as pd
import numpy as np
import os

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath, pdb, offset, chains, atom):
    """
    Convert xTable into a list format that can be used as
    input for the xWalk standalone programme.

    Format is:

    index \t pdb-file \t RESIDUE-NO--ATOM \t RESIDUE-NO--ATOM

    As xWalk can only validate one protein at a time, the function
    generates several oouput files for all intra-protein cross-links

    Args:
        xtable (pandas.DataFrame): data table structure
        pdb (str): PDB-file name
        offset (list or str): shift between PDB AA indices and the xTable
        chains: (dict or str) comma separated list protein:chain allocations
        atom (str): Atom identifier (e.g. CB)
        outpath (str): path to write file
    """

    
    def AA_from_sequence(pepseq, xlink):
        """
        Return the 3-character amino acid label of the cross-linked AA
        from a peptide sequence
        and the relative position of the cross-linker in the sequence
        
        Args:
            pepseq (str): peptide sequence
            xlink (int): position of the cross-link within the sequence
        Returns:
            str: 3-letter amino acid code for the cross-linked amino acid
        """
    
        aa_dict = {'R': 'ARG',
               'H': 'HIS',
               'K': 'LYS',
               'D': 'ASP',
               'E': 'GLU',
               'S': 'SER',
               'T': 'THR',
               'N': 'ASN',
               'Q': 'GLN',
               'C': 'CYS',
               'U': 'SEC',
               'G': 'GLY',
               'P': 'PRO',
               'A': 'ALA',
               'V': 'VAL',
               'I': 'ILE',
               'L': 'LEU',
               'M': 'MET',
               'F': 'PHE',
               'Y': 'TYR',
               'W': 'TRP'
               }
    
        AA = aa_dict[pepseq[int(xlink)-1].upper()]
    
        return AA

    pdbBase = os.path.basename(pdb)

    if not pdbBase.endswith('.pdb'):
        raise Exception('Please provide a valid PDB file')

    xtable['File name'] = pdbBase

    if len(atom.strip()) > 4:
        raise Exception('[xWalk Write] Please provide PDB atom code with at most 4 characters')

    xtable['atom'] = str(atom).upper()

    chainDict = dict()
    
    if not isinstance(chains, dict):
        if isinstance(chains, str):
            chains = [x.strip() for x in chains.split(',')]
        try:
            for annotation in chains:
                protein, chain = annotation.strip().split(':')
                # by generating a list of the string, all characters will be represented
                # as single chain identifiers
                chainDict[protein] = list(chain.upper())
        except:
            raise Exception('[xWalk Write] Please specify protein:chain in an comma-separated list from the GUI or as a dict')

    # drop duplicates on the cross-link position as only the absolute position
    # is relevant to xWalk
    xtable.drop_duplicates(subset=['xpos1', 'xpos2'], keep='first', inplace=True)

    # remove rows that contain NaN in prot1 or prot2 i.e. monolinks
    xtable.dropna(subset=['prot1', 'prot2'], inplace=True)

    # set the 3-character code for the cross-linked amino acids
    xtable['linked_aa1'] = np.vectorize(AA_from_sequence)\
        (xtable['pepseq1'],
         xtable['xlink1'])

    xtable['linked_aa2'] = np.vectorize(AA_from_sequence)\
        (xtable['pepseq2'],
         xtable['xlink2'])

    allChainTables = list()

    for proteinA in chainDict.keys():
        for chainA in chainDict[proteinA]:
            for proteinB in chainDict.keys():
                for chainB in chainDict[proteinB]:
                    thisXTable = xtable[(xtable['prot1'] == proteinA) &\
                                        (xtable['prot2'] == proteinB)][['File name',
                                                                        'atom',
                                                                        'pepseq1',
                                                                        'pepseq2',
                                                                        'xpos1',
                                                                        'xpos2',
                                                                        'prot1',
                                                                        'prot2',
                                                                        'linked_aa1',
                                                                        'linked_aa2']]
                    thisXTable['chain1'] = chainA
                    thisXTable['chain2'] = chainB

                    allChainTables.append(thisXTable)

    xWalkTable = pd.concat(allChainTables)

    # to assign offsets to every protein, a single integer (one for all) or a 
    # dict mapping protein names to offsets is required
    if not isinstance(offset, dict):
        # convert the offset user-input into an integer as requried for pandas below
        try:
            # the input is a single integer
            offset = int(offset)
        except:
            # the input is a list of protein:offset pair strings
            try:
                offset = [x.strip() for x in offset.split(',')]
            except:
                raise Exception('[xWalk Write] Please provide an integer offset for all chains or a list of protein:offset assignments!')
    
    # if the offset is an integer, use it for all protein positions
    if isinstance(offset, int):
        xWalkTable['xpos1'] += offset
        xWalkTable['xpos2'] += offset
    else:
        # if it is a list (see above) the protein:offset pairs are parted
        if isinstance(offset, list):
            offsetDict = dict()
            try:
                for annotation in offset:
                    protein, offset = annotation.strip().split(':')
                    offsetDict[protein] = int(offset)
            except:
                raise Exception('[xWalk Write] Please specify protein:offset in an comma-separated list from the GUI or as a dict')
        # if it si a dict, it can directly be used
        if isinstance(offset, dict):
            offsetDict = offset
    
        try:
            print(offsetDict)
            for pr, of in offsetDict.items():
                xWalkTable.loc[xWalkTable['prot1'] == pr, 'xpos1'] += of
                xWalkTable.loc[xWalkTable['prot2'] == pr, 'xpos2'] += of
        except:
            raise Exception('[xWalk Write] error during assignment of offsets to proteins')

    atomInfo1 = list()
    atomInfo2 = list()

    for idx, row in xWalkTable.iterrows():
        atomInfo1.append('-'.join([str(row['linked_aa1']),
                                   str(int(row['xpos1'])),
                                   str(row['chain1']),
                                   str(row['atom'])]))

        atomInfo2.append('-'.join([str(row['linked_aa2']),
                                   str(int(row['xpos2'])),
                                   str(row['chain2']),
                                   str(row['atom'])]))

    xWalkTable['Atom Info 1'] = atomInfo1
    xWalkTable['Atom Info 2'] = atomInfo2

    # Remove those amino acids interacting with itself (distance = 0)
    xWalkTable = xWalkTable[xWalkTable['Atom Info 1'] != xWalkTable['Atom Info 2']]

    xWalkTable.reset_index(inplace=True)

    xWalkTable.loc[:, ['File name', 'Atom Info 1', 'Atom Info 2']]\
        .to_csv('{}_{}.tsv'.format(hf.FSCompatiblePath(outpath), 'xWalk'),
                                   index = True,
                                   index_label = 'Index',
                                   header = True,
                                   sep='\t')

if __name__ == '__main__':
    from xTable import Read

    pdb = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\xWalk\1aqf.pdb'
    atom = 'CB'
    out = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\xWalk\xWalk'
    chains = 'SPA_STAAU:AD, IgG4_heavy:BC'
    offset = 'SPA_STAAU:1, IgG4_heavy:0'

    xtable = Read(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\xWalk\pLink1_xtable_xTable_to_xTable.xlsx')

    xtable = Write(xtable=xtable, outpath=out, pdb=pdb, offset=offset, chains=chains, atom=atom)