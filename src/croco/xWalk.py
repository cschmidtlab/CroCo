# -*- coding: utf-8 -*-

"""
Functions to write xWalk data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd
import numpy as np
import os

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def AA_from_sequence(pepseq, xlink):
    """
    Return the 3-character amino acid label of the cross-linked AA
    from a peptide sequence
    and the relative position of the cross-linker in the sequence
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


def Write(xtable, outpath, pdb, offset, chains, atom):
    """
    Converts xTable into a list format that can be used as
    input for the xWalk web-server at http://www.xwalk.org/cgi-bin/home.cgi

    Format is:

    index \t pdb-file \t RESIDUE-NO--ATOM \t RESIDUE-NO--ATOM

    As xWalk can only validate one protein at a time, the function
    generates several oouput files for all intra-protein cross-links

    Args:
        xtable: data table structure
        pdb (string): PDB-file name
        offset (int): shift between PDB AA indices and the xTable
        chains: (string) comma separated list protein:chain allocations
        atom (string): Atom identifier (e.g. CB)
        outpath: path to write file
    """


    print(outpath, pdb, offset, chains, atom)

    pdbBase = os.path.basename(pdb)

    if not pdbBase.endswith('.pdb'):
        raise Exception('Please provide a valid PDB file')

    xtable['File name'] = pdbBase

    if len(atom.strip()) > 4:
        raise Exception('[xWalk Write] Please provide PDB atom code with at most 4 characters')

    xtable['atom'] = str(atom).upper()

    chainDict = dict()
    try:
        for annotation in chains.split(','):
            protein, chain = annotation.strip().split(':')
            # by generating a list of the string, all characters will be represented
            # as single chain identifiers
            chainDict[protein] = list(chain.upper())
    except:
        raise Exception('[xWalk Write] Please specify protein:chain in an comma-separated list')

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

    # convert the offset user-input into an integer as requried for pandas below
    try:
        offset = int(offset)
    except:
        raise Exception('[xWalk Write] Please provide an integer offset!')

    xWalkTable['xpos1'] += offset
    xWalkTable['xpos2'] += offset

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
    out = r'C:\Users\User\Downloads\test'
    chains = 'SPA_STAAU:AD, IgG4_heavy:BC'

    xtable = Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')

    xtable = Write(xtable=xtable, outpath=out, pdb=pdb, offset=0, chains=chains, atom=atom)