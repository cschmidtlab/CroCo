# -*- coding: utf-8 -*-

"""
Functions to write xWalk data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd
import numpy as np
import os

from . import HelperFunctions as HeFn

def xWalk_from_xTable(pepseq, xlink, xpos, atom):
    """
    Calculate a string like LYS-367--CB from xtable columns
    
    Args:
        pepseq: a peptide sequence
        xlink: position of lxink within peptide sequence
        xpos: absolute position of xlink in protein
        atom: atoms to indicate xWalk to link (e.g. CB)
        
    Returns:
        xWalk-string
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
    
    return AA + '-' + str(xpos) + '--' + str(atom)

def Write(xtable, outpath, pdb, atom):
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
        atom (string): Atom identifier (e.g. CB)
        outpath: path to write file  
    """
    
    pdb = os.path.basename(pdb)
    
    if not pdb.endswith('.pdb'):
        raise Exception('Please provide a valid PDB file')
    
    # filter intra-protein cross-links
    xtable = xtable[xtable['prot1'] == xtable['prot2']].copy()

    xtable['pdb'] = pdb
    xtable['atom'] = atom

    for prot in set(xtable['prot1']):
        thisXtable = xtable[xtable['prot1'] == prot]\
            .loc[:, ['pepseq1',
                     'pepseq2',
                     'xpos1',
                     'xpos2',
                     'xlink1',
                     'xlink2',
                     'prot1',
                     'prot2',
                     'pdb',
                     'atom']].copy()
        thisXtable[['xpos1', 'xpos2']].drop_duplicates()
                
        thisXtable['xWalk1'] = np.vectorize(xWalk_from_xTable)\
            (thisXtable['pepseq1'],
             thisXtable['xlink1'],
             thisXtable['xpos1'],
             thisXtable['atom'])
        
        thisXtable['xWalk2'] = np.vectorize(xWalk_from_xTable)\
            (thisXtable['pepseq2'],
             thisXtable['xlink2'],
             thisXtable['xpos2'],
             thisXtable['atom'])

        outpath = os.path.splitext(outpath)[0]

        thisXtable.loc[:, ['pdb', 'xWalk1', 'xWalk2']]\
            .to_csv('{}_{}.tsv'.format(outpath, HeFn.alphanum_string(prot)),
                                       index = True,
                                       header = False,
                                       sep='\t')

if __name__ == '__main__':
    import xTable
    
    pdb = '1aqf.pdb'
    atom = 'CB'
    out = r'C:\Users\User\Downloads\test.tsv'
    
    xtable = xTable.Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')
    
    Write(xtable, out, pdb, atom)