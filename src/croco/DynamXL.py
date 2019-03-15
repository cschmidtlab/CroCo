# -*- coding: utf-8 -*-

"""
    Converts xTable data into cross-link information
    file for the dynamXL analysis software:

    K189 NZ K192 NZ 18
    K189 NZ K196 NZ 5
    K200 NZ K192 NZ 1
    K236 NZ K189 NZ 3

    dynamxl.chem.ox.ac.uk

"""

import pandas as pd
import numpy as np

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath):
    """
    Convert xTable to DynamXL input file.

    Args:
        xtable: data table structure
        outpath: path to write file
    """
    def aa_and_pos(row):
        """
        Calculate tne amino acid and the absolute position
        of a crosslink based on the sequence of an xl-peptide,
        the relative and the absolute position of the cross-link
        
        Args:
            row (object): xTable row containing pepseq, xlink, and xpos
        Returns:
            id1 (str)
            id2 (str)
        """
        # prevent error from calculation of AA from entries without xlink
        if hf.isNaN(row['xlink1']):
            id1 = np.nan
        else:
            aa1 = row['pepseq1'][int(row['xlink1'])-1]
            id1 = str(aa1) + str(int(row['xpos1']))

        if hf.isNaN(row['xlink2']):
            id2 = np.nan
        else:
            aa2 = row['pepseq2'][int(row['xlink2'])-1]
            id2 = str(aa2) + str(int(row['xpos2']))

        return id1, id2

    def xlink_atom_from_AA(ID):
        """
        Return the typical cross-linked atom in PDB code
        for a specific amino-acid
        
        Args:
            ID (str): DynamXL ID e.g. K27
        Returns:
            PDB code of cross-linked atom
        """

        if hf.isNaN(ID):
            return np.nan
        else:
            if ID[0] == 'K':
                return 'NZ'
            else:
                return 'CA'

    print('Converting to dynamXL input file format')

    # init xtable with column containing lists of rawfile, scanno, prec_ch
    dynamxl = pd.DataFrame(xtable['score'])

    dynamxl['ID1'], dynamxl['ID2'] = zip(*xtable.apply(aa_and_pos, axis=1))

    dynamxl['atom1'] = dynamxl['ID1'].apply(xlink_atom_from_AA)

    dynamxl['atom2'] = dynamxl['ID2'].apply(xlink_atom_from_AA)

    # reorder df
    dynamxl = dynamxl[['ID1', 'atom1', 'ID2', 'atom2', 'score']]

    dynamxl.to_csv(hf.FSCompatiblePath(outpath + '.txt'),
                   sep = '\t',
                   header=False,
                   float_format='%.3f',
                   index=False)