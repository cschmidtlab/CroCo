
"""
Functions to write-out cross-link information from different sources.

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

import pandas as pd

def WritePercolator(xtable, outdir):
    """
    writes an xtable data structure to file (in Percolator format)
    
    :params: xtable: data table structure
    :params: outdir: filepath to write Kojak-formatted file to
    
    :returns: True if successfull
    """

def WriteXtable(xtable, outpath):
    """
    writes an xtable data structure to file (in xtable format)
    
    :params: xtable: data table structure
    :params: outpath to write file (w/o file extension!)
    """
    xtable.to_excel(outpath + '.xlsx',
                    index=False)

def WritexVis(xtable, outpath):
    """
    Converts xtable data structure to cross-link
    data file for xVis data 
    visualisation tool
    (https://xvis.genzentrum.lmu.de/CrossVisNoLogin.php)
    
    :params: xtable: data table structure
    :params: outpath: path to write file
    """
    xvis = xtable.loc[:,['prot1','prot2', 'xpos1', 'xpos2', 'score']]
    
    rename_dict = {'prot1':'Protein1',
                   'prot2':'Protein2',
                   'xpos1': 'AbsPos1',
                   'xpos2': 'AbsPos2',
                   'score': 'Id-Score'}
    xvis.rename(index=str,
                columns=rename_dict,
                inplace=True)
    
    xvis.to_csv(outpath + '.csv',
                index=False)

def WritexiNET(xtable, outpath):
    """
    Converts xtable data structure to cross-link
    data file for xiNET data 
    visualisation tool
    (http://crosslinkviewer.org)
    
    :params: xtable: data table structure
    :params: outpath: path to write file
    """
    xinet = xtable.loc[:,['prot1',
                          'pos1',
                          'pepseq1',
                          'xlink1',
                          'prot2',
                          'pos2',
                          'pepseq2',
                          'xlink2',
                          'score',
                          'ID']]
    
    rename_dict = {'prot1':'Protein1',
                   'prot2':'Protein2',
                   'pos1': 'PepPos1',
                   'pos2': 'PepPos2',
                   'pepseq1': 'PepSeq1',
                   'pepseq2': 'PepSeq2',
                   'xlink1': 'LinkPos1',
                   'xlink2': 'LinkPos2',
                   'score': 'Score',
                   'ID': 'Id'}

    xinet.rename(index=str,
                 columns=rename_dict,
                 inplace=True)
    
    xinet.to_csv(outpath + '.csv',
                 index=False)   

def WriteDynamXL(xtable, outpath):
    """
    Converts xTable data into cross-link information
    file for the dynamXL analysis software:
    
    K189 NZ K192 NZ 18
    K189 NZ K196 NZ 5
    K200 NZ K192 NZ 1
    K236 NZ K189 NZ 3
    
    dynamxl.chem.ox.ac.uk
    
    :params: xtable data table structure
    :params: outpath_ path to write file
    """
    def aa_and_pos(row):
        """
        Calculate tne amino acid and the absolute position
        of a crosslink based on the sequence of an xl-peptide,
        the relative and the absolute position of the cross-link
        """        
        aa1 = row['pepseq1'][int(row['xlink1'])-1]
        id1 = str(aa1) + str(row['xpos1'])
        
        aa2 = row['pepseq2'][int(row['xlink2'])-1]
        id2 = str(aa2) + str(row['xpos2'])
        
        return id1, id2
        
    def xlink_atom_from_AA(AA):
        """
        Return the typical cross-linked atom in PDB code
        for a specific amino-acid
        """
        
        if AA == 'K':
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
    
    dynamxl.to_csv(outpath + '.txt',
                   sep = '\t',
                   header=False,
                   float_format='%.3f',
                   index=False) 