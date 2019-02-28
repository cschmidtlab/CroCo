# -*- coding: utf-8 -*-

"""
Functions to read Xi processed crosslink data filtered with xiFDR.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import os

if __name__ == '__main__':
    import HelperFunctions as hf
    import Xi as xi
else:
    from . import Xi as xi
    from . import HelperFunctions as hf

def Read(xi_file, xifdr_linksPSM_file=None, xifdr_linearPSM_file=None, modstring=None, col_order=None, compact=False):
    """
    Collects data from Xi spectrum search filtered by xiFDR and returns an xtable data array.

    Args:
        xi_file: path to percolated Kojak file
        xifdr_linksPSM_file: xlink results file from xiFDR (contains PSM_xiFDR)
        xifdr_linearPSM_file: linear peptide results file from xiFDR (contains PSM_xiFDR)
        modstring (str): a string (like 'bs3') that is part of the monolink modification name in Xi
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        keep (bool): Whether to keep the columns of the original dataframe or not

    Returns:
        xtable: xtable data table
    """

    def inferMonoLinks(mod, modpos):
        """
        Take a list of modifications with their positions and return those modpositions
        that match a given modification name for the monolinker

        Args:
            mod (str): a ;-separated list of modification names
            modpos (str): a also ;-separated list of the corresponding positions

        Returns:
            position (int or np.nan)
        """
        if hf.isNaN(mod):
            return np.nan

        mods = mod.split(';')
        modposns = modpos.split(';')

        for idx, m in enumerate(mods):
            if modstring in m:
                return int(modposns[idx])
            else:
                return np.nan

    ### Get the standard xi-file ###

    if isinstance(xi_file, list):
        if len(xi_file) > 1:
            raise Exception('[xiFDR Read] Sorry! Only one xi-file file per conversion is allowed to unambiguously relate it to a xi-file')
        xi_file = xi_file[0]

    print('[xiFDR Read] Reading xi-file: {}'.format(xi_file))

    xi_dtypes = {'Scan': pd.Int64Dtype(),
                 'PrecoursorCharge': pd.Int64Dtype(),
                 'BasePeptide1': str,
                 'ProteinLink1': pd.Int16Dtype(),
                 'BasePeptide2': str,
                 'ProteinLink2': pd.Int16Dtype(),
                 'Protein1': str,
                 'Protein2': str,
                 'Start1': pd.Int32Dtype(),
                 'Start2': pd.Int32Dtype(),
                 'Link1': pd.Int16Dtype(),
                 'Link2': pd.Int16Dtype(),
                 'match score': float
                 }

    try:
        xiraw = pd.read_csv(hf.FSCompatiblePath(xi_file), delimiter=',', dtype=xi_dtypes)
    except:
        raise Exception('[xTable Read] Failed opening file: {}'.format(xi_file))

    # list to collect the tables from the FDR processed files
    allFDR = list()

    ### get the xiFDR file containing the cross-links ###
    if xifdr_linksPSM_file != None:
        if not 'Links_PSM_xiFDR' in xifdr_linksPSM_file:
            raise Exception('[xiFDR Read] The string "Links_PSM_xiFDR" is missing in your input file. Did you choose the right file?')

        print('[xiFDR Read] Reading xiFDR-file: {}'.format(xifdr_linksPSM_file))

        try:
            linkXiFDR = pd.read_csv(hf.FSCompatiblePath(xifdr_linksPSM_file), delimiter=',')

            linkXiFDR.rename(columns={'run': 'Run',
                                      'scan': 'Scan',
                                      'Protein1': 'Protein1_FDR',
                                      'Protein2': 'Protein2_FDR'}, inplace=True)
        except Exception as e:
            raise Exception('[xiFDR Read] Error while reading the xiFDR file: {}'.format(e))

        # Merge with left join (only keys that are in the percolated DF will be re-
        # tained)
        s = pd.merge(linkXiFDR, xiraw, on=['Run', 'Scan'], how='left')

        allFDR.append(s)

    ### get the xiFDR file containing the monlinks ###

    if xifdr_linearPSM_file != None:

        if not 'Links_Linear_PSM_xiFDR' in xifdr_linearPSM_file:
            raise Exception('[xiFDR Read] The string "Links_Linear_PSM_xiFDR" is missing in your input file. Did you choose the right file?')

        print('[xiFDR Read] Reading xiFDR-file: {}'.format(xifdr_linearPSM_file))

        try:
            linaerXiFDR = pd.read_csv(hf.FSCompatiblePath(xifdr_linearPSM_file), delimiter=',')

            linaerXiFDR.rename(columns={'run': 'Run',
                                        'scan': 'Scan',
                                        'Protein1': 'Protein1_FDR',
                                        'Protein2': 'Protein2_FDR'}, inplace=True)
        except Exception as e:
            raise Exception('[xiFDR Read] Error while reading the xiFDR file: {}'.format(e))

        # Merge with left join (only keys that are in the percolated DF will be re-
        # tained)
        s = pd.merge(linaerXiFDR, xiraw, on=['Run', 'Scan'], how='left')

        allFDR.append(s)

    # only continue if at least one of the FDR-processed files has been provided
    if len(allFDR) == 0:
        raise Exception('[xiFDR] YOu must provide either a linear or a PSM FDR file!')
    else:
        data = pd.concat(allFDR)

    ### Process the data to comply to xTable format
    xtable = data.rename(columns={'Scan': 'scanno',
                                  'PrecoursorCharge': 'prec_ch',
                                  'BasePeptide1': 'pepseq1',
                                  'ProteinLink1': 'xpos1',
                                  'BasePeptide2': 'pepseq2',
                                  'ProteinLink2': 'xpos2',
                                  'ModificationMasses1': 'modmass1',
                                  'ModificationMasses2': 'modmass2',
                                  'Modifications1': 'mod1',
                                  'Modifications2': 'mod2',
                                  'Protein1': 'prot1',
                                  'Protein2': 'prot2',
                                  'Start1': 'pos1',
                                  'Start2': 'pos2',
                                  'Link1': 'xlink1',
                                  'Link2': 'xlink2',
                                  'ModificationPositions1': 'modpos1',
                                  'ModificationPositions2': 'modpos2',
                                  'match score': 'score'
                                  })

    xtable['rawfile'] = xtable['Source'].apply(xi.rawfile_from_source)

    print('[xiFDR Read] Rawfile from source')

    if modstring != None:
        xtable.loc[xtable['xlink1'].isnull(), 'xlink1'] =\
            np.vectorize(inferMonoLinks)(xtable.loc[xtable['xlink1'].isnull(), 'mod1'],
                                         xtable.loc[xtable['xlink1'].isnull(), 'modpos1'])

        xtable.loc[xtable['xlink1'].notnull(), 'xpos1'] =\
            xtable.loc[xtable['xlink1'].notnull(), 'pos1'].astype(int) +\
            xtable.loc[xtable['xlink1'].notnull(), 'xlink1'].astype(int)

        print('[xiFDR Read] Inferred monolinks from modifications')

    # assign cateogries of cross-links based on identification of prot1 and prot2
    xtable['type'] = xtable[['prot1', 'prot2', 'xlink1', 'xlink2']].apply(\
        xi.assign_type, axis=1)

    print('[xiFDR Read] assigned Type')

    # generate an ID for every crosslink position within the protein(s)
    xtable['ID'] =\
        pd.Series(np.vectorize(hf.generateID,
                               otypes=['object'])(xtable['type'],
                                                  xtable['prot1'],
                                                  xtable['xpos1'],
                                                  xtable['prot2'],
                                                  xtable['xpos2']),
                 index=xtable.index).replace('nan', np.nan)

    print('[xiFDR Read] generated ID')

    if len(xtable[xtable['type'] == 'inter']) > 0:
        # Reassign the type for inter xlink to inter/intra/homomultimeric
        onlyInter = xtable['type'] == 'inter'
        xtable.loc[onlyInter, 'type'] =\
            np.vectorize(hf.categorizeInterPeptides)(xtable[onlyInter]['prot1'],
                                                     xtable[onlyInter]['pos1'],
                                                     xtable[onlyInter]['pepseq1'],
                                                     xtable[onlyInter]['prot2'],
                                                     xtable[onlyInter]['pos2'],
                                                     xtable[onlyInter]['pepseq1'])
        print('[xQuest Read] categorized inter peptides')
    else:
        print('[xQuest Read] skipped inter peptide categorization')

    xtable['xtype'] = np.nan

    xtable['search_engine'] = 'XiSearch and XiFDR'

    xtable = hf.applyColOrder(xtable, col_order, compact)

    return xtable
if __name__ == '__main__':

    # defines the column headers required for xtable output
    col_order = [ 'rawfile', 'scanno', 'prec_ch',
                  'pepseq1', 'xlink1',
                  'pepseq2', 'xlink2', 'xtype',
                  'modmass1', 'modpos1', 'mod1',
                  'modmass2', 'modpos2', 'mod2',
                  'prot1', 'xpos1', 'prot2',
                  'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

    os.chdir(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\PK\Xi')
    xi_file = r'XI_results_XiVersion1.6.739.csv'
    xifdr_linksPSM_file = r'XiFDR_5_FDR_Links_PSM_xiFDR1.0.22.csv'
    xifdr_linearPSM_file = r'XiFDR_5_FDR_Links_Linear_PSM_xiFDR1.0.22.csv'
    xtable = Read(xi_file, xifdr_linearPSM_file=xifdr_linearPSM_file, xifdr_linksPSM_file=xifdr_linksPSM_file, modstring='bs3', compact=True)
