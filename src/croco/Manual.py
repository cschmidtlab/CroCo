# -*- coding: utf-8 -*-

"""
Functions to read Kojak data.

This script is part of the CroCo cross-link converter project
"""

import numpy as np
import pandas as pd

import re

def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

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

    :returns: list of rawfile, scanno, prec_chManual.init()

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