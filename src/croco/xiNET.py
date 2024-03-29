# -*- coding: utf-8 -*-

"""
Functions to write data for the xiNET data
    visualisation tool (http://crosslinkviewer.org)
"""

import pandas as pd
import numpy as np

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath):
    """
    Convert xtable data structure to xiNET
    data file 

    Args:
        xtable: data table structure
        outpath: path to write file
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

    # remove mono-links
    xinet = xinet[xinet['xlink2'].notnull()]

    # sort by score before dropping duplicates
    xinet.sort_values(by='score',
                      inplace=True,
                      ascending=False)
    # drop duplicates
    xinet.drop_duplicates(inplace=True,
                          keep='first',
                          subset=['prot1','prot2', 'pos1', 'pos2', 'xlink1', 'xlink2'])

    xinet['score'] = -np.log10(xinet['score'])

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

    if outpath.endswith('.csv'):
        xinet.to_csv(hf.compatible_path(outpath), index=False)
    else:
        xinet.to_csv(hf.compatible_path(outpath) + '.csv',
                    index=False)