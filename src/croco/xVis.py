# -*- coding: utf-8 -*-

"""
Functions to write data structures as input for the xVis webserver (https://xvis.genzentrum.lmu.de/CrossVisNoLogin.php).

"""

import pandas as pd

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath):
    """
    Convert xtable data structure to cross-link
    data file for xVis data visualisation tool

    Args:
        xtable (pandas.DataFrame): data table structure
        outpath (str): path to write file
    """
    xvis = xtable.loc[:,['prot1','prot2', 'xpos1', 'xpos2', 'score']]

    # remove mono-links
    xvis = xvis[xvis['xpos2'].notnull()]

    # sort by score before dropping duplicates
    xvis.sort_values(by='score',
                     inplace=True,
                     ascending=False)
    # drop duplicates
    xvis.drop_duplicates(inplace=True,
                         keep='first',
                         subset=['prot1','prot2', 'xpos1', 'xpos2'])

    rename_dict = {'prot1':'Protein1',
                   'prot2':'Protein2',
                   'xpos1': 'AbsPos1',
                   'xpos2': 'AbsPos2',
                   'score': 'Id-Score'}
    xvis.rename(index=str,
                columns=rename_dict,
                inplace=True)

    if outpath.endswith('.csv'):
        xvis.to_csv(hf.compatible_path(outpath), index=False)
    else:
        xvis.to_csv(hf.compatible_path(outpath) + '.csv',
                index=False)