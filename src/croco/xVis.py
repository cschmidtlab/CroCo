# -*- coding: utf-8 -*-

"""
Functions to write xWalk data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath):
    """
    Converts xtable data structure to cross-link
    data file for xVis data
    visualisation tool
    (https://xvis.genzentrum.lmu.de/CrossVisNoLogin.php)

    :params: xtable: data table structure
    :params: outpath: path to write file
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
        xvis.to_csv(hf.FSCompatiblePath(outpath), index=False)
    else:
        xvis.to_csv(hf.FSCompatiblePath(outpath) + '.csv',
                index=False)