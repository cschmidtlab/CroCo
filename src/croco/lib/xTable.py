# -*- coding: utf-8 -*-

"""
Functions to write xTable data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd
from . import HelperFunctions as hf

def Write(xtable, outpath):
    """
    writes an xtable data structure to file (in xtable format)
    
    Args:
        xtable: data table structure
        outpath to write file (w/o file extension!)
    
    """
    
    def convert_list(entry):
        if isinstance(entry, list):
            return ';'.join([str(x) for x in entry])
        else:
            return entry
    
    xtable = xtable.applymap(convert_list)
    
    xtable.to_excel(outpath + '.xlsx',
                    index=False)


def Read(inpath):
    """
    Read an xTable data structure form file
    
    Args:
        inpath: path to the xtable file
    
    Returns:
        xtable: xTable dataframe object
    """
    
    xtable = pd.read_excel(inpath)
    
    # convert only those columns to lists where lists are expected
    xtable[['modmass1','modmass2']] = xtable[['modmass1', 'modmass2']]\
        .applymap(lambda x: hf.convertToListOf(x, float))

    xtable[['modpos1', 'modpos2']] = xtable[['modpos1' ,'modpos2']]\
        .applymap(lambda x: hf.convertToListOf(x, int))
        
    xtable[['mod1', 'mod2']] = xtable[['mod1', 'mod2']]\
        .applymap(lambda x: hf.convertToListOf(x, str))

    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    return xtable
    
if __name__ == '__main__':
    xtable = Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')