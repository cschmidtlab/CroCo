# -*- coding: utf-8 -*-

"""
Functions to write xTable data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd

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
    
    def convert_separated_string(input):
        if isinstance(input, str):

            string = str(input)
            if ';' in string:
                return string.split(';')
            else:
                return string
        else:
            return input

    xtable = pd.read_excel(inpath)
    xtable = xtable.applymap(convert_separated_string)
    
    return xtable
    
    
if __name__ == '__main__':
    xtable = Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')