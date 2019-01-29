# -*- coding: utf-8 -*-

"""
Functions to read and write xTable data.

"""

import pandas as pd

if __name__ == '__main__' or __name__ =='xTable':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath, col_order=None, compact=False):
    """
    writes an xtable data structure to file (in xlsx format)

    Args:
        xtable: data table structure
        outpath to write file (w/o file extension!)
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to compact the xTable to only those columns listed in col_order
    """

    def convert_list(entry):
        if isinstance(entry, list):
            return ';'.join([str(x) for x in entry])
        else:
            return entry

    xtable = xtable.applymap(convert_list)

    xtable = hf.applyColOrder(xtable, col_order, compact)

    xtable.to_excel(hf.FSCompatiblePath(outpath) + '.xlsx',
                    index=False)


def Read(xTable_files, col_order=None, compact=False):
    """
    Read an xTable data structure from file

    Args:
        xTable_files: path to the xtable file(s)
        col_order (list): List of xTable column titles that are used to sort and compress the resulting datatable
        compact (bool): Whether to compact the xTable to only those columns listed in col_order
    Returns:
        xtable: xTable dataframe object
    """

    # convert to list if the input is only a single path
    if not isinstance(xTable_files, list):
        xTable_files = [xTable_files]
    
    allData = list()
    
    for file in xTable_files:
        try:
            s = pd.read_excel(hf.FSCompatiblePath(file))
            allData.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))
    
    xtable = pd.concat(allData)
    # convert only those columns to lists where lists are expected
    xtable[['modmass1','modmass2']] = xtable[['modmass1', 'modmass2']]\
        .applymap(lambda x: hf.convertToListOf(x, float))

    xtable[['modpos1', 'modpos2']] = xtable[['modpos1' ,'modpos2']]\
        .applymap(lambda x: hf.convertToListOf(x, int))

    xtable[['mod1', 'mod2']] = xtable[['mod1', 'mod2']]\
        .applymap(lambda x: hf.convertToListOf(x, str))

    xtable = hf.applyColOrder(xtable, col_order, compact)

    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    return xtable

if __name__ == '__main__':
    xtable = Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')