# -*- coding: utf-8 -*-

"""
Functions to read and write xTable data.

"""

import pandas as pd

if __name__ == '__main__' or __name__ =='xTable':
    from HelperFunctions import convertToListOf, FSCompatiblePath
else:
    from .HelperFunctions import convertToListOf, FSCompatiblePath


def init(this_order):
    """
    Set required variables for conversion
    """
    global col_order
    col_order = this_order

def Write(xtable, outpath, keep=True):
    """
    writes an xtable data structure to file (in xlsx format)

    Args:
        xtable: data table structure
        outpath to write file (w/o file extension!)
        keep (optional): whether to restrict columns to the defined xTable cols
    """

    def convert_list(entry):
        if isinstance(entry, list):
            return ';'.join([str(x) for x in entry])
        else:
            return entry

    xtable = xtable.applymap(convert_list)

    if keep is True:
        # reorder columns to start with the xtable columns
        all_cols = list(xtable.columns.values)
        remaining_cols = [x for x in all_cols if x not in col_order]
        new_order = col_order + remaining_cols

        xtable = xtable[new_order]
    elif keep is False:
        xtable = xtable[col_order]

    xtable.to_excel(FSCompatiblePath(outpath) + '.xlsx',
                    index=False)


def Read(xTable_files):
    """
    Read an xTable data structure from file

    Args:
        xTable_files: path to the xtable file(s)

    Returns:
        xtable: xTable dataframe object
    """

    # convert to list if the input is only a single path
    if not isinstance(xTable_files, list):
        xTable_files = [xTable_files]
    
    allData = list()
    
    for file in xTable_files:
        try:
            s = pd.read_excel(FSCompatiblePath(file))
            allData.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))
    
    xtable = pd.concat(allData)
    # convert only those columns to lists where lists are expected
    xtable[['modmass1','modmass2']] = xtable[['modmass1', 'modmass2']]\
        .applymap(lambda x: convertToListOf(x, float))

    xtable[['modpos1', 'modpos2']] = xtable[['modpos1' ,'modpos2']]\
        .applymap(lambda x: convertToListOf(x, int))

    xtable[['mod1', 'mod2']] = xtable[['mod1', 'mod2']]\
        .applymap(lambda x: convertToListOf(x, str))

    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    return xtable

if __name__ == '__main__':
    xtable = Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')