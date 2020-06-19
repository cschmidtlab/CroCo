# -*- coding: utf-8 -*-

"""
Functions to read and write xTable data.

"""

import pandas as pd

if __name__ == '__main__' or __name__ =='xTable':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def _join_list_by_semicolon(entry):
    if isinstance(entry, list):
        return ';'.join([str(x) for x in entry])
    else:
        return entry

def _retain_topn(xtable, group, scoring, n, direction):
    """
    Return an xTable that contains only the n highest/lowest entries in the
    scoring column

    Args:
        group(str): Column name to group by (only topN PSMs per group will be returned)
        scoring(str): Column name to score (scoring will define the order in the groups)
        n(int): Number of rows retained
        direction(str): 'lowest' or 'highest'. Return the lowest or highest scoring rows

    Returns:
        xtable(pd.dataframe)
    """

    print('entering retain_topn')

    if n != None:
        try:
            n = int(n)
        except:
            raise Exception('[xTable Write] Please provide integer of rows to retain. Use 0 for all')

    if isinstance(group, str):
        group_list = [x.strip() for x in group.split(',')]
    elif isinstance(group, list):
        group_list = group
    else:
        raise Exception('[xTable Write] Group by must be a comma separated string or a list')
    for g in group_list:
        if g not in xtable.columns:
            raise Exception('[xTable Write] groupby string not found in column names')

    if scoring not in xtable.columns:
        raise Exception('[xTable Write] scoring string not found in column names')
    if direction not in ['lowest', 'highest']:
        raise Exception('[xTable Write] Direction string must be "lowest" or "highest"')

    if (n is None) or (n == 0):
        pass
    elif direction == 'lowest':
        print('lowest')
        xtable = xtable.sort_values(scoring, axis=0).groupby(group_list).head(n)
    elif direction == 'highest':
        print('highest')
        xtable = xtable.sort_values(scoring, axis=0).groupby(group_list).tail(n)
    else:
        raise Exception('[xTable Write] Something went wrong during filtering.')

    return xtable


def Write(xtable, outpath, do_filter=False, group='ID, rawfile', scoring='score', n=None, direction='lowest'):
    """
    writes an xtable data structure to file (in csv format)

    Args:
        xtable: data table structure
        outpath to write file (w/o file extension!)
        do_filter: Whether to filter the xTable or not
        group(str): Column name to group by (only topN PSMs per group will be returned)
        scoring(str): Column name to score (scoring will define the order in the groups)
        n(int): Number of rows retained if filtering is active
        direction(str): 'lowest' or 'highest'. Return the lowest or highest scoring rows
    """
    
    if do_filter:
        print('[xTable Write] Size before filtering: {}'.format( xtable.size))
        xtable = _retain_topn(xtable, group, scoring, n, direction)
        print('[xTable Write] Size after filtering: {}'.format( xtable.size))
    
    # select only object dtypes as lists will anyways be found only in those
    # and applymap struggles with nullable int64 dtype
    xtable.loc[:,xtable.dtypes == 'object'] = xtable.loc[:,xtable.dtypes == 'object'].applymap(_join_list_by_semicolon)

    xtable.to_csv(hf.compatible_path(outpath) + '.csv',
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
            s = pd.read_csv(hf.compatible_path(file))
            allData.append(s)
        except:
            raise Exception('[xTable Read] Failed opening file: {}'.format(file))

    xtable = pd.concat(allData, sort=False)
    # Remove rows that contain no values (may be caused by Excel saving routine for csv files)

    xtable.dropna(axis=0, how='all', inplace=True)
    # convert only those columns to lists where lists are expected
    xtable[['modmass1','modmass2']] = xtable[['modmass1', 'modmass2']]\
        .applymap(lambda x: hf.convert_to_list_of(x, float))

    xtable[['modpos1', 'modpos2']] = xtable[['modpos1' ,'modpos2']]\
        .applymap(lambda x: hf.convert_to_list_of(x, int))

    xtable[['mod1', 'mod2']] = xtable[['mod1', 'mod2']]\
        .applymap(lambda x: hf.convert_to_list_of(x, str))

    xtable = hf.order_columns(xtable, col_order, compact)

    xtable = xtable.apply(pd.to_numeric, errors = 'ignore')

    return xtable

if __name__ == '__main__':
    xtable = Read(r'C:\Users\User\Documents\03_software\python\CroCo\testdata\ExampleData\output\all_merged_xTable.csv')

    Write(xtable, r'C:\Users\User\Documents\03_software\python\CroCo\testdata\ExampleData\output\all_merged_xTable2', do_filter=True, n=2)
