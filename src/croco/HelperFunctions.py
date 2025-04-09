# -*- coding: utf-8 -*-
"""
HelperFunctions: Function that are used in multiple modules of CroCo
"""

import numpy as np
import os

### variables of repeated use that are centrally stored

regexDict = {'mgfTITLE': r'(.+?)\.\d+\.(\d+)\.(\d+)\.*\d*'}
default_col_order = ['rawfile', 'scanno', 'prec_ch',
                     'pepseq1', 'xlink1',
                     'pepseq2', 'xlink2',
                     'modmass1', 'modpos1', 'mod1',
                     'modmass2', 'modpos2', 'mod2',
                     'prot1', 'xpos1', 'prot2',
                     'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']


### Functions that are repeatedly used
def compatible_path(raw_path, encoding=None):
    """
    Convert paths on Win to overcome the win32 pathlength limit
    """
    if (not isinstance(raw_path, str) and
            encoding is not None):
        raw_path = raw_path.decode(encoding)
    path = os.path.abspath(raw_path)
    if (os.name == 'nt') and (len(path) > 255):
        print('[compatible_path] converting to Windows extended path')
        if path.startswith(u"\\\\"):
            return u"\\\\?\\UNC\\" + path[2:]
        return u"\\\\?\\" + path
    else:
        return path


def categorize_inter_peptides(prot1, pos1, pepseq1, prot2, pos2, pepseq2):
    """
    Categorizes cross-linked peptides into inter, intra, homomultimeric and 
    sequential peptides
    """
    pepend1 = int(pos1) + len(pepseq1) - 1
    pepend2 = int(pos2) + len(pepseq2) - 1

    if prot1 != prot2:
        return 'inter'

    else:
        if (pos2 <= pepend1) and (pepend1 <= pepend2):
            return 'homomultimeric'
        elif (pos1 <= pepend2) and (pepend2 <= pepend1):
            return 'homomultimeric'
        elif (pepend1 + 1 == pos2) or (pepend2 + 1 == pos1):
            return 'sequential'
        else:
            return 'intra'


def order_columns(xtable, col_order, compact):
    """
    Sort columns of xtable by col_order and return the whole xtable including
    columns mentioned in col_order if keep is true. Otherwise, only return a
    minimal xTable
    """
    compact = bool(compact)
    if col_order is not None:
        if compact is False:
            # reorder columns to start with the xtable columns
            all_cols = list(xtable.columns.values)
            remaining_cols = [x for x in all_cols if x not in col_order]
            new_order = col_order + sorted(remaining_cols)

            xtable = xtable[new_order]
        elif compact is True:
            try:
                xtable = xtable[col_order]
            except Exception:
                raise Exception(
                    "[order_columns] Couldn't apply col order. Did you pass compact=True and a list of column-titles?")
        else:
            raise Exception('Compact argument passed to order_columns must be either True or False')

    return xtable


def generate_id(xl_type, prot1, xpos1, prot2, xpos2):
    """
    Return a link ID based on the type of the xlink
    """

    if xl_type in ['mono', 'loop']:
        xpos1 = int(xpos1)
        return '-'.join([str(prot1), str(xpos1)])
    elif xl_type in ['inter', 'intra', 'homomultimeric']:
        xpos1 = int(xpos1)
        xpos2 = int(xpos2)
        if xpos1 < xpos2:
            return '-'.join([str(prot1), str(xpos1), str(prot2), str(xpos2)])
        elif xpos1 == xpos2:
            prot_list = sorted([str(prot1), str(prot2)])
            return '-'.join([prot_list[0], str(xpos1), prot_list[1], str(xpos2)])
        else:
            return '-'.join([str(prot2), str(xpos2), str(prot1), str(xpos1)])
    else:
        return np.nan


def isnan(num):
    return num != num


def convert_to_list_of(input_str, typefunc, delimiter=';'):
    """
    Take an object that is not NaN, check if it contains a delimiter, split
    by delimiter and return list of elements of type typefunc
    
    Args:
        input_str: input object
        typefunc: e.g. Python int, str, or float
        delimiter (optional): string to split on
    Returns:
        List of objects of type typefunc
    """
    if not isnan(input_str):
        if not isinstance(input_str, str):
            return [input_str]
        else:
            input_list = input_str.split(delimiter)
            return [typefunc(x) for x in input_list]
    else:
        return input_str


def split_concatenated_lists(dataframe, where, delimiter=';'):
    """
    Splits each row of a dataframe that contains a delimiter-separated
    string into two columns with each element of the string in each row.
    
    Args:
        dataframe: dataframe to operate on
        where (list): column-name in which to find the strings
        delimiter: (optional) the delimiter-character to look for
    
    Returns:
        dataframe: Modified dataframe
    """

    if not isinstance(where, list):
        raise Exception('Please specify a list as where')

    dataframe['split_entry'] = False

    # recalculate the df-splitting for every column that should be split
    for w in where:

        rows = []  # store split rows
        rows2drop = []  # store rows to remove

        # iterate over all rows in the input df
        for idx, row in dataframe.iterrows():
            # skip rows containing NaN
            if isnan(row[w]):
                continue
            # if the delimiter is found at the specified column
            if delimiter in row[w]:
                # strip the delimiter from the end of the entry
                # split the string in the column
                elements = row[w].strip(delimiter).split(delimiter)
                # check if this was a split row
                if len(elements) > 1:
                    # append the original row to delete later
                    rows2drop.append(idx)
                    for element in elements:
                        # create working copy of row
                        mod_row = row.copy()
                        # replace the original string with one of its constituents
                        mod_row[w] = element + delimiter
                        # add an identifier for split entries
                        mod_row['split_entry'] = True
                        # append the new row at the bottom of the df
                        rows.append(mod_row)

        # drop the original rows
        dataframe.drop(dataframe.index[rows2drop], inplace=True)
        # append the new rows
        for row in rows:
            dataframe = dataframe.append(row)
        # reset the index
        dataframe = dataframe.sort_index()
        # reset the index (recount from 0 to N)
        dataframe = dataframe.reset_index(drop=True)

    return dataframe

def assign_type(row):
    """
    Assign mono, loop, inter and intra link
    based on prot1, prot2, xlink1 and xlink2 entries

    Args:
        row (Series): a series or list containing prot1, prot2, xlink1, xlink2
    Returns:
        str or np.nan: type of cross-link (inter, intra, loop, mono)
    """
    prot1, prot2, xlink1, xlink2 = row

    prot1 = str(prot1)
    prot2 = str(prot2)
    xlink1 = str(xlink1)
    xlink2 = str(xlink2)

    if prot2 != 'nan' and prot1 == prot2:
        t = 'intra'
    elif prot2 != 'nan':
        t = 'inter'
    elif prot2 == 'nan' and xlink2 != 'nan':
        t = 'loop'
    elif prot1 != 'nan' and prot2 == 'nan' and xlink1 != 'nan':
        t = 'mono'
    else:
        t = np.nan
    return t