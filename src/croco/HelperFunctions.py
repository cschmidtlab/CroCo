# -*- coding: utf-8 -*-
"""
HelperFunctions: Function that are used in multiple modules of CroCo
"""

import pandas as pd
import numpy as np
import os, re

### variables of repeated use that are centrally stored

regexDict = {'mgfTITLE': r'(.+?)\.\d+\.(\d+)\.(\d+)\.*\d*'}

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
    Categorizes cross-linked peptides into inter, intra, homomultimeric
    """
    pepend1 = int(pos1) + len(pepseq1) - 1
    pepend2 = int(pos2) + len(pepseq2) - 1
    
    if prot1 != prot2:
        return 'inter'
    
    else:
        if (pos2 <= pepend1 <= pepend2) or (pos1 <= pepend2 <= pepend1):
            return 'homomultimeric'
        else:
            return 'intra'
            
def order_columns(xtable, col_order, compact):
    """
    Sort columns of xtable by col_order and return the whole xtable including
    columns mentioned in col_order if keep is true. Otherwise only return a
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
            except Exception as e:
                raise Exception('[order_columns] Couldnt apply col order. Did you pass compact=True and a list of column-titles?')
        else:
            raise Exception('Compact argument passed to order_columns must be either True or False')
        
    return xtable

def generate_id(type, prot1, xpos1, prot2, xpos2):
    """
    Return a link ID based on the type of the xlink
    """

    if type in ['mono', 'loop']:
        xpos1 = int(xpos1)
        return '-'.join([str(prot1), str(xpos1)])
    elif type in ['inter', 'intra', 'homomultimeric']:
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

def convert_to_list_of(input, typefunc, delimiter=';'):
    """
    Take an object that is not NaN, check if it contains a delimiter, split
    by delimiter and return list of elements of type typefunc
    
    Args:
        input: input object
        typefunc: e.g. Python int, str, or float
        delimiter (optional): string to split on
    Returns:
        List of objects of type typefunc
    """
    if not isnan(input):
        if not isinstance(input, str):
            return [input]
        else:
            inputList = input.split(delimiter)
            return [typefunc(x) for x in inputList]
    else:
        return input

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

    # recalculate the df-splitting for every column that should be splitted
    for w in where:

        rows = []  # store splitted rows
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
                # check if this was a splitted row
                if len(elements) > 1:
                    # append the original row to delete later
                    rows2drop.append(idx)
                    for element in elements:
                        # create working copy of row
                        mod_row = row.copy()
                        # replace the original string with one of its constituents
                        mod_row[w] = element + delimiter
                        # add an identifier for splitted entries
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