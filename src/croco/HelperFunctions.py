# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 14:18:32 2018

@author: User
"""

import pandas as pd
import numpy as np
import os

def FSCompatiblePath(raw_path, encoding=None):
    """
    Convert paths on Win to overcome the win32 pathlength limit
    """
    if encoding is not None:
        raw_path = raw_path.decode(encoding)
    path = os.path.abspath(raw_path)
    if os.name == 'nt':
        if path.startswith(u"\\\\"):
            return u"\\\\?\\UNC\\" + path[2:]
        return u"\\\\?\\" + path 
    else:
        return path

def categorizeInterPeptides(prot1, pos1, pepseq1, prot2, pos2, pepseq2):
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
            
def testVariable(prot1, pos1, pepseq1, prot2, pos2, pepseq2):
    
    return ', '.join([str(prot1), str(pos1), str(pepseq1), str(prot2), str(pos2), str(pepseq2)])

def applyColOrder(xtable, col_order, compact):
    """
    Sort columns of xtable by col_order and return the whole xtable including
    columns mentioned in col_order if keep is true. Otherwise only return a
    minimal xTable
    """
    compact = bool(compact)
    
    if compact is False:
        # reorder columns to start with the xtable columns
        all_cols = list(xtable.columns.values)
        remaining_cols = [x for x in all_cols if x not in col_order]
        new_order = col_order + sorted(remaining_cols)
    
        xtable = xtable[new_order]
    elif compact is True:
        xtable = xtable[col_order]
    else:
        raise Exception('Compact argument passed to applyColOrder must be either True or False')
    
    return xtable

def generateID(type, prot1, xpos1, prot2, xpos2):
    """
    Return a link ID based on the type of the xlink
    """

    if type in ['mono', 'loop']:
        return '-'.join([str(prot1), str(xpos1)])
    else:
        xpos1 = int(xpos1)
        xpos2 = int(xpos2)
        if xpos1 > xpos2:
            return '-'.join([str(prot1), str(xpos1), str(prot2), str(xpos2)])
        elif xpos1 == xpos2:
            prot_list = sorted([str(prot1), str(prot2)])
            return '-'.join([prot_list[0], str(xpos1), prot_list[1], str(xpos2)])
        else:
            return '-'.join([str(prot2), str(xpos2), str(prot1), str(xpos1)])

def toList(strorList):
    """
    take lists, floats or strings as input and return either the list or
    a one-element list of the string or float
    """
    if isinstance(strorList, list):
        return strorList

    elif isinstance(strorList, float):
        if not np.isnan(strorList):
            return [strorList]
    else:
        return [strorList]

def isNaN(num):
    return num != num

def convertToListOf(input, typefunc, delimiter=';'):
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
    if not isNaN(input):
        if not isinstance(input, str):
            return [input]
        else:
            inputList = input.split(delimiter)
            return [typefunc(x) for x in inputList]
    else:
        return input

def castIfNotNan(input, typefunc):
    if not isNaN(input):
        return typefunc(input)
    else:
        return input

def alphanum_string(s):
    """
    Method to clean strings from incorrect characters for file output
    """
    import re
    # new compiler that finds non-alphanumeric characters
    rex = re.compile(r'\W')
    # actually replace the strings
    result = rex.sub('_', s)

    # remove double occurences of _ in the string
    # initialize
    old_char = ''
    new_result = ''
    for char in result:
        if old_char != char:
            new_result += char
        else:
            # prevent removal of double occurences of other strings than _
            if old_char != '_':
                new_result += char
        old_char = char

    return new_result

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