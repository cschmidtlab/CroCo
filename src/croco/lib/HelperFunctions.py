# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 14:18:32 2018

@author: User
"""

import pandas as pd
import numpy as np

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
                        mod_row[w] = element
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

if __name__ == '__main__':
        
    d = {'one' : ['ProtA', 'ProtB;ProtE;', 'ProtC', 'ProtD'],
         'two' : [1, 2, 3, 4],
         'three' : ['ProtA', 'ProtB;ProtE;', 'ProtC', 'ProtD']}
    
    my_df = pd.DataFrame(d)
    
    print(my_df)
    
    print(split_concatenated_lists(my_df, ['one', 'three']))