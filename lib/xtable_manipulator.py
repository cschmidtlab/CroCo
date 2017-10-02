"""
Manipulate xtable dataframes

"""

import pandas as pd

def Compare(xtable1, xtable2, xtable_names, compare_by):
    """
    Return a new dataframe where rows from two different dataframes are merged
    if identical (in compare_by values) and the respective entries for compare_what
    values are returned as distinct columns for plotting.
    
    :params: xtable1
    :params: xtable2
    :params: xtable_names
    :params: compare_by
    :params: compare_what
    
    :returns: xcomp_table
    """
    
    suffixes = [' ' + x for x in xtable_names]
   
    print('Comparing ' + ' & '.join(suffixes))  
    
    xcomp_table = pd.merge(xtable1,
                           xtable2,
                           on=compare_by,
                           how='inner',
                           suffixes=suffixes)
    
    return xcomp_table