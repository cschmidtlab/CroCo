# -*- coding: utf-8 -*-

"""
Functions to write xTable data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd

def Write(xtable, outpath):
    """
    writes an xtable data structure to file (in xtable format)
    
    :params: xtable: data table structure
    :params: outpath to write file (w/o file extension!)
    """
    xtable.to_excel(outpath + '.xlsx',
                    index=False)