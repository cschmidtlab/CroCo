# -*- coding: utf-8 -*-

"""
Functions to write crosslink data as defined by a user-provided template.

"""

import pandas as pd
if __name__ == '__main__':
    import sys
    sys.path.append(r'C:\Users\User\Documents\03_software\python\CroCo\src')
    import HelperFunctions as hf
    import croco
else:
    from . import HelperFunctions as hf

import re

def _from_match(match, row):
    toReplace = match.group(1)[1:-1]
    try:
        return str(row[toReplace])
    except Exception as e:
        raise Exception('Could not resolve string from xtable column header: {}'.format(e))


def Write(xtable, outpath, customTemplatePath):
    """
    writes an xtable data structure to file (in xtable format) based on a
    user-provided template file
    
    Args:
        xtable: data table structure
        outpath to write file (w/o file extension!)
        customTemplatePath: Path to template file defining the output structure
    
    """
    
    with open(hf.compatible_path(customTemplatePath), 'r') as tmpl:
        
        headerTemp = ''
        dataTemp = ''
        footerTemp = ''
        Templates = [headerTemp, dataTemp, footerTemp]
        pointer = None
        
        for line in tmpl.readlines():
            line = line.strip()
                        
            if line.startswith('[header]'):
                pointer = 0
            elif line.startswith('[data]'):
                pointer = 1
            elif line.startswith('[footer]'):
                pointer = 2
            elif pointer != None:
                Templates[pointer] += str(line) + '\n'
    
    substituteMatcher = re.compile(r'(\[.*?\])')
    
    with open(hf.compatible_path(outpath + '.csv'), 'w') as out:
        print('Writing to {}'.format(outpath + '.csv'))
        # write the header
        out.write(Templates[0])
        # write the data
        if Templates[1] != '':
            for idx, row in xtable.iterrows():
                out.write(substituteMatcher.sub(lambda match, row=row: _from_match(match, row), Templates[1]))
        # write footer
        out.write(Templates[2])
    
if __name__ == '__main__':
    xtable = croco.xTable.Read(r'C:\Users\User\Documents\02_experiments\05_croco_dataset\002_20180425\crosslink_search\pLink2_reports_xtable.xlsx')
    outpath = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\customTable_test'
    customTemplatePath = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\custom_table.txt'
    
    Write(xtable, outpath, customTemplatePath)