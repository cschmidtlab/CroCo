# -*- coding: utf-8 -*-

"""
Functions to write pLabel data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd
import os
import numpy as np

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath, mgfPath, xlinker):
    """
    Converts xtable data structure to (multiple) input file(s)
    for the pLabel cross-link annotation tool

    :params: xtable: data table structure
    :params: outpath: path to write file
    """

    def uniqueMods(modlist):
        """
        Go through a list containing lists of modifications, strings of
        modifications or NaN and extract all occuring unique mod-strings

        Args:
            modlist: list of modifications form xtable
        Returns:
            Lisr of unique modificaitons
        """
        alist = []
        for element in modlist:
            if isinstance(element, list):
                for e in element:
                    if e not in alist:
                        alist.append(e)

            elif isinstance(element, float):
                if not np.isnan(element):
                    if element not in alist:
                        alist.append(str(element))
            else:
                if element not in alist:
                    alist.append(element)

        return list(set(alist))

    def GeneratepLabelPepString(xtype, xlink1, xlink2, pepseq1, pepseq2, score, mod1, mod2, modpos1, modpos2, mods2num):
        """
        pep1=3 7 2 VFLLPDKK 22.7634 KKFETK 1 4,1
        1: 3-1 = cross, loop, mono
        2: xlink2
        3: xlink2
        4: pep1
        5: score?
        6: pep2
        7: unknown, mostly 1
        8: modification (e.g. 4,1 = modification 1 at position 4)

        """

        if xtype == 'inter':
            typeno = 3
        elif xtype == 'intra':
            typeno = 3
        elif xtype == 'loop':
            typeno = 2
        elif xtype == 'mono':
            typeno = 1

        if typeno > 1:
            pepStringElements = [typeno, hf.castIfNotNan(xlink1, int),
                                 hf.castIfNotNan(xlink2, int), pepseq1,
                                 '{:.4f}'.format(score), pepseq2, '1']
        else:
            pepStringElements = [typeno, hf.castIfNotNan(xlink1, int),
                                 pepseq1, '1']

        pepStringElements = [str(x) for x in pepStringElements]

        modlabels = []
        mods = []
        modposs = []

        if not hf.isNaN(modpos2):
            
            # increment the modpos2 position to fit pLabel numbering
            # add position for: nterm1, cterm1, xlink, nterm2
            modpos2 = [x + (len(pepseq1) + 4) for x in modpos2]

            if not hf.isNaN(modpos1):
                mods.extend(hf.toList(mod1))
                modposs.extend(hf.toList(modpos1))
                      
            mods.extend(hf.toList(mod2))             
            modposs.extend(hf.toList(modpos2))
            
        elif not hf.isNaN(modpos1):
            
            mods.extend(hf.toList(mod1))
            modposs.extend(hf.toList(modpos1))
            
        for mod, pos in zip(mods, modposs):
            modlabels.append('{},{}'.format(pos, mods2num[mod]))

        pepStringElements.extend(modlabels)

        return ' '.join(pepStringElements)

    def GeneratepLabelNameString(rawfile, scanno, prec_ch, rawfiles2titles):
        """
        2018_05_04_JB_MINCA5.16826.16826.3.0.DTA
        """
        
        rawfile = str(rawfile)
        scanno = str(int(scanno))
        prec_ch = str(int(prec_ch))
        
        allTitles = rawfiles2titles[rawfile]
    
        for title in allTitles:
            if title.startswith('.'.join([rawfile,
                                scanno,
                                scanno,
                                prec_ch])):
                return title.upper()


    def ExtractMGFTitles(filenames, mgfPath):
        """
        Parse all mgf files matching to a list of filenames and extract all
        TITLE arguments as list
        """
        
        rawfiles2titles = {}
        
        for f in filenames:
            if os.path.isfile(os.path.join(mgfPath, f + '.mgf')):
                rawfiles2titles[f] = []
                with open(os.path.join(mgfPath, f + '.mgf')) as inf:
                    for line in inf.readlines():
                        if line.startswith('TITLE='):
                            rawfiles2titles[f].append(line[6:].strip())
            else:
                raise Exception('MGF file not found: {}'.format(f + '.mgf'))

        return rawfiles2titles

    rawfiles = xtable['rawfile'].unique().tolist()

    rawfiles2titles = ExtractMGFTitles(rawfiles, mgfPath)

    # separate by rawfile
    for rf in rawfiles:
        xtablePerRawfile = xtable[xtable['rawfile'] == rf].copy()
        types = xtablePerRawfile['type'].unique().tolist()
        # separate per type within rawfile
        for t in types:
            xtablePerType = xtablePerRawfile[xtablePerRawfile['type'] == t].copy()
            outfile = os.path.join(outpath, '_'.join([rf, t + '.pLabel']))
            with open(outfile, 'w') as out:
                out.write('[FilePath]\n')
                out.write('File_Path=' + os.path.join(mgfPath, rf + '.mgf\n'))

                modifications = uniqueMods(xtablePerType['mod2'].tolist() +\
                                           xtablePerType['mod1'].tolist())

                mods2num = {} # dict mapping mod names to indices
                out.write('[Modification]\n')
                for idx, mod in enumerate(modifications):
                    out.write('{}={}\n'.format(idx+1, mod))
                    mods2num[mod] = idx+1

                out.write('[xlink]\n')
                out.write('xlink={}\n'.format(xlinker))

                out.write('[Total]\n')
                out.write('total={}\n'.format(len(xtablePerType.index)))

                idx = 1
                for _, row in xtablePerType.iterrows():
                    
                    out.write('[Spectrum{}]\n'.format(idx))
                    idx += 1
                    # name for the spectrum
                    out.write('name={}\n'.format(\
                                  GeneratepLabelNameString(row['rawfile'],
                                                           row['scanno'],
                                                           row['prec_ch'],
                                                           rawfiles2titles)))

                    out.write('pep1={}\n'.format(GeneratepLabelPepString(row['type'],
                                                                         row['xlink1'],
                                                                         row['xlink2'],
                                                                         row['pepseq1'],
                                                                         row['pepseq2'],
                                                                         row['score'],
                                                                         row['mod1'],
                                                                         row['mod2'],
                                                                         row['modpos1'],
                                                                         row['modpos2'],
                                                                         mods2num)))


if __name__ == '__main__':
    import sys
    sys.path.append(r'C:\Users\User\Documents\03_software\python\CroCo\src')
    
    import croco
    
    infile = r'C:\Users\User\Downloads\test\test_xTable.xlsx'
    xTable = croco.xTable.Read(infile)

    outpath = r'C:\Users\User\Downloads\test'
    mgfPath = r'C:\Users\User\Downloads\test'

    Write(xTable, outpath, mgfPath, 'BS3')