# -*- coding: utf-8 -*-

"""
Functions to write pLabel data.

This script is part of the CroCo cross-link converter project
"""

import pandas as pd
import os
import numpy as np
import re

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def Write(xtable, outpath, mgfDir, xlinker, mergepLabel = False):
    """
    Converts xtable data structure to (multiple) input file(s)
    for the pLabel cross-link annotation tool

    Args:
        xtable: data table structure
        outpath: path to write file (w/o file extension!)
        xlinker: xlinker as given to pLabel
        mergepLabel (bool): Whether to generate a new MGF and single pLabel file

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
        elif xtype == 'homomultimeric':
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
            # add position for: cterm1, xlink, nterm2
            modpos2 = [x + (len(pepseq1) + 3) for x in modpos2]

            if not hf.isNaN(modpos1):
                mods.extend(hf.toList(mod1))
                modposs.extend(hf.toList(modpos1))

            mods.extend(hf.toList(mod2))
            modposs.extend(hf.toList(modpos2))

        elif not hf.isNaN(modpos1):

            mods.extend(hf.toList(mod1))
            modposs.extend(hf.toList(modpos1))

        for mod, pos in zip(mods, modposs):
            modlabels.append('{},{}'.format(int(pos), mods2num[mod]))

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

    def ParseMGF(filenames, mgfDir):
        """
        Parse all mgf files matching to a list of filenames and extract all
        TITLE arguments as list

        Args:
            filenames: List of rawfile basenames to look for
            mgfDir: Path to mgf files
        Returns:
            titles2mgfoffset: Dictionary mapping mgf header titles to
                                    the resp rawfile and the position within
                                    the mgf-file
        """

        titles2mgfoffset = {}
        localMGFFiles = []

        # collect mgf file names in mgfDir
        for file in os.listdir(mgfDir):
            if file.endswith('.mgf'):
                localMGFFiles.append(file)

        # check which files referenced in xtable are present in the dir
        mgfToOpen = []
        mgfNotFound = []
        for file in filenames:
            if file + '.mgf' in localMGFFiles:
                mgfToOpen.append(file + '.mgf')
            else:
                mgfNotFound.append(file + '.mgf')

        # raise error if file is missing
        if len(mgfNotFound) > 0:
            raise Exception('The following mgf files were not found at the ' +
                            'specified directory: {}'.format(', '.join(mgfNotFound)))

        pattern = re.compile(r'TITLE=([^\.]+\.\d+\.\d+\.\d+\.\d+.*$)')

        # parse the mgf files for titles
        for f in mgfToOpen:
            mgfFile = os.path.join(mgfDir, f)
            with open(mgfFile) as inf:
                offset_last = 0
                offset_before_last = 0
                for line in inf.readlines():
                    if line.startswith('TITLE='):
                        # in case of pXtract:
                        # TITLE=2017_08_04_SVs_BS3_16.2419.2419.4.dta
                        # for MSConvert with TPP compatibility:
                        # TITLE=2017_08_18_SK_3.1093.1093.2 File:"2017_08_18_SK_3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=1093"
                        # MSConvert w/o TPP:
                        # TITLE=2017_08_18_SK_3.1093.1093.2
                        if pattern.match(line):
                            m = pattern.match(line)
                            title = m.group(1)
                        else:
                            raise(Exception('Title not found'))

                        titles2mgfoffset[title.upper()] = mgfFile, offset_before_last
                    offset_before_last = offset_last
                    offset_last += len(line) + 1

        return titles2mgfoffset

    rawfiles = xtable['rawfile'].unique().tolist()

    titles2mgfoffset = ParseMGF(rawfiles, mgfDir)
    allTitles = list(set(titles2mgfoffset.keys()))

    if not mergepLabel:

        # separate by rawfile
        for rf in rawfiles:
            xtablePerRawfile = xtable[xtable['rawfile'] == rf].copy()
            # separate per type within rawfile
            outfile = os.path.join(outpath + '_' + rf + '.pLabel')
            print('Opening {} to write'.format(outfile))
            with open(outfile, 'w') as out:
                out.write('[FilePath]\n')
                out.write('File_Path=' + os.path.join(mgfDir, rf + '.mgf\n'))

                modifications = uniqueMods(xtablePerRawfile['mod2'].tolist() +\
                                           xtablePerRawfile['mod1'].tolist())

                mods2num = {} # dict mapping mod names to indices
                out.write('[Modification]\n')
                for idx, mod in enumerate(modifications):
                    out.write('{}={}\n'.format(idx+1, mod))
                    mods2num[mod] = idx+1

                out.write('[xlink]\n')
                out.write('xlink={}\n'.format(xlinker))

                out.write('[Total]\n')
                out.write('total={}\n'.format(len(xtablePerRawfile.index)))

                idx = 1
                for _, row in xtablePerRawfile.iterrows():

                    out.write('[Spectrum{}]\n'.format(idx))
                    idx += 1

                    scanno = str(int(row['scanno']))
                    prec_ch = str(int(row['prec_ch']))
                    
                    title = ''
                    for t in allTitles:
                        # add the scanno twice to the search string to avoid
                        # matching of substrings e.g. 2516 to 25164
                        if '.'.join([rf, scanno, scanno, prec_ch]).upper() in t:
                            title = t
                    
                    # Generate the spectrum title as used by pLabel from
                    # rawfile name, scanno and precursor charge
                    out.write('name={}\n'.format(title))

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

    elif mergepLabel:

        # Write only one pLabel file
        outfile = outpath + '.pLabel'
        outMGF = outpath + '.mgf'

        filesWithOffsetToCopy = []

        print('Opening {} to write'.format(outfile))
        with open(outfile, 'w') as plabel:

            plabel.write('[FilePath]\n')
            plabel.write('File_Path=' + outMGF + '\n')

            modifications = uniqueMods(xtable['mod2'].tolist() +\
                                       xtable['mod1'].tolist())

            mods2num = {} # dict mapping mod names to indices
            plabel.write('[Modification]\n')
            for idx, mod in enumerate(modifications):
                plabel.write('{}={}\n'.format(idx+1, mod))
                mods2num[mod] = idx+1

            plabel.write('[xlink]\n')
            plabel.write('xlink={}\n'.format(xlinker))

            plabel.write('[Total]\n')
            plabel.write('total={}\n'.format(len(xtable.index)))

            idx = 1
            for rf in rawfiles:

                xtablePerRawfile = xtable[xtable['rawfile'] == rf].copy()

                for _, row in xtablePerRawfile.iterrows():

                    toWrite = ''

                    toWrite += ('[Spectrum{}]\n'.format(idx))
                    idx += 1

                    scanno = str(int(row['scanno']))
                    prec_ch = str(int(row['prec_ch']))

                    title = ''
                    for t in allTitles:
                        # add the scanno twice to the search string to avoid
                        # matching of substrings e.g. 2516 to 25164
                        if '.'.join([rf, scanno, scanno, prec_ch]).upper() in t:
                            title = t

                    # Generate the spectrum title as used by pLabel from
                    # rawfile name, scanno and precursor charge
                    toWrite += ('name={}\n'.format(title))

                    filesWithOffsetToCopy.append(titles2mgfoffset[title])

                    toWrite += ('pep1={}\n'.format(GeneratepLabelPepString(row['type'],
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
                    plabel.write(toWrite)


        # Generate merged MGF file containing only the matching spectra
        print('Opening {} to write'.format(outMGF))
        with open(outMGF, 'w') as mgf:
            templates = set([file for file, offset in filesWithOffsetToCopy])
            for template in templates:
                with open(template, 'r') as t:
                    print('Opening {} to read'.format(template))

                    offsets = []
                    for file, offset in filesWithOffsetToCopy:
                        if file == template:
                            offsets.append(offset)

                    for o in offsets:
                        t.seek(o, 0)
                        while True:
                            # loop through the lines of the spectrum until end-signa
                            line = t.readline()
                            if line.startswith('END IONS'):
                                # leave loop if the current spectrum ends
                                mgf.write(line)
                                break
                            else:
                                mgf.write(line)

if __name__ == '__main__':
    import sys
    sys.path.append(r'C:\Users\User\Documents\03_software\python\CroCo\src')

    import croco

    infile = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\pLabel\test_xTable.xlsx'
    xTable = croco.xTable.Read(infile)

    outpath = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\pLabel\test_xTable'
    mgfDir = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\pLabel'

    Write(xTable, outpath, mgfDir, 'BS3', mergepLabel = True)