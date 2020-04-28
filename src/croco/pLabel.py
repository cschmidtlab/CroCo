# -*- coding: utf-8 -*-

"""
Functions to write pLabel data.
"""

import pandas as pd
import os
import numpy as np
import re

if __name__ == '__main__':
    import HelperFunctions as hf
else:
    from . import HelperFunctions as hf

def _unique_mods(modlist):
    """
    Go through a list containing lists of modifications, strings of
    modifications or NaN and extract all occuring unique mod-strings

    Args:
        modlist: list of modifications from xtable
    Returns:
        List of unique modifications
    """
    alist = []
    for element in modlist:
        # When multiple modifications occur on a peptide, the element is
        # a list
        if isinstance(element, list):
            for e in element:
                if e not in alist:
                    alist.append(str(e))

        # If only a single modification occurs, the element can be directly
        # added to the unique list if it is not NaN
        elif isinstance(element, float):
            if np.isnan(element):
                continue

        else:
            if element not in alist:
                alist.append(str(element))

    return list(set(alist))

def _generate_plabel_pepstring(xtype, xlink1, xlink2, pepseq1, pepseq2, score, mod1, mod2, modpos1, modpos2, mods2num):
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
        pepStringElements = [typeno, _cast_if_not_nan(xlink1, int),
                             _cast_if_not_nan(xlink2, int), pepseq1,
                             '{:.4f}'.format(score), pepseq2, '1']
    else:
        pepStringElements = [typeno, _cast_if_not_nan(xlink1, int),
                             pepseq1, '1']

    pepStringElements = [str(x) for x in pepStringElements]

    modlabels = []
    mods = []
    modposs = []

    if not hf.isnan(modpos2):

        # increment the modpos2 position to fit pLabel numbering
        # add position for: cterm1, xlink, nterm2
        modpos2 = [x + (len(pepseq1) + 3) for x in modpos2]

        if not hf.isnan(modpos1):
            mods.extend(_make_list(mod1))
            modposs.extend(_make_list(modpos1))

        mods.extend(_make_list(mod2))
        modposs.extend(_make_list(modpos2))

    elif not hf.isnan(modpos1):

        mods.extend(_make_list(mod1))
        modposs.extend(_make_list(modpos1))

    for mod, pos in zip(mods, modposs):
        modlabels.append('{},{}'.format(int(pos), mods2num[mod]))

    pepStringElements.extend(modlabels)

    # There is one space at the end of pLabel pepstrins
    return ' '.join(pepStringElements) + ' '

def _parse_mgf(filenames, mgfDir):
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
    for file in os.listdir(hf.compatible_path(mgfDir)):
        if file.endswith('.mgf'):
            localMGFFiles.append(file)

    # check which files referenced in xtable are present in the dir
    mgfToOpen = []
    mgfNotFound = []
    for file in filenames:
        if file + '.mgf' in localMGFFiles:
            mgfToOpen.append(file + '.mgf')
        # pXtract usually adds the fragmentation method after conversion
        # allow Orbitrap files to be recognised
        elif file + '_HCDFT' + '.mgf' in localMGFFiles:
            mgfToOpen.append(file + '_HCDFT' + '.mgf')
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
        with open(hf.compatible_path(mgfFile)) as inf:
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


def _make_list(strorList):
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

def _cast_if_not_nan(input, typefunc):
    if not hf.isnan(input):
        return typefunc(input)
    else:
        return input


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

    rawfiles = xtable['rawfile'].unique().tolist()

    if xlinker == '':
        raise Exception('Please provide a name for the cross-linker')

    titles2mgfoffset = _parse_mgf(rawfiles, mgfDir)
    allTitles = list(titles2mgfoffset.keys())

    if not mergepLabel:

        # separate by rawfile
        for rf in rawfiles:
            xtablePerRawfile = xtable[xtable['rawfile'] == rf].copy()
            # separate per type within rawfile
            outfile = os.path.join(outpath + '_' + rf + '.pLabel')
            print('Opening {} to write'.format(outfile))
            with open(hf.compatible_path(outfile), 'w') as out:
                out.write('[FilePath]\n')
                out.write('File_Path=' + os.path.join(mgfDir, rf + '.mgf\n'))

                modifications = _unique_mods(xtablePerRawfile['mod2'].tolist() +\
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
                    nothingFound = True

                    for idx, t in enumerate(allTitles):
                        # add the scanno twice to the search string to avoid
                        # matching of substrings e.g. 2516 to 25164
                        if '.'.join([rf, scanno, scanno, prec_ch]).upper() in t:
                            title = t
                            # set the variable to check if any matching title
                            # was found for a row
                            nothingFound = False
                            # remove the title from the list to avoid setting
                            # the same title twice
                            del allTitles[idx]
                            # leave the loop once the title has been removed
                            break

                    if nothingFound:
                        raise Exception('[pLabel writer] couldnt find a matching spectrum for {}. If converting an xTable that was not generated from pLink input searched with the same mgf-file, please activate the merge-pLabel-option'.\
                                            format('.'.join([rf, scanno, scanno, prec_ch]).upper()))

                    # Generate the spectrum title as used by pLabel from
                    # rawfile name, scanno and precursor charge
                    out.write('name={}.DTA\n'.format(title.upper()))

                    out.write('pep1={}\n'.format(_generate_plabel_pepstring(row['type'],
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
        # a list with new mgf spectrum titles to integrate non-pLink results
        # into the pLabel viewer
        new_titles_and_charges_for_copy = []

        print('[pLabel] Opening {} to write'.format(outfile))
        with open(hf.compatible_path(outfile), 'w') as plabel:

            plabel.write('[FilePath]\n')
            plabel.write('File_Path=' + outMGF + '\n')

            modifications = _unique_mods(xtable['mod2'].tolist() +\
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

            plabel_specno = 1
            for rf in rawfiles:

                xtablePerRawfile = xtable[xtable['rawfile'] == rf].copy()

                for _, row in xtablePerRawfile.iterrows():

                    toWrite = ''

                    toWrite += ('[Spectrum{}]\n'.format(plabel_specno))
                    plabel_specno += 1

                    scanno = str(int(row['scanno']))
                    prec_ch = str(int(row['prec_ch']))

                    nothingFound = True
                    title = ''
                    for idx, t in enumerate(allTitles):
                        # add the scanno twice to the search string to avoid
                        # matching of substrings e.g. 2516 to 25164
                        # the charge is not considered here as charge assignment
                        # can vary between different prorgammes
                        if '.'.join([rf, scanno, scanno]).upper() in t:
                            # generate a new mgf-spectrum title unique for this
                            # entry (pLabel cannot take a spectrum twice)
                            counter = 0
                            while '.'.join([rf, scanno, scanno, prec_ch, str(counter), '.dta']) in new_titles_and_charges_for_copy:
                                counter +=1
                            title = '.'.join([rf, scanno, scanno, prec_ch, str(counter)])
                            new_titles_and_charges_for_copy.append((title, prec_ch))

                            # save the position of each title in the MGF file
                            # for MGF-file merging
                            filesWithOffsetToCopy.append(titles2mgfoffset[t])
                            # set the variable to check if any matching title
                            # was found for a row
                            nothingFound = False
                            # leave the loop
                            break

                    if nothingFound:
                        raise Exception('[pLabel writer] couldnt find a matching spectrum for {}'.\
                                            format('.'.join([rf, scanno, scanno, prec_ch]).upper()))

                    # Generate the spectrum title as used by pLabel from
                    # rawfile name, scanno and precursor charge
                    toWrite += ('name={}.DTA\n'.format(title.upper()))

                    toWrite += ('pep1={}\n'.format(_generate_plabel_pepstring(row['type'],
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

        print('[pLabel] Merging MGF files')
        # Generate merged MGF file containing only the matching spectra
        print('Opening {} to write'.format(outMGF))
        with open(hf.compatible_path(outMGF), 'w') as mgf:
            # sequentially open all MGF-files to copy from
            templates = set([file for file, offset in filesWithOffsetToCopy])
            for template in templates:
                with open(hf.compatible_path(template), 'r') as t:
                    print('Opening {} to read'.format(template))

                    offsets = []
                    new_titles_and_charges = []
                    for idx, (file, offset) in enumerate(filesWithOffsetToCopy):
                        if file == template:
                            # parts to read from that file
                            offsets.append(offset)
                            # new titles to generate for each part read
                            new_titles_and_charges.append(new_titles_and_charges_for_copy[idx])

                    for idx, o in enumerate(offsets):
                        # move to the part of the file where the spectrum is stored
                        t.seek(o, 0)
                        new_title_and_charge = new_titles_and_charges[idx]
                        while True:
                            # loop through the lines of the spectrum until end-signa
                            line = t.readline()
                            if line.startswith('END IONS'):
                                # leave loop if the current spectrum ends
                                mgf.write(line)
                                break
                            elif line.startswith('TITLE'):
                                # change the title line
                                mgf.write('TITLE={}.DTA\n'.format(new_title_and_charge[0].upper()))
                            elif line.startswith('CHARGE'):
                                # change the charge line
                                mgf.write('CHARGE={}+\n'.format(new_title_and_charge[1]))
                            else:
                                mgf.write(line)

if __name__ == '__main__':
    import sys
    sys.path.append(r'C:\Users\User\Documents\03_software\python\CroCo\src')

    import croco

    infile = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\ExampleData\output\all_merged_xTable_intra.xlsx'
    xTable = croco.xTable.Read(infile)

    outpath = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\ExampleData\output\xTable_to_vis\pLabel'
    mgfDir = r'C:\Users\User\Documents\03_software\python\CroCo\testdata\ExampleData'

    Write(xTable, outpath, mgfDir, 'BS3', mergepLabel = True)