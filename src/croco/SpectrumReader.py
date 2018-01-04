# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:37:37 2017

@author: User
"""
import re
import os
import pickle

def IndexMGF(mgf_file_handler):
    """
    Read mgf-file and return a dictionary of spectrum#: line pairs.
    Necessary for directly accessing spectra during run
    
    :params: mgf_file: file-handler of an opened mgf-file in tpp or std style
    
    :returns: spactrum2offset: Dict object with the scan numbers as keys and the
    corresponding offset of the mgf-file as values
    """
    
    # try restoring the index from a previously saved file
    indexpath = os.path.splitext(mgf_file_handler.name)[0] + '.idx'
    try:
        with open(indexpath, 'rb') as f:
            print('Reading existing indexfile at {}...'.format(os.path.abspath(indexpath)))
            return pickle.load(f)
    # if the indexfile was not present, generate a new index
    except:
        std_pattern = re.compile(r'TITLE=[^\.]+\.(\d+)\.\d+\.\d+')
        tpp_pattern = re.compile(r'TITLE=.*scan=(\d+)')
        
        print('Indexing...')
        
        spectrum2offset = {}
        offset = 0
        for line in mgf_file_handler:
            if line.startswith('TITLE='):
                # in case of pXtract:
                # TITLE=2017_08_04_SVs_BS3_16.2419.2419.4.dta
                # for MSConvert with TPP compatibility:
                # TITLE=2017_08_18_SK_3.1093.1093.2 File:"2017_08_18_SK_3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=1093"
                # MSConvert w/o TPP:
                # TITLE=2017_08_18_SK_3.1093.1093.2
                if std_pattern.match(line):
                    m = std_pattern.match(line)
                    spectrum = m.group(1)
                elif tpp_pattern.match(line):
                    m = tpp_pattern.match(line)
                    spectrum = m.group(1)
                else:
                    raise(Exception('Title not found'))
                spectrum2offset[int(spectrum)] = offset
            offset += len(line) + 1
    
        with open(indexpath, 'wb') as f:
            pickle.dump(spectrum2offset, f, pickle.HIGHEST_PROTOCOL)
    
        return spectrum2offset

def ReadSpectrum(scanno, mgf_file_handler, spectrum2offset):
    """
    Grab a spectrum by scan number from an indexed mgf-file
    
    :params: mgf_file_handler: handler for the input file
    :params: spectrum2offset: dictionary linking scan # and line #
    """
    mgf_file_handler.seek(spectrum2offset[scanno], 0)
    print('Starting at offset {}'.format(spectrum2offset[scanno]))
    mz_list, intens_list = [], []
    while True:
        # loop through the lines of the spectrum until end-signal
        line = mgf_file_handler.readline()
        if line[0] in '0123456789':
            # skip lines not starting with a number
            try:
                [mz, intens] = line.strip().split(' ')
            except:
                continue
            mz_list.append(float(mz))
            intens_list.append(float(intens))
        elif line.startswith('END IONS'):
            # leave loop if the current spectrum ends
            break
            
    mz2intens = dict(zip(mz_list, intens_list))
    
    return mz2intens
