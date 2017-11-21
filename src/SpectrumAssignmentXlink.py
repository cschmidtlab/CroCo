# -*- coding: utf-8 -*-
"""
Annotate spectra with results from cross-link identification searches.
Part of the CroCo project.
Created on Fri Oct 13 12:59:06 2017
@author: aretaon

"""

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import csv
import numpy as np

def IonsFromSequence(sequence, ion_types, charge, mass_type, mods):
    """
    Calculates theoretical ions for an amino acid
    sequence with modifications.
    
    :params: sequence: Peptide sequence
    :params: ion_types: list of which ions to return (a/b/y)
    :params: charge state of the returned ions
    :params: mass_type: mono or average
    :params: mods: modifications to consider
    """
    
    # table contains amino acid masses - H2O for better calculation
    with open('../data/config/aa_masses.txt', mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        if mass_type == 'mono':
            aa_dict = {rows[0]:float(rows[1]) for rows in reader}
        elif mass_type == 'average':
            aa_dict = {rows[0]:float(rows[2]) for rows in reader}
        else:
            print('No ro wrong mass_type specified, assuming monoisotopic masses')

    def nterm_peptides(parent):
        length = len(parent)
        nterm = [parent[0:j+1] for j in range(length-1)]
        return nterm
  
    def cterm_peptides(parent):
        length = len(parent)
        cterm = [parent[i+1:length+1] for i in range(length-1)]
        return cterm

    def bion_mass(peptide, aa_dict, charge, mods):
        mass = 0
        for aa in peptide:
            mass += aa_dict[aa]
        # plus H20 from condensation - OH- loss from b-ion formation
        # proton charges = total charge - 1 charge from b-ion formation
        mass += 18.01528 - 17.00734 + (charge-1) * 1.00794
        mass /= charge # charging effect
        desc = (peptide, 'b_{}'.format(len(peptide)), charge, mods)
        return mass, desc

    def aion_mass(peptide, aa_dict, charge, mods):
        mass = 0
        for aa in peptide:
            mass += aa_dict[aa]
        mass += -28.0101 + charge * 1.00794 # M + nH+
        mass /= charge # charging effect
        desc = (peptide, 'a_{}'.format(len(peptide)), charge, mods)
        return mass, desc
    
    def yion_mass(peptide, aa_dict, charge, mods):
        mass = 0
        for aa in peptide:
            mass += aa_dict[aa]
        # plus H2O from condensation plus proton mass per charge
        mass += 18.01528 + charge * 1.00794 # M + nH+
        mass /= charge # charging effect
        desc = (peptide, 'y_{}'.format(len(peptide)), charge, mods)
        return mass, desc

    def parention_mass(peptide, aa_dict, charge, mods):
        mass = 0
        for aa in sequence:
            mass += aa_dict[aa]
        mass += 18.01528 + charge * 1.00794 # M + nH+
        mass /= charge
        desc = (peptide, '(M+{}H)'.format(charge), charge, mods)
        return mass, desc

    all_ions = []
    all_desc = []

    for pep in nterm_peptides(sequence):
        if 'a' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = aion_mass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)
        if 'b' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = bion_mass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)

    for pep in cterm_peptides(sequence):
        if 'y' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = yion_mass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)

    for c in range(charge[0], charge[1]+1):
        mass, desc = parention_mass(sequence, aa_dict, c, mods)
        all_ions.append(mass)
        all_desc.append(desc)
    return all_ions, all_desc


def IonsFromXlinkSequence(peptide1, xlink1, peptide2, xlink2, xlinker_mod,
                          ion_types, charge, mass_type, mods, max_mass):
    """
    Calculates theoretical ions for an amino acid
    sequence with modifications.
    
    :params: peptide1: Peptide1 sequence
    :params: xlink1: Relative cross-link position in peptide1
    :params: peptde2: Peptide2 sequence
    :params: xlink2: Relative cross-link possition peptide 2
    :params: xlinker_mod: Modification mass of the crosslinker
    :params: ion_types: list of which ions to return (a/b/y)
    :params: charge state of the returned ions
    :params: mass_type: mono or average
    :params: mods: modifications to consider
    :params: max_mass: maximum allowed mass (e.g. max m/z of quad)
    """
    
    # table contains amino acid masses - H2O for better calculation
    with open('../data/config/aa_masses.txt', mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        if mass_type == 'mono':
            aa_dict = {rows[0]:float(rows[1]) for rows in reader}
        elif mass_type == 'average':
            aa_dict = {rows[0]:float(rows[2]) for rows in reader}
        else:
            print('No ro wrong mass_type specified, assuming monoisotopic masses')

    def nterm_peptides(parent, start=1, end=None):
        """
        return all n-terminal peptides from start (including)
        to end (non-including)
        """
        if end == None:
            end = len(parent)

        nterm = [parent[0:j] for j in range(start, end)]
        return nterm
  
    def nterm_xpeptides(peptide1, xlink1, peptide2, xlink2):
        # truncation of peptide1
        nterm_unmod1 = nterm_peptides(peptide1, end=xlink1)
        nterm_mod1 = nterm_peptides(peptide1, start=xlink1)
        
        # truncation fo peptide2
        nterm_unmod2 = nterm_peptides(peptide2, end=xlink2)
        nterm_mod2 = nterm_peptides(peptide2, start=xlink2)
        
        return nterm_unmod1, nterm_mod1, nterm_unmod2, nterm_mod2
        
    def cterm_peptides(parent, start=1, end=None):
        if end == None:
            end = len(parent) 
        cterm = [parent[i:len(parent)] for i in range(start, end)]
        return cterm
    
    def cterm_xpeptides(peptide1, xlink1, peptide2, xlink2):
        # truncation of peptide 1
        cterm_unmod1 = cterm_peptides(peptide1, start=xlink1)
        cterm_mod1 = cterm_peptides(peptide1, end=xlink1)
        
        # truncation of peptide2
        cterm_unmod2 = cterm_peptides(peptide2, start=xlink2)
        cterm_mod2 = cterm_peptides(peptide2, end=xlink2)

        return cterm_unmod1, cterm_mod1, cterm_unmod2, cterm_mod2

    def bion_mass(peptide, aa_dict, charge, mods, peptide_type):
        ion_string = 'B' if peptide_type == 'alpha' else 'b'
        mass = 0
        for aa in peptide:
            mass += aa_dict[aa]
        # plus H20 from condensation - OH- loss from b-ion formation
        # proton charges = total charge - 1 charge from b-ion formation
        mass += 18.01528 - 17.00734 + (charge-1) * 1.00794
        mass /= charge # charging effect
        desc = [peptide, ion_string, len(peptide), charge,
                peptide_type, mods]
        return mass, desc

    def aion_mass(peptide, aa_dict, charge, mods, peptide_type):
        ion_string = 'A' if peptide_type == 'alpha' else 'a'
        mass = 0
        for aa in peptide:
            mass += aa_dict[aa]
        mass += -28.0101 + charge * 1.00794 # M + nH+
        mass /= charge # charging effect
        desc =  [peptide, ion_string, len(peptide), charge,
                 peptide_type, mods]
        return mass, desc
    
    def yion_mass(peptide, aa_dict, charge, mods, peptide_type):
        ion_string = 'Y' if peptide_type == 'alpha' else 'y'
        mass = 0
        for aa in peptide:
            mass += aa_dict[aa]
        # plus H2O from condensation plus proton mass per charge
        mass += 18.01528 + charge * 1.00794 # M + nH+
        mass /= charge # charging effect
        desc = [peptide, ion_string, len(peptide), charge,
                peptide_type, mods]
        return mass, desc

    def xparention_mass(pep1_mass, pep2_mass, peptide1, xlink1, peptide2,
                        xlink2, charge, mods):
        mass = sum((pep1_mass, pep2_mass))
        mass += charge * 1.00794 # M + nH+
        mass /= charge
        desc = ['{}-{}-{}-{}'.format(peptide1, xlink1, peptide2, xlink2),\
                '(M+{}H)'.format(charge), charge, mods]
        return mass, desc

    #%% Start of IonsFromXlinkSequence

    all_ions = []
    all_desc = []

    # collect sequences from modified and unmodified peptides (n-term of CID)
    nterm_unmod1, nterm_mod1, nterm_unmod2, nterm_mod2 =\
        nterm_xpeptides(peptide1, xlink1, peptide2, xlink2)
    # C-terminal of CID break
    cterm_unmod1, cterm_mod1, cterm_unmod2, cterm_mod2 =\
        cterm_xpeptides(peptide1, xlink1, peptide2, xlink2)

    print('Unmodified nterm peptides: {}'.format((', '.join(nterm_unmod1 +
                                                            nterm_unmod2))))

    print('Unmodified cterm peptides: {}'.format((', '.join(cterm_unmod1 +
                                                            cterm_unmod2))))

    print('Modified nterm peptides: {}'.format((', '.join(nterm_mod1 +
                                                            nterm_mod2))))

    print('Modified cterm peptides: {}'.format((', '.join(cterm_mod1 +
                                                            cterm_mod2))))


    ######################################################
    # Nterminal non xlinked peptides
    ######################################################
    def CalcABions(peptides, peptide_type, ion_types, charge, aa_dict,
                   moddesc, modmass, mods):
        for pep in peptides:
            if 'a' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = aion_mass(pep, aa_dict, c, mods, peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + modmass/c) # treat second peptide as constant adduct
                        desc[0] += moddesc
                        all_desc.append(desc)
                        
            if 'b' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = bion_mass(pep, aa_dict, c, mods, peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + modmass/c) # treat second peptide as constant adduct
                        desc[0] += moddesc
                        all_desc.append(desc)
                        
    CalcABions(nterm_unmod1, 'alpha', ion_types, charge, aa_dict,
               '', # moddesc
               0, # modmass
               mods)
    CalcABions(nterm_unmod2, 'beta', ion_types, charge, aa_dict,
               '', # moddesc
               0, # modmass
               mods)

    ######################################################
    # Nterminal xlinked peptides
    ######################################################
    
    # calculate mass of modification with peptide2 and xlinker
    pep2_mass = sum([aa_dict[aa] for aa in peptide2]) + 18.01528
    modmass1 = xlinker_mod + pep2_mass
    print('Mass of modification for peptide 1: {}'.format(modmass1))
    # same with peptide1
    pep1_mass = sum([aa_dict[aa] for aa in peptide1]) + 18.01528
    modmass2 = xlinker_mod + pep1_mass
    print('Mass of modification for peptide 2: {}'.format(modmass2))
    # define descriptions for peptide2 bound to (truncated) peptide1
    moddesc1 = '-{}-{}-{}'.format(xlink1, peptide2, xlink2)
    # same for peptide1 bound to (truncated) peptide2
    moddesc2 = '-{}-{}-{}'.format(xlink2, peptide1, xlink1)

    CalcABions(nterm_mod1, 'alpha', ion_types, charge, aa_dict,
               moddesc1, # moddesc
               modmass1, # modmass
               mods)
    
    CalcABions(nterm_mod2, 'beta', ion_types, charge, aa_dict,
               moddesc2, # moddesc
               modmass2, # modmass
               mods)

    ######################################################
    # Cterminal non-xlinked peptides
    ######################################################

    def CalcYions(peptides, peptide_type, ion_types, charge, aa_dict,
                   moddesc, modmass, mods):
        for pep in peptides:
            if 'y' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = yion_mass(pep, aa_dict, c, mods, peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + modmass/c)
                        desc[0] += moddesc
                        all_desc.append(desc)
    
    CalcYions(cterm_unmod1, 'alpha', ion_types, charge, aa_dict,
               '', # moddesc
               0, # modmass
               mods)    
        
    CalcYions(cterm_unmod2, 'beta', ion_types, charge, aa_dict,
               '', # moddesc
               0, # modmass
               mods) 

    ######################################################
    # Cterminal xlinked peptides
    ######################################################      

    CalcYions(cterm_mod1, 'alpha', ion_types, charge, aa_dict,
               moddesc1, # moddesc
               modmass1, # modmass
               mods)  

    CalcYions(cterm_mod2, 'beta', ion_types, charge, aa_dict,
               moddesc2, # moddesc
               modmass2, # modmass
               mods)  

    ######################################################
    # Parent ions
    ######################################################   

    for c in range(charge[0], charge[1]+1):
        mass, desc = xparention_mass(pep1_mass, pep2_mass, peptide1, xlink1,\
                                     peptide2, xlink2, c, mods)
        if mass < max_mass:
            all_ions.append(mass)
            all_desc.append(desc)
    return all_ions, all_desc


#%%

import re

def IndexMGF(mgf_file_handler):
    """
    Read mgf-file and return a dictionary of spectrum#: line pairs.
    Necessary for directly accessing spectra during run
    
    :params: mgf_file: file-handler of an opened mgf-file in tpp or std style
    
    :returns: spactrum2offset: Dict object with the scan numbers as keys and the
    corresponding offset of the mgf-file as values
    """
    
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
    return mz_list, intens_list

#%%
    
fig, ax = plt.subplots(2)

ppm = 20
assignment_error = []

ions, desc =  IonsFromXlinkSequence('ADALQAGASQFETSAAKLK', # peptid1
                                    17, # xlink1
                                    'VNVDKVLER', # peptide2
                                    5, # xlink2
                                    138.068, #mass of xlinker
                                    ['a', 'b', 'y'], # ion types
                                    [1,2], # min, max charge
                                    'mono', # mass type
                                    None, # modifications
                                    2000) # max mass

with open('peptides.log', 'w') as f:
    for idx, i in enumerate(ions):
        f.write('{}\t{}\n'.format(i, '\t'.join([str(i) for i in desc[idx]])))
                                    
with open('../testdata/SV_plink/2017_08_04_SVs_BS3_16.mgf', 'r') as f:
    spectrum2offset = IndexMGF(f)
    
    mz, intens = ReadSpectrum(43712, # spectrum
                              f, # file handle
                              spectrum2offset) # spectrum dict
         
    for idx, m in enumerate(mz):
        ax[0].plot([m, m], # x1, x2
                   [0, intens[idx]], # y1, y2
                   'k-', # style
                   lw=1 # line width
                   )
        for jdx, n in enumerate(ions):
            if n >= m * (1-ppm/10**6):
                if n <= m * (1+ppm/10**6):
                    # add the relative error of the assignments to the list
                    assignment_error.append((n-m)/n * 10**6)
                    ax[0].plot([n, n],
                               [0, intens[idx]],
                               'r-' if desc[jdx][4] == 'alpha' else 'b-',
                               lw=2)
                    ax[0].text(n,
                              intens[idx]+1000,
                              '${0}_{{{1}}}^{{{2}+}}$'.format(desc[jdx][1],
                                                              desc[jdx][2],
                                                              desc[jdx][3]),
                              {'ha': 'left', 'va': 'bottom'},
                              rotation=90)
            
    ax[0].set_xlim(min(mz),max(mz))
    ax[0].set_xlabel('m/z')
    ax[0].set_ylabel('Intensity')
    
   
    n, bins, patches = ax[1].hist(assignment_error,
                         50, # bins
                         normed=1,
                         alpha=0.8)
    mean = np.mean(assignment_error)
    variance = np.var(assignment_error)
    sigma = np.sqrt(variance)
    
    # add a 'best fit' line
    y = mlab.normpdf(bins, mean, sigma)
    ax[1].plot(bins, y, 'r--', linewidth=1)
    ax[1].axvline(mean, color='r', linewidth=1)

    ax[1].set_xlabel('Error in ppm')
    ax[1].set_ylabel('Probability Density')
    
    plt.show()