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

import SpectrumCalculations as Calc
import SpectrumReader as Reader

#%% Start of IonsFromSequence
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
            print('No mass_type specified, assuming monoisotopic masses')

    all_ions = []
    all_desc = []

    for pep in Calc.NTermPeptides(sequence):
        if 'a' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = Calc.AIonMass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)
        if 'b' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = Calc.BIonMass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)

    for pep in Calc.NTermPeptides(sequence):
        if 'y' in ion_types:
            for c in range(charge[0], charge[-1]+1):
                mass, desc = Calc.YIonMass(pep, aa_dict, c, mods)
                all_ions.append(mass)
                all_desc.append(desc)

    for c in range(charge[0], charge[1]+1):
        mass, desc = Calc.ParentionMass(sequence, aa_dict, c, mods)
        all_ions.append(mass)
        all_desc.append(desc)
    return all_ions, all_desc

#%% Start of IonsFromXlinkSequence
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

    all_ions = []
    all_desc = []

    # collect sequences from modified and unmodified peptides (n-term of CID)
    nterm_unmod1, nterm_mod1, nterm_unmod2, nterm_mod2 =\
        Calc.NTermXPeptides(peptide1, xlink1, peptide2, xlink2)
    # C-terminal of CID break
    cterm_unmod1, cterm_mod1, cterm_unmod2, cterm_mod2 =\
        Calc.CTermXPeptides(peptide1, xlink1, peptide2, xlink2)

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
                    mass, desc = Calc.AIonMass(pep, aa_dict, c,
                                                mods, peptide_type)
                    if mass < max_mass:
                        all_ions.append(mass + modmass/c) # treat second peptide as constant adduct
                        desc[0] += moddesc
                        all_desc.append(desc)
                        
            if 'b' in ion_types:
                for c in range(charge[0], charge[-1]+1):
                    mass, desc = Calc.BIonMass(pep, aa_dict, c,
                                                mods, peptide_type)
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
                    mass, desc = Calc.YIonMass(pep, aa_dict, c,
                                                mods, peptide_type)
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
        mass, desc = Calc.XParentionMass(pep1_mass, pep2_mass,
                                          peptide1, xlink1,\
                                          peptide2, xlink2, c, mods)
        if mass < max_mass:
            all_ions.append(mass)
            all_desc.append(desc)
    return all_ions, all_desc


#%%

def AssignAndPlotPSM(spec_mz, spec_intens, theo_mz, theo_desc, ppm, ax=None):
    """
    Annotate a given spectrum with a given sequence of theoretical
    mz and descriptions
    
    :params: spec_mz: list of mz values for the centroided spectrum
    :params: spec_intens: list of corresponding intensities
    :params: theo_mz: list of theoretical peaks to be ound in the spectrum
    :params: the_desc: corresponding description-list
    :params: ppm: error allowed during assign in parts-per-million
    :params: axis (optional): axis to plot the figure on, default currrent axis
    
    :return: assignment_error: list of relative errors in ppm for the assignment
    """
    
    if ax is None:
        ax = plt.gca()
    for idx, m in enumerate(mz):
        # store axis-object in variable for return
        ax.plot([m, m], # x1, x2
                    [0, intens[idx]], # y1, y2
                    'k-', # style
                    lw=1 # line width
                    )
        for jdx, n in enumerate(ions):
            if n >= m * (1+ppm[0]/10**6):
                if n <= m * (1+ppm[1]/10**6):
                    # add the relative error of the assignments to the list
                    assignment_error.append((n-m)/n * 10**6)
                    ax.plot([n, n],
                               [0, intens[idx]],
                               'r-' if desc[jdx][4] == 'alpha' else 'b-',
                               lw=2)
                    ax.text(n,
                              intens[idx]+1000,
                              '${0}_{{{1}}}^{{{2}+}}$'.format(desc[jdx][1],
                                                              desc[jdx][2],
                                                              desc[jdx][3]),
                              {'ha': 'left', 'va': 'bottom'},
                              rotation=90,
                              color='r' if desc[jdx][4] == 'alpha' else 'b')
            
    ax.set_xlim(min(mz),max(mz))
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')

    return assignment_error

def PlotHist(data, limits, ax=None):
    """
    Plots a histogram with a fitted gaussion distribution over a set
    of data on a given axis
    
    :params: data: a list of values to plot a histogram from
    :params: limits: range of the plot as list [lower, upper]
    :params: ax: axis to plot the data on
    """
    if ax is None:
        ax = plt.gca()
    n, bins, patches = ax.hist(assignment_error,
                         50, # bins
                         normed=1,
                         alpha=0.8)
    mean = np.mean(assignment_error)
    variance = np.var(assignment_error)
    sigma = np.sqrt(variance)
    
    # add a 'best fit' line
    y = mlab.normpdf(bins, mean, sigma)
    ax.plot(bins, y, 'r--', linewidth=1)
    ax.axvline(mean, color='r', linewidth=1)

    ax.set_xlabel('Error in ppm')
    ax.set_ylabel('Probability Density')
    ax.set_xlim(limits)

#%%

fig, ax = plt.subplots(2)

ppm = [-30, 30]

assignment_error = []

ions, desc =  IonsFromXlinkSequence('YHPDKNPDNPEAADKFK', # peptid1
                                    15, # xlink1
                                    'KLALK', # peptide2
                                    1, # xlink2
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
    spectrum2offset = Reader.IndexMGF(f)
    
    mz, intens = Reader.ReadSpectrum(17079, # spectrum
                              f, # file handle
                              spectrum2offset) # spectrum dict

    assignment_error = AssignAndPlotPSM(mz, intens, ions, desc, ppm, ax=ax[0])
    
    PlotHist(assignment_error, ppm, ax[1])
    
    plt.show()