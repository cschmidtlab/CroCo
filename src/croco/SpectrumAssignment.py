# -*- coding: utf-8 -*-
"""
Annotate spectra with results from cross-link identification searches.
Part of the CroCo project.
Created on Fri Oct 13 12:59:06 2017
@author: aretaon

"""
import os

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import csv
import numpy as np

if __name__ == '__main__':
    import SpectrumCalculations as Calc
    import SpectrumReader as Reader
    # set the start for relative paths if script run in standalone
    binPath = '../../bin'

else:
    import croco.SpectrumCalculations as Calc
    import croco.SpectrumReader as Reader
    # current path is bin-path if called from bin -folder
    binPath = '.'

#%% Start of IonsFromSequence
def IonsFromSequence(sequence, ion_types, charge, mass_type, mods):
    """
    Calculates theoretical ions for an amino acid
    sequence with modifications.

    :params: sequence: Peptide sequence
    :params: ion_types: list of which ions to return (a/b/y)
    :params: charge state of the returned ions
    :params: mass_type: mono or average
    :params: mods: modifications to consider (of type [[[mas1s_pep1, pos1_pep1], [mass2_pep1, pos2_pep1]], [[mass_pep2, pos_pep2]]])
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
    :params: peptide2: Peptide2 sequence
    :params: xlink2: Relative cross-link possition peptide 2
    :params: xlinker_mod: Modification mass of the crosslinker
    :params: ion_types: list of which ions to return (a/b/y)
    :params: charge state of the returned ions
    :params: mass_type: monoisotopic or average
    :params: mods: variable modifications to consider
    :params: max_mass: maximum allowed mass (e.g. max m/z of quad)
    """

    #%% CalcABIons
    def CalcABions(peptides, peptide_type, ion_types, charge, aa_dict,
                   moddesc, modmass, mods):
        """
        Calculate theoretical masses and corresponding descriptions
        for set of peptides sequences.

        :params: peptides: List of peptide sequences
        :params: peptide_type: string "alpha" or "beta" (first or last peptide truncated respectively)
        :params: ion_types: list of ion types to consider (e.g. ['b', 'y'])
        :params: charge: list of [min_charge, max_charge]
        :params: aa_dict: dict linking one-character ID to mass
        :params: moddesc: description to add to the xlinked peptides
        :params: modmass: mass of the cross-linker
        :params: mods: [[[mass1_1, mass2_1], [pos1_1, pos2_1]], [[mass2_1,], [pos2_1]]]
        """
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


    #%% CalcYions
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

    #%% Continue IonsFromXlinkSequence

    # table contains amino acid masses - H2O for better calculation
    aaPath = os.path.join(binPath, '../data/config/aa_masses.txt')
    with open(aaPath, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        if mass_type == 'monoisotopic':
            aa_dict = {rows[0]:float(rows[1]) for rows in reader}
        elif mass_type == 'average':
            aa_dict = {rows[0]:float(rows[2]) for rows in reader}
        else:
            print('No or wrong mass_type specified, assuming monoisotopic masses')

    all_ions = []
    all_desc = []

    # collect sequences from modified and unmodified peptides (n-term of CID)
    nterm_unmod1, nterm_mod1, nterm_unmod2, nterm_mod2, nterm_positions=\
        Calc.XPeptides(peptide1, xlink1, peptide2, xlink2, CTerm=False)
    # C-terminal of CID break
    cterm_unmod1, cterm_mod1, cterm_unmod2, cterm_mod2, cterm_positions =\
        Calc.XPeptides(peptide1, xlink1, peptide2, xlink2, CTerm=True)

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

    nterm_unmod1 = Calc.NTermPeptides(peptide1, end=xlink1)
    nterm_unmod2 = Calc.NTermPeptides(peptide2, end=xlink2)

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
    xpepmass1 = xlinker_mod + pep2_mass
    print('Mass of modification for peptide 1: {}'.format(xpepmass1))
    # same with peptide1
    pep1_mass = sum([aa_dict[aa] for aa in peptide1]) + 18.01528
    xpepmass2 = xlinker_mod + pep1_mass
    print('Mass of modification for peptide 2: {}'.format(xpepmass2))
    # define descriptions for peptide2 bound to (truncated) peptide1
    moddesc1 = '-{}-{}-{}'.format(xlink1, peptide2, xlink2)
    # same for peptide1 bound to (truncated) peptide2
    moddesc2 = '-{}-{}-{}'.format(xlink2, peptide1, xlink1)

    CalcABions(nterm_mod1, 'alpha', ion_types, charge, aa_dict,
               moddesc1, # moddesc
               xpepmass1, # modmass
               mods)

    CalcABions(nterm_mod2, 'beta', ion_types, charge, aa_dict,
               moddesc2, # moddesc
               xpepmass2, # modmass
               mods)

    ######################################################
    # Cterminal non-xlinked peptides
    ######################################################

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
               xpepmass1, # modmass
               mods)

    CalcYions(cterm_mod2, 'beta', ion_types, charge, aa_dict,
               moddesc2, # moddesc
               xpepmass2, # modmass
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

    ions2desc = dict(zip(all_ions, all_desc))

    return ions2desc


#%%

def AssignAndPlotPSM(mz2intens, ions2desc, ppm, ax=None):
    """
    Annotate a given spectrum with a given sequence of theoretical
    mz and descriptions

    :params: mz2intens: dict mapping experimental mz to intensity
    :params: ion2desc: dict mapping theoretical mz to description
    :params: ppm: error allowed during assign in parts-per-million
    :params: axis (optional): axis to plot the figure on, default currrent axis

    :return: assignment_error: list of relative errors in ppm for the assignment
    """

    mz_list = sorted(mz2intens.keys())
    ion_list = sorted(ions2desc.keys())

    if ax is None:
        ax = plt.gca()

    last = 0

    def ColorFromDesc(desc_list):
        """
        Return the plink color scheme for a given ion description
        """
        if desc_list[1] == 'A':
            return 'yellow'
        elif desc_list[1] == 'a':
            return 'brown'
        elif desc_list[1] == 'B':
            return 'green'
        elif desc_list[1] == 'b':
            return 'orange'
        elif desc_list[1] == 'Y':
            return 'red'
        elif desc_list[1] == 'y':
            return 'purple'
        elif 'M' in  desc_list[1]:
            return 'cyan'

    assignment_error = []

    for mz in mz_list:
        # store axis-object in variable for return
        ax.plot([mz, mz], # x1, x2
                [0, mz2intens[mz]], # y1, y2
                'k-', # style
                lw=1 # line width
                )
        for idx, ion in enumerate(ion_list[last:]):
            if ion >= mz * (1+ppm[0]/10**6):
                if ion <= mz * (1+ppm[1]/10**6):
                    # add the relative error of the assignments to the list
                    assignment_error.append((ion-mz)/ion * 10**6)
                    ax.plot([ion, ion],
                               [0, mz2intens[mz]],
                               ColorFromDesc(ions2desc[ion]),
                               lw=2)
                    ax.text(ion,
                              mz2intens[mz]+1000,
                              '${0}_{{{1}}}^{{{2}+}}$'.format(ions2desc[ion][1],
                                                              ions2desc[ion][2],
                                                              ions2desc[ion][3]),
                              {'ha': 'left', 'va': 'bottom'},
                              rotation=90,
                              color=ColorFromDesc(ions2desc[ion]))
                    last += idx

    ax.set_xlim(min(mz_list),max(mz_list))
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

if __name__ == '__main__':

    fig, ax = plt.subplots(2)

    ppm = [-50, 50]

    assignment_error = []

#    ions2desc =  IonsFromXlinkSequence('DQKLSELDDR', # peptid1
#                                    3, # xlink1
#                                    'KICEGFR', # peptide2
#                                    1, # xlink2
#                                    138.068, #mass of xlinker
#                                    ['a', 'b', 'y'], # ion types
#                                    [1,2], # min, max charge
#                                    'monoisotopic', # mass type
#                                    [[[], []], [[57.021464,], [2,]]], # modifications (on pep1 and pep2)
#                                    2000) # max mass


    ions2desc =  IonsFromXlinkSequence('KAWGNNQDGVVASQPAR', # peptid1
                                    1, # xlink1
                                    'LKSSDAYK', # peptide2
                                    2, # xlink2
                                    138.068, #mass of xlinker
                                    ['a', 'b', 'y'], # ion types
                                    [1,3], # min, max charge
                                    'monoisotopic', # mass type
                                    None, # modifications
                                    2000) # max mass

    with open('peptides.log', 'w') as f:
        f.write('\t'.join(['Mass', 'Peptide', 'IonType', 'Length', 'Charge', 'PepType', 'Mods']) + '\n')
        for idx, i in enumerate(ions2desc.keys()):
            f.write('{}\t{}\n'.format(i, '\t'.join([str(i) for i in ions2desc[i]])))

    mgfPath = os.path.join(binPath,
                           '../testdata/SV_plink/2017_08_04_SVs_BS3_16.mgf')
    with open(mgfPath, 'r') as f:
        spectrum2offset = Reader.IndexMGF(f)

        mz2intens = Reader.ReadSpectrum(18452, # spectrum
                                        f, # file handle
                                        spectrum2offset) # spectrum dict

        assignment_error = AssignAndPlotPSM(mz2intens, ions2desc, ppm, ax=ax[0])

        PlotHist(assignment_error, ppm, ax[1])

        plt.show()