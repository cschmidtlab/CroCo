# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 12:59:06 2017

@author: User
"""

import csv
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
    with open('../config/aa_masses.txt', mode='r') as infile:
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

ions, desc =  IonsFromSequence('PEPTIDE',
                               ['a', 'b', 'y'],
                               [1,2],
                               'mono',
                               None)
   
for idx, ion in enumerate(ions):
    print(ion,': ', desc[idx])