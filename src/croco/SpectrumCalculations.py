# -*- coding: utf-8 -*-
"""
Library of functions associated with ion mass calculation fram a sequence
Created on Wed Nov 22 09:26:33 2017

@author: User
"""

############################################
# Calculation of fragment peptide sequences
############################################

def NTermPeptides(parent, start=1, end=None):
    """
    return all n-terminal peptides from start (including)
    to end (non-including)
    """
    if end == None:
        end = len(parent)

    nterm = [parent[0:j] for j in range(start, end)]
    return nterm


def CTermPeptides(parent, start=1, end=None):
    """
    return all c-terminal peptides from start (non-including)
    to end (including)
    """
    if end == None:
        end = len(parent) 
    cterm = [parent[i:len(parent)] for i in range(start, end)]
    return cterm

def XPeptides(peptide1, xlink1, peptide2, xlink2, CTerm=False):
    """
    Calculate all fragmentation ion sequences for a cross-link peptide
    and return them sorted by whether they contain the cross-link or not
    
    :params: peptide1: sequence of the alpha peptide of the cross-link
    :params: xlink1: relative position of the cross-link within peptide1
    :params: peptide2: sequence of the beta peptide of the cross-link ion
    :params: xlink2:  relative position of the cross-link within peptide2
    :params: CTerm: Calc Cterminal peptides if true, else Nterminal
    """
    if CTerm:
        """
        Return all c-terminal peptides from a cross-linked sequence
        """
        # truncation of peptide 1
        unmod1 = CTermPeptides(peptide1, start=xlink1)
        mod1 = CTermPeptides(peptide1, end=xlink1)
        # return the indices of the beginning and the end of the substrings
        unmod1_pos = [xlink1+1, len(peptide1)]
        mod1_pos = [2, len(peptide1)]
        
        # truncation of peptide2
        unmod2 = CTermPeptides(peptide2, start=xlink2)
        mod2 = CTermPeptides(peptide2, end=xlink2)
        unmod2_pos = [xlink1+1, len(peptide2)]
        mod2_pos = [2, len(peptide2)]

    else:
        """
        Return all n-terminal peptides from a cross-linked sequence
        """
        # truncation of peptide1
        unmod1 = NTermPeptides(peptide1, end=xlink1)
        mod1 = NTermPeptides(peptide1, start=xlink1)
        # return the indices of the beginning and the end of the substrings
        unmod1_pos = [1, xlink1-1]
        mod1_pos = [1, len(peptide1)-1]
        
        # truncation of peptide2
        unmod2 = NTermPeptides(peptide2, end=xlink2)
        mod2 = NTermPeptides(peptide2, start=xlink2)
        unmod2_pos = [1, xlink2-1]
        mod2_pos = [xlink2, len(peptide2)]
        
    # return positions as single list of lists
    positions = [unmod1_pos,
                 mod1_pos,
                 unmod2_pos,
                 mod2_pos]
    
    return unmod1, mod1, unmod2, mod2, positions

############################################
# Calculation of masses of sequences
############################################

def BIonMass(peptide, aa_dict, charge, mods, peptide_type=None):
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

def AIonMass(peptide, aa_dict, charge, mods, peptide_type=None):
    ion_string = 'A' if peptide_type == 'alpha' else 'a'
    mass = 0
    for aa in peptide:
        mass += aa_dict[aa]
    mass += -28.0101 + charge * 1.00794 # M + nH+
    mass /= charge # charging effect
    desc =  [peptide, ion_string, len(peptide), charge,
             peptide_type, mods]
    return mass, desc

def YIonMass(peptide, aa_dict, charge, mods, peptide_type=None):
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

def ParentionMass(peptide, aa_dict, charge, mods, peptide_type=None):
    ion_string = 'M+{}H'.format(charge)
    mass = 0
    for aa in peptide:
        mass += aa_dict[aa]
    mass += 18.01528 + charge * 1.00794 # M + nH+
    mass /= charge
    desc = [peptide, ion_string, '', charge,
            peptide_type, mods]
    return mass, desc

def XParentionMass(pep1_mass, pep2_mass, peptide1, xlink1, peptide2,
                    xlink2, charge, mods, peptide_type=None):
    ion_string = 'M+{}H'.format(charge)
    mass = sum((pep1_mass, pep2_mass))
    mass += charge * 1.00794 # M + nH+
    mass /= charge
    peptide = '-'.join([str(x) for x in [peptide1, xlink1, peptide2, xlink2]])
    desc = [peptide, ion_string, '', charge,
            peptide_type, mods]
    return mass, desc
