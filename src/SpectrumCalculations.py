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
 
def NTermXPeptides(peptide1, xlink1, peptide2, xlink2):
    """
    Return all n-terminal peptides from a cross-linked sequence
    """
    # truncation of peptide1
    nterm_unmod1 = NTermPeptides(peptide1, end=xlink1)
    nterm_mod1 = NTermPeptides(peptide1, start=xlink1)
    
    # truncation fo peptide2
    nterm_unmod2 = NTermPeptides(peptide2, end=xlink2)
    nterm_mod2 = NTermPeptides(peptide2, start=xlink2)
    
    return nterm_unmod1, nterm_mod1, nterm_unmod2, nterm_mod2
 

def CTermPeptides(parent, start=1, end=None):
    """
    return all c-terminal peptides from start (including)
    to end (non-including)
    """
    if end == None:
        end = len(parent) 
    cterm = [parent[i:len(parent)] for i in range(start, end)]
    return cterm

def CTermXPeptides(peptide1, xlink1, peptide2, xlink2):
    """
    Return all c-terminal peptides from a cross-linked sequence
    """
    # truncation of peptide 1
    cterm_unmod1 = CTermPeptides(peptide1, start=xlink1)
    cterm_mod1 = CTermPeptides(peptide1, end=xlink1)
    
    # truncation of peptide2
    cterm_unmod2 = CTermPeptides(peptide2, start=xlink2)
    cterm_mod2 = CTermPeptides(peptide2, end=xlink2)

    return cterm_unmod1, cterm_mod1, cterm_unmod2, cterm_mod2

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

def ParentionMass(peptide, aa_dict, charge, mods):
    mass = 0
    for aa in peptide:
        mass += aa_dict[aa]
    mass += 18.01528 + charge * 1.00794 # M + nH+
    mass /= charge
    desc = (peptide, '(M+{}H)'.format(charge), charge, mods)
    return mass, desc

def XParentionMass(pep1_mass, pep2_mass, peptide1, xlink1, peptide2,
                    xlink2, charge, mods):
    mass = sum((pep1_mass, pep2_mass))
    mass += charge * 1.00794 # M + nH+
    mass /= charge
    desc = ['{}-{}-{}-{}'.format(peptide1, xlink1, peptide2, xlink2),\
            '(M+{}H)'.format(charge), charge, mods]
    return mass, desc
