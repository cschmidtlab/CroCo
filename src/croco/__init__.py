#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions to read and write cross-link information from different sources.

===============================
Definition of the xtable format
===============================

rawfile - name of the corresponding rawfile
scanno - Scan the cross-link was identified from
prec_ch - charge of the precursor ion

prot1 - fasta-header of the protein corresponding to the first peptide
prot2 - same for the second peptide

pepseq1 - Sequence of the longer (alpha) peptide
pepseq2 - Sequence of the shorter (beta) peptide

xlink1 - Position of the cross-linker within the longer peptide
xlink2 - Position of the cross-linker within the shorter peptide

pos1 - Absolute postiton of the first AA of the first peptide
pos2 - Absolute postiton of the first AA of the second peptide

xpos1 - Absolute position of the cross-linker of the longer peptide
xpos2 - Absolute position of the cross-linker of the shorter peptide (only if interlink)

mod1 - name of the modification(s) of the first peptide (;-delimited string)
mnod2 - name of the modification(s) of the second pepetide (;-delimited string)

modpos1 - relative position of a modification within peptide 1 (;-delimited string)
modpos2 - same for peptide 2
modmass1 - mass of modification 1 (;-delimited string)
modmass2 - mass of modification 2

ID - Identifier for the position of a cross-link between two proteins
decoy - true or false

score - not normalised score
type - inter, intra, loop or mono
xtype - heavy or light label

==============================
Definition of the xinfo format
==============================

xlinker - Name of the used reagent
link_dist - nominally bridged distance in Angstrom

"""

import timeit

starttime = timeit.default_timer()


# defines the column headers required for xtable output
col_order = [ 'rawfile', 'scanno', 'prec_ch',
              'pepseq1', 'xlink1',
              'pepseq2', 'xlink2', 'xtype',
              'modmass1', 'modpos1', 'mod1',
              'modmass2', 'modpos2', 'mod2',
              'prot1', 'xpos1', 'prot2',
              'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

plink_col_order = ['Order',] + col_order

# all conversion scripts are imported as modules and initialised
from .lib import pLink1
pLink1.init(plink_col_order)
print('Finished with pLink1: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from .lib import pLink2
pLink2.init(plink_col_order)
print('Finished with pLink2: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from .lib import Kojak
from .lib import KojakPercolator
Kojak.init(col_order)
KojakPercolator.init(col_order)
print('Finished with Kojak: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from .lib import Xi
from .lib import XiSearchFDR
Xi.init(col_order)
XiSearchFDR.init(col_order)
print('Finished with Xi: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from .lib import xQuest
xQuest.init(col_order)
print('Finished with xQuest: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from .lib import StavroX
StavroX.init(col_order)
print('Finished with StavroX: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from .lib import Manual
Manual.init(col_order)

from .lib import DynamXL

from .lib import xiNET

from .lib import xTable

from .lib import xVis

from .lib import xWalk
