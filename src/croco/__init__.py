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
from . import pLink1
pLink1.init(plink_col_order)
print('Finished with pLink1: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import pLink2
pLink2.init(plink_col_order)
print('Finished with pLink2: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import Kojak
from . import KojakPercolator
Kojak.init(col_order)
KojakPercolator.init(col_order)
print('Finished with Kojak: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import Xi
from . import XiSearchFDR
Xi.init(col_order)
XiSearchFDR.init(col_order)
print('Finished with Xi: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import xQuest
xQuest.init(col_order)
print('Finished with xQuest: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import StavroX
StavroX.init(col_order)
print('Finished with StavroX: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import HelperFunctions

from . import DynamXL

from . import xiNET

from . import xTable
xTable.init(col_order)

from . import xVis

from . import xWalk

from . import pLabel

from . import customTable
