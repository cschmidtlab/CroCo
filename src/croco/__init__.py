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

# all conversion scripts are imported as modules and initialised
from . import pLink1
print('Finished with pLink1: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import pLink2
print('Finished with pLink2: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import Kojak
from . import KojakPercolator
print('Finished with Kojak: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import Xi
from . import XiSearchFDR
print('Finished with Xi: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import xQuest
print('Finished with xQuest: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import StavroX
print('Finished with StavroX: {:4.2f} ms'.format((timeit.default_timer() - starttime) * 1000))
starttime = timeit.default_timer()

from . import HelperFunctions

from . import DynamXL

from . import xiNET

from . import xTable

from . import xVis

from . import xWalk

from . import pLabel

from . import customTable
