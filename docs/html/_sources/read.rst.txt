Modules to Read cross-link datafiles
************************************

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Kojak
=====
CroCo contains modules for reading data from the direct output of the Kojak search engine (http://www.kojak-ms.org) as well as for FDR-filtered results obtained by using Percolator.
A detailed workflow how to obtain the FDR-filtered table can be found at http://www.kojak-ms.org/docs/validation.html.

Kojak only
----------

.. automodule:: croco.Kojak
   :members:

Kojak and Percolator
--------------------
   
.. automodule:: croco.KojakPercolator
   :members:

Kojak Helper Functions
----------------------

.. automodule:: croco.KojakFunctions
   :members:

   
pLink
=====

CroCo can read data from the search engine pLink2 (http://pfind.ict.ac.cn/software/pLink1/index.html) as well as from its predecessor pLink1. Even though the latter is deprecated since 2018, we provide the parser for compatibility to old analyses.

pLink1
------
.. automodule:: croco.pLink1
   :members:

   
pLink2
------
.. automodule:: croco.pLink2
   :members:

StavroX
=======

Currently CroCo can only handle data from the cross-link search engine StavroX (https://www.stavrox.com/) that cannot target MS-cleavable cross-linkers.
Support for MeroX (the analogue of StavroX that targets cleavable cross-linker) may be supported in future versions of CroCo.

.. automodule:: croco.StavroX
   :members:

Xi
==

Results from the Xi cross-link search engine (http://rappsilberlab.org/rappsilber-laboratory-home-page/tools/) can also be parsed by CroCo.

.. automodule:: croco.Xi
   :members:

.. automodule:: croco.XiSearchFDR
   :members:

xQuest
======

xQuest (http://prottools.ethz.ch/orinner/public/htdocs/xquest/) results are also supported by CroCo.

.. automodule:: croco.xQuest
   :members:
