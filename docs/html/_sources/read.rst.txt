.. _crocoread:

Functions documentation for read functions
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Kojak
-----
CroCo contains modules for reading data from the direct output of the Kojak search engine (http://www.kojak-ms.org) as well as for FDR-filtered results obtained by using Percolator.
A detailed workflow how to obtain the FDR-filtered table can be found at http://www.kojak-ms.org/docs/validation.html.

Kojak only
~~~~~~~~~~

-  **Load file(s)**: e.g. \ ``FILENAME.kojak.txt``
-  OptionsWindow: Rawfile title (e.g. ``FILENAME.raw``)

Kojak does not store the rawfile inside the output but it is required inside the xTable. 

.. automodule:: croco.Kojak
   :members:

Kojak and Percolator
~~~~~~~~~~~~~~~~~~~~

For this script to work, the unpercolated Kojak file
(e.g. ``FILENAME.kojak.txt``) has to be in the same directory as the
percolated file.

-  **Load file(s)**: e.g. \ ``FILENAME.validated.txt``
-  OptionsWindow: Rawfile title (e.g. ``FILENAME.raw``)

.. automodule:: croco.KojakPercolator
   :members:

Kojak Helper Functions
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: croco.KojakFunctions
   :members:

   
pLink
-----

CroCo can read data from the search engine pLink2 (http://pfind.ict.ac.cn/software/pLink1/index.html) as well as from its predecessor pLink1. Even though the latter is deprecated since 2018, we provide the parser for compatibility to old analyses.

pLink1
~~~~~~
.. automodule:: croco.pLink1
   :members:

   
pLink2
~~~~~~
.. automodule:: croco.pLink2
   :members:

StavroX
-------

Currently CroCo can only handle data from the cross-link search engine StavroX (https://www.stavrox.com/) that cannot target MS-cleavable cross-linkers.
Support for MeroX (the analogue of StavroX that targets cleavable cross-linker) may be supported in future versions of CroCo.

.. automodule:: croco.StavroX
   :members:

Xi
--

Results from the Xi cross-link search engine (http://rappsilberlab.org/rappsilber-laboratory-home-page/tools/) can also be parsed by CroCo.

Xi only
~~~~~~~

-  **Load file(s)**: Path to Xi results file
   (e.g. ``FILENAME_XiVersion1.6.739.csv``)


.. automodule:: croco.Xi
   :members:

Xi & XiFDR
~~~~~~~~~~

-  **Load file(s)**: Path to xiFDR file
   (e.g. ``FILENAME_5_FDR_PSM_xiFDR1.0.22.csv``)
-  OptionsWindow: Path to corresponding Xi results file
   (e.g. ``FILENAME_XiVersion1.6.739.csv``)

   
.. automodule:: croco.XiSearchFDR
   :members:

xQuest
------

xQuest (http://prottools.ethz.ch/orinner/public/htdocs/xquest/) results are also supported by CroCo.


-  **Load file(s)**: xQuest results file exported as csv
   (e.g. ``FILENAME_xquest.csv``)

.. automodule:: croco.xQuest
   :members:
