**CroCo** converts multiple data format from cross-linking mass
spectrometry software tools to xTable format (in csv format). Currently
the following input formats are supported:

-  `Kojak`_
-  `Kojak with Percolator`_
-  `StavroX`_
-  `Xi`_
-  Xi + XiFDR
-  pLink1 (discontinued)
-  `pLink2`_
-  `xQuest`_
-  xTable

Input files can be converted into different formats typically used for
analysis of cross-linking data (e.g. visualisation, spectra annotation,
…). The following formats are supported:

-  `DynamXL`_
-  customTable
-  `pLabel`_
-  xTable
-  `xVis`_
-  `xWalk`_
-  `xiNet`_

CroCo is distributed as graphical program to be run from an executable
and as a Python module to be integrated into workflows.

System requirements
-------------------

*For the GUI*: Windows 10

*For the Python module*:  Python3 with the following modules installed
  - pandas
  - numpy
  - re

Usage of the GUI
----------------

For the conversion of data of every input program, a slightly different
usage is required for gathering all data that are required for xTable.
In general, information that is not present in the input files will be
asked from the user.

For detailed information how to use the GUI, see :ref:`introGUI`.

.. _Kojak: http://www.kojak-ms.org/
.. _Kojak with Percolator: https://github.com/percolator
.. _StavroX: https://www.stavrox.com/
.. _Xi: https://github.com/Rappsilber-Laboratory/XiSearch
.. _pLink2: http://pfind.ict.ac.cn/software/pLink/index.html
.. _xQuest: http://proteomics.ethz.ch/cgi-bin/xquest2_cgi/
.. _DynamXL: https://degiacomi.org/software/dynamxl/
.. _pLabel: http://pfind.ict.ac.cn/software/pLabel/index.html
.. _xVis: https://xvis.genzentrum.lmu.de/login.php
.. _xWalk: http://www.xwalk.org/cgi-bin/home.cgi
.. _xiNet: http://crosslinkviewer.org/