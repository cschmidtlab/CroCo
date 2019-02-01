.. _introGUI:

The graphical User-Interface
==================================

A portable GUI to easily transfer cross-link data files and help users focusing on experimental work is the main idea behind CroCo.
However, usage of the modules inside Python-scripts is possible (and documented at :ref:`crocoread` and :ref:`crocowrite`).

General usage
=============

 - Download the compiled GUI from https://github.com/cschmidtlab/CroCo/releases/latest
    - If your OS is not supported you can still try `Starting the CroCo GUI from Python`_.
 - Start the programme
    - If Windows complains about the non-verified origin of the programme, acknowledge and continue
 - You should see something like the following

 .. figure:: _static/img/crocoGUI.png
    :align: center
    :alt: The CroCO GUI
	
    The CroCo GUI

 - In the Input and output dropdown menus, choose a format and select a directory or one or multiple files (see `Detailed descriptions for the formats`_.
    - the Start button will only be available once both formats are set and files are selected
 - If you selected multiple files for input, you can merge them by activating the ``Merge files`` checkbox
 - CroCo transfers all information found in the input files into the generated xTable files. If you want to restrict the xTable headers to the minimum columns necessary, you can choose ``Compact xTable``
 - Click on Start
 - Depending on your choice of input and output formats, CroCo may require additional information. This will pop up in a second window.
    - Detailed information on how to enter format-specific information can be found at `Detailed descriptions for the formats`_.



Detailed descriptions for the formats
=====================================

Input formats
~~~~~~~~~~~~~

Kojak
^^^^^

-  **Load file(s)**: e.g. \ ``FILENAME.kojak.txt``
-  Provide: Rawfile title (e.g. ``FILENAME.raw``)

Kojak does not store the rawfile inside the output but it is required inside the xTable. 

Kojak & Percolator
^^^^^^^^^^^^^^^^^^

For this script to work, the unpercolated Kojak file
(e.g. ``FILENAME.kojak.txt``) has to be in the same directory as the
percolated file.

-  **Load file(s)**: e.g. \ ``FILENAME.validated.txt``
-  Provide: Rawfile title (e.g. ``FILENAME.raw``)

StavroX
^^^^^^^

-  **Load file(s)**: StavroX results file (e.g. FILENAME.csv)
-  Provide: Path to SSF-file

Xi
^^

-  **Load file(s)**: Path to Xi results file
   (e.g. ``FILENAME_XiVersion1.6.739.csv``)

Xi & XiFDR
^^^^^^^^^^

-  **Load file(s)**: Path to xiFDR file
   (e.g. ``FILENAME_5_FDR_PSM_xiFDR1.0.22.csv``)
-  Provide: Path to corresponding Xi results file
   (e.g. ``FILENAME_XiVersion1.6.739.csv``)

pLink1
^^^^^^

-  **Load file(s)** (folder): sample folder within the pLink results dir
   (e.g. ``2.report\sample1``)

pLink2
^^^^^^

-  **Load file(s)** (folder): reports folder within the pLink results
   (e.g. ``pLink_task_2018.06.12.09.33.10\reports``)

xQuest
^^^^^^

-  **Load file(s)**: xQuest results file exported as csv
   (e.g. ``FILENAME_xquest.csv``)

Output formats
~~~~~~~~~~~~~~

DynamXL
^^^^^^^

-  **Write to**: Directory in which to save the DynamXL file

customTable
^^^^^^^^^^^

-  **Write to**: Directory in which to save the customTable csv file
-  Provide: customTable template file

customTable Format
''''''''''''''''''

::

   [header]
   Protein 1, Protein 2
   [data]
   [prot1], [prot2]
   [footer]
   This is the footer of the file

-  Everything between [header] and [data] is considered the header and
   is printed once on top of the output file
-  In the [data] block columns identified by their xTable header are
   written in substitution of the header name.

   -  E.g. instead of ``[prot1], [prot2]`` a line like
      ``SPA_STAAU,  IgG4_heavy`` is written for every line in the xTable
      file
   -  Everything **not** enclosed in brackets will be written as is
   -  If an invalid header name is given in brackets, the program will
      stop

pLabel
^^^^^^

-  **Write to**: Directory in which to save the pLabel file
-  Provide: Directory containing the corresponding mgf-files

   -  The mgf filenames **must** match the rawfile names given in the
      xTable

xTable
^^^^^^

-  **Write to**: Directory in which to save xTable file

xVis
^^^^

-  **Write to**: Directory in which to save xVis file

xWalk
^^^^^

-  **Write to**: Directory in which to save xWalk file
-  Provide:

   -  PDB to map xlinks: Name of a PDB file that should be analysed by
      xWalk
   -  PDB Atom code (Text): A PDB atom code (e.g. CB) that should be
      used for distance calculation

xiNet
^^^^^

-  **Write to**: Directory in which to save xiNet file

Version History
---------------

**0.5**

-  customTable support

.. _section-1:

**0.4**


-  Added pLabel support
-  Options Window on GUI

.. _pythonGUI:

Starting the CroCo GUI from Python
==================================
 - Clone the CroCo repository from https://github.com/cschmidtlab/CroCo.git
    - If you're not familiar with github, download and extract the latest version of the CroCo source code from https://github.com/cschmidtlab/CroCo/archive/master.zip
 - Install Python 3+
    - The easiest way is usually to use anaconda Python (https://www.anaconda.com/distribution/)
    - Install the required python packages for CroCo (see requirements.txt in the CroCo root directory) either via ``pip install PACKAGE`` or via ``conda install PACKAGE``
 - Open a Python-aware command-line
    - On Windows, use Anaconda Prompt or cmd.exe if you exported Python to your OS PATH-Variable
 - Navigate to ``CroCo/src``
 - Start the GUI by typing ``python croco_wx.py``