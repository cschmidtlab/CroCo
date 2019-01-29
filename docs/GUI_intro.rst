.. _introGUI:

The CroCo graphical User-Interface
==================================

Input formats
~~~~~~~~~~~~~

Kojak
^^^^^

-  **Load file(s)**: e.g. \ ``FILENAME.kojak.txt``
-  Provide: Rawfile title (e.g. ``FILENAME.raw``)

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

\``\` [header] Protein 1, Protein 2 [data] [prot1], [pr
