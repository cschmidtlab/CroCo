.. _crocowrite:

Functions documentation for write functions
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

xTable
------


-  **Write to**: Directory in which to save xTable file

.. automodule:: croco.xTable
   :members:

xWalk
-----


-  **Write to**: Directory in which to save xWalk file
-  OptionsWindow:

   -  PDB to map xlinks: Name of a PDB file that should be analysed by
      xWalk
   -  Offset: shift between PDB AA indices and the xTable
   -  Chains: comma separated list of protein:chain allocations, e.g. ProteinA:AB, ProteinB:C
   -  PDB Atom code (Text): A PDB atom code (e.g. CB) that should be used for distance calculation

.. automodule:: croco.xWalk
   :members:

xVis
----
-  **Write to**: Directory in which to save xVis file

.. automodule:: croco.xVis
   :members:

customTable
-----------

-  **Write to**: Directory in which to save the customTable csv file
-  OptionsWindow: customTable template file

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

.. automodule:: croco.customTable
   :members:

DynamXL
-------


-  **Write to**: Directory in which to save the DynamXL file

.. automodule:: croco.DynamXL
   :members:

xiNET
-----


-  **Write to**: Directory in which to save xiNet file


.. automodule:: croco.xiNET
   :members:

pLabel
------


-  **Write to**: Directory in which to save the pLabel file
-  OptionsWindow: Directory containing the corresponding mgf-files
   - The mgf filenames **must** match the rawfile names given in the xTable
- OptionsWindow: ``Merge mgf files`` will generate a new mgf-files containing only those spectra that were references in the xTable. It will adopt all paths in the pLabel files to match the location of the merged mgf-file.


.. automodule:: croco.pLabel
   :members:
