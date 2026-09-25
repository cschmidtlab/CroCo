# CroCo

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub tag](https://img.shields.io/github/tag/cschmidtlab/croco.svg)](https://GitHub.com/cschmidtlab/croco/tags/)
[![Github all releases](https://img.shields.io/github/downloads/cschmidtlab/croco/total.svg)](https://GitHub.com/cschmidtl/releases)

**CroCo** converts multiple data format from cross-linking mass spectrometry software tools to xTable format (in csv format).


## Important links

You can find the **latest release** including compiled versions at https://github.com/cschmidtlab/CroCo/releases/latest/.

The full documentation for Croco is available at https://cschmidtlab.github.io/CroCo/.

## Supported formats

Currently the following input formats are supported:

  - [Kojak](http://www.kojak-ms.org/)
  - [Kojak with Percolator](https://github.com/percolator)
  - [StavroX](https://www.stavrox.com/)
  - [Xi](https://github.com/Rappsilber-Laboratory/XiSearch)
  - Xi + XiFDR
  - pLink1 (discontinued)
  - pLink2 (discontinued)
  - [pLink3](https://pfind.ict.ac.cn/se/plink/) (via the pLink2 module)
  - [xQuest](http://proteomics.ethz.ch/cgi-bin/xquest2_cgi/)
  - xTable

Input files can be converted into different formats typically used for analysis of cross-linking data (e.g. visualisation, spectra annotation, ...).
The following formats are supported:

  - [DynamXL](https://degiacomi.org/software/dynamxl/)
  - customTable
  - [pLabel](http://pfind.ict.ac.cn/software/pLabel/index.html)
  - xTable
  - [xVis](https://xvis.genzentrum.lmu.de/login.php)
  - [xWalk](http://www.xwalk.org/cgi-bin/home.cgi)
  - [xiNet](http://crosslinkviewer.org/)

CroCo is distributed as graphical program to be run from an executable and as a Python module to be integrated into workflows.

## System requirements
*For the GUI*:
  * Windows 10

*For the Python module*:
  * Python3 with the following modules installed
    * pandas
    * numpy
    * re
  * Install the dependencies declared in [`pyproject.toml`](pyproject.toml), e.g. `pip install .`

## Usage
For the conversion of data of every input program, a slightly different usage is required for gathering all data that are required for xTable.
In general, information that is not present in the input files will be asked from the user.

## Compiling the Windows binaries
The standalone Windows GUI executable (`croco_wx.exe`) is built from `src/croco_wx.py` with [PyInstaller](https://pyinstaller.org/).
The build must be performed **on Windows** (wxPython does not ship wheels for other platforms).

  - **Quick build:** run `build_windows.bat` in the repository root. It creates a virtual environment, installs the `gui` and `build` extras declared in `pyproject.toml` and produces `dist/croco_wx.exe`.
  - **Manual build:** see [docs/build_windows.rst](docs/build_windows.rst), e.g. `pip install ".[gui,build]"` followed by `pyinstaller croco_wx_single.spec`; a reproducible conda environment is provided in `environment.yaml`.

## Versioning
The version is derived automatically from the git tags with [hatch-vcs](https://github.com/ofek/hatch-vcs); there is no version number to maintain by hand.
Builds between tags get a development version such as `0.7.2.dev25` (the tag `0.7.1` plus 25 commits). To release a new version, create a tag, e.g. `git tag 0.8.0`.
The GUI About dialog and the documentation read the version from the package (`croco.__version__`), so it stays consistent everywhere.

### Input formats

#### Kojak
  * **Load file(s)**: e.g. `FILENAME.kojak.txt`
  * Provide: Rawfile title (e.g. `FILENAME.raw`)

#### Kojak & Percolator
For this script to work, the unpercolated Kojak file (e.g. `FILENAME.kojak.txt`) has to be in the same directory as the percolated file.

  * **Load file(s)**: e.g. `FILENAME.validated.txt`
  * Provide: Rawfile title (e.g. `FILENAME.raw`)

#### StavroX
  * **Load file(s)**: StavroX results file (e.g. FILENAME.csv)
  * Provide: Path to SSF-file

#### Xi
  * **Load file(s)**: Path to Xi results file (e.g. `FILENAME_XiVersion1.6.739.csv`)

#### Xi & XiFDR
  * **Load file(s)**: Path to xiFDR file (e.g. `results_Xi1.8.7_CSM_xiFDR2.3.5.csv`)
  * Provide: Path to corresponding Xi config file (e.g. `results_Xi1.8.7.config`)

#### pLink1
  * **Load file(s)** (folder): sample folder within the pLink results dir (e.g. `2.report\sample1`)

#### pLink2
  * **Load file(s)** (folder): reports folder within the pLink results (e.g. `pLink_task_2018.06.12.09.33.10\reports`)

#### xQuest
  * **Load file(s)**: xQuest results file exported as csv (e.g. `FILENAME_xquest.csv`)

### Output formats

#### DynamXL
  * **Write to**: Directory in which to save the DynamXL file

#### customTable
  * **Write to**: Directory in which to save the customTable csv file
  * Provide: customTable template file

##### customTable Format
```
[header]
Protein 1, Protein 2
[data]
[prot1], [prot2]
[footer]
This is the footer of the file
```
  * Everything between [header] and [data] is considered the header and is printed once on top of the output file
  * In the [data] block columns identified by their xTable header are written in substitution of the header name.
    * E.g. instead of `[prot1], [prot2]` a line like `SPA_STAAU,  IgG4_heavy` is written for every line in the xTable file
    * Everything **not** enclosed in brackets will be written as is
    * If an invalid header name is given in brackets, the program will stop

#### pLabel
  * **Write to**: Directory in which to save the pLabel file
  * Provide: Directory containing the corresponding mgf-files
    * The mgf filenames **must** match the rawfile names given in the xTable

#### xTable
  * **Write to**: Directory in which to save xTable file

#### xVis
  * **Write to**: Directory in which to save xVis file

#### xWalk
  * **Write to**: Directory in which to save xWalk file
  * Provide:
    * PDB to map xlinks: Name of a PDB file that should be analysed by xWalk
    * PDB Atom code (Text): A PDB atom code (e.g. CB) that should be used for distance calculation

#### xiNet
  * **Write to**: Directory in which to save xiNez file


## Version History 

#### 0.5
  * customTable support

#### 0.4
  * Added pLabel support
  * Options Window on GUI
