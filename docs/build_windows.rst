.. _buildwindows:

Compiling the CroCo Windows binaries
====================================

CroCo's graphical user interface is distributed as a standalone Windows
executable (``croco_wx.exe``) built with PyInstaller. This page explains
how to build that executable from the source code.

.. note::

   The build must happen **on Windows** (64-bit). wxPython does not ship
   wheels for other operating systems and PyInstaller does not
   cross-compile, so building a Windows binary on Linux or macOS is not
   supported.

Prerequisites
-------------

* 64-bit Windows 10 or Windows 11
* Either uv installed from https://docs.astral.sh/uv/ (recommended: it can
  provide Python itself), or Python 3.11 or 3.12 installed from
  https://www.python.org/ (tick ``Add python.exe to PATH`` during installation)
* Git (optional) for cloning the repository
* A working internet connection to download the dependencies

Quick build (recommended)
-------------------------

#. Get the source code::

       git clone https://github.com/cschmidtlab/CroCo.git
       cd CroCo

   or download and extract the latest source archive from
   https://github.com/cschmidtlab/CroCo/archive/master.zip

#. Double-click ``build_windows.bat`` (or run it from a command prompt).
   The script uses ``uv`` when it is installed (letting uv provide Python
   3.12), otherwise it falls back to a locally installed Python 3.11/3.12.
   It then creates a virtual environment (``.venv``), installs the ``gui``
   and ``build`` extras from ``pyproject.toml`` and runs PyInstaller.

#. Take the finished binary from::

       dist\croco_wx.exe

Manual build
------------

Create an isolated environment and install the dependencies. Either use
``venv`` + pip::

    py -3.12 -m venv .venv
    .venv\Scripts\activate
    python -m pip install -U pip
    pip install ".[gui,build]"

or conda (using the ``conda-forge`` channel)::

    conda env create -f environment.yaml
    conda activate croco_env

Then bundle the GUI with PyInstaller::

    pyinstaller --noconfirm --clean croco_wx_single.spec

The single-file executable is written to ``dist\croco_wx.exe``.

Optional: faster installs with uv
---------------------------------

If `uv <https://docs.astral.sh/uv/>`_ is installed, the ``venv`` + pip
step can be replaced with::

    uv sync --extra gui --extra build
    uv run pyinstaller --noconfirm --clean croco_wx_single.spec

What the build does
-------------------

* ``croco_wx_single.spec`` bundles ``src/croco_wx.py`` and the ``src/croco``
  package into a **one-file** executable. Data files (e.g.
  ``src/croco/data/croco_logo.ico``) are packed into the exe and unpacked
  to a temporary directory at runtime (``sys._MEIPASS``).
* Dependencies are declared in ``pyproject.toml``. The ``gui`` extra
  provides wxPython and the ``build`` extra provides PyInstaller, so
  ``pip install ".[gui,build]"`` prepares everything needed for the bundle.
* pandas and openpyxl are collected automatically by the hooks that
  ship with PyInstaller >= 6.
* Since the executable is unsigned, Windows SmartScreen will warn on the
  first start. Acknowledge and continue, or sign the binary with a
  code-signing certificate before distributing it.

Notes on the modernized build
-----------------------------

The build originally required a conda environment pinned to Python 3.6,
a pre-release wxPython and PyInstaller 3.x, and the spec file contained
absolute paths of a single developer machine. The current procedure:

* Uses Python 3.12 with released wxPython 4.2+ wheels and PyInstaller >= 6.
* Ships a spec file whose paths are resolved relative to the spec file
  location (``SPECPATH``), so it builds on any machine without editing.
* Disables UPX compression by default to avoid antivirus false positives;
  set ``upx=True`` in ``croco_wx_single.spec`` to re-enable it.
* Keeps the dependency list in ``pyproject.toml`` (extras ``gui`` and
  ``build``) as the single source of truth; ``environment.yaml`` is still
  available for conda users, and ``build_windows.bat`` provides a
  one-command build.

Troubleshooting
---------------

* ``'pyinstaller' is not recognized`` - the environment was not activated;
  run the PyInstaller command from the activated ``.venv`` / ``croco_env``.
* ``Python was not found; run without arguments to install from the
  Microsoft Store`` - no system Python is available (or the Store alias is
  enabled). Install Python from python.org, or install ``uv`` so that
  ``build_windows.bat`` provisions Python automatically.
* wxPython fails to install - make sure you use Python 3.11 or 3.12
  (64-bit), the versions for which official wxPython wheels exist.
* The executable does not start on a target machine - install the
  Microsoft Visual C++ Redistributable (x64), on which some bundled
  libraries depend.
* A console window opens next to the GUI - set ``console=False`` in
  ``croco_wx_single.spec`` to start the GUI without the console window.