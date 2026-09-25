# -*- mode: python -*-
#
# PyInstaller spec file for the CroCo Windows GUI binary.
#
# Path-independent: always run PyInstaller from the repository root, i.e.
#     pyinstaller --noconfirm croco_wx_single.spec
# The variable SPECPATH is provided by PyInstaller and points to the directory
#that contains this spec file, so no hard-coded machine paths are needed.

import os
import sys

sys.setrecursionlimit(5000)

ROOT = SPECPATH

block_cipher = None

a = Analysis(['src/croco_wx.py'],
             pathex=[os.path.join(ROOT, 'src')],
             binaries=[],
             datas=[(os.path.join(ROOT, 'src', 'croco', 'data'), 'data')],
             # pandas and openpyxl submodules are collected automatically by the
             # hooks shipped with modern PyInstaller (>=6.0), so no hidden
             # imports have to be listed here any more.

             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=['PyQt5', 'matplotlib'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
           a.scripts,
           a.binaries,
           a.zipfiles,
           a.datas,
           name='croco_wx',
           icon=os.path.join(ROOT, 'artwork', 'croco_logo.ico'),
           debug=False,
           strip=False,
           # UPX is disabled because it triggers antivirus false positives;
           # set upx=True to make the binary smaller.

           upx=False,
           runtime_tmpdir=None,
           # set console=False to hide the console window of the GUI.
           console=True)