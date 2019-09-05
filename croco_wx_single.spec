# -*- mode: python -*-
import sys
sys.setrecursionlimit(5000) # or more

block_cipher = None


a = Analysis(['src\\croco_wx.py'],
             pathex=['C:\\Users\\User\\Documents\\03_software\\python\\CroCo'],
             binaries=[],
             datas=[('src\\croco\\data', 'data')],
             hiddenimports=['pandas._libs.tslibs.timedeltas', 'pandas._libs.skiplist', 'openpyxl'],
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
          icon=r'C:\Users\User\Documents\03_software\python\CroCo\artwork\croco_logo.ico',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=True )
