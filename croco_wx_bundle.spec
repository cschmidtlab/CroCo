# -*- mode: python -*-

block_cipher = None


a = Analysis(['src\\croco_wx.py'],
             pathex=['C:\\Users\\User\\Documents\\03_software\\python\\CroCo'],
             binaries=[],
             datas=[('src\\croco\\data', 'data')],
             hiddenimports=['pandas._libs.tslibs.timedeltas', 'openpyxl'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['PyQt5'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='croco_wx',
          debug=False,
          strip=False,
          upx=True,
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='croco_wx')
