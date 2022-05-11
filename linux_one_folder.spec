# -*- mode: python ; coding: utf-8 -*-

add_binaries = []
add_data = [('include_linux/freeda_logo.png','.'), ('include_linux/bedtools', 'bedtools'), ('include_linux/lib','.'), ('include_linux/bin', '.')]
add_imports = ['PIL._tkinter_finder']

block_cipher = None


a = Analysis(['freeda_pipeline_GUI.py'],
             pathex=[],
             binaries=add_binaries,
             datas=add_data,
             hiddenimports=add_imports,
             hookspath=[],
             hooksconfig={},
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts, 
          [],
          exclude_binaries=True,
          name='freeda_pipeline_GUI',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas, 
               strip=False,
               upx=True,
               upx_exclude=[],
               name='freeda_pipeline_GUI')
