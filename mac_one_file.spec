# -*- mode: python ; coding: utf-8 -*-

add_binaries = []

add_data = [('include_mac_M1/entitlements.plist','.'),
			('include_mac_M1/freeda_logo.png', '.'), 
			('include_mac_M1/bedtools', 'bedtools'),
			('include_mac_M1/lib', '.'),
			('include_mac_M1/bin', '.')]
			
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
          a.binaries,
          a.zipfiles,
          a.datas,  
          [],
          name='freeda_pipeline_GUI',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=False,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=False,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None )
app = BUNDLE(exe,
             name='freeda_pipeline_GUI.app',
             icon='include_mac_M1/freeda_img4_icon.icns',
             bundle_identifier=None)
