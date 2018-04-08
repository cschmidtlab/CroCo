Rem Batch script to call all programmes for exe-compilation

Rem TCL and TK libraries have to be manually set when using Anaconda
set TCL_LIBRARY=C:\ProgramData\Anaconda3\tcl\tcl8.6
set TK_LIBRARY=C:\ProgramData\Anaconda3\tcl\tk8.6

Rem Actual call of the cx_freeze compiler
C:\ProgramData\Anaconda3\python.exe C:\ProgramData\Anaconda3\Scripts\cxfreeze^
	--base-name=Win32GUI^
	--include-path="src"^
	--target-dir="dist/win"^
	--include-modules="numpy.core._methods"^
	bin\croco_qt.py

Rem --base-name=Win32GUI^
