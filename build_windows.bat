@ECHO OFF
REM ============================================================================
REM  CroCo  -  build the CroCo Windows GUI binary with PyInstaller
REM
REM  Creates or reuses a local virtual environment and bundles the GUI into a
REM  single executable:  dist\croco_wx.exe
REM
REM  Requires Python 3.11/3.12 OR uv (https://docs.astral.sh/uv/). When uv is
REM  available it is used and provisions its own Python, so a system Python is
REM  not needed.
REM
REM  Usage:  double-click this file or run  build_windows.bat  from cmd
REM ============================================================================

setlocal EnableExtensions
pushd %~dp0

REM ---- prefer uv: it can provision a suitable Python itself --------------------
where uv >nul 2>nul
if %errorlevel%==0 goto :build_with_uv

REM ---- fall back to an installed Python 3 interpreter --------------------------
set "PY_CMD="
py -3.12 -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY_CMD=py -3.12"
if defined PY_CMD goto :have_python
py -3.11 -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY_CMD=py -3.11"
if defined PY_CMD goto :have_python
py -3 -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY_CMD=py -3"
if defined PY_CMD goto :have_python
python -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY_CMD=python"
if defined PY_CMD goto :have_python

ECHO ERROR: no working Python 3 interpreter was found.
ECHO.
ECHO Install Python 3.11/3.12 from https://www.python.org/downloads/
ECHO or install uv from https://docs.astral.sh/uv/ and run this script again.
popd
exit /b 1

:have_python
if not exist ".venv\Scripts\python.exe" (
    ECHO [pip 1/3] Creating virtual environment in .venv ...
    %PY_CMD% -m venv .venv
    if errorlevel 1 goto :error
)
call ".venv\Scripts\activate.bat"

ECHO [pip 2/3] Installing dependencies (gui + build extras from pyproject.toml)
python -m pip install --upgrade pip
if errorlevel 1 goto :error
pip install ".[gui,build]"
if errorlevel 1 goto :error

ECHO [pip 3/3] Building dist\croco_wx.exe with PyInstaller ...
pyinstaller --noconfirm --clean croco_wx_single.spec
if errorlevel 1 goto :error
goto :success

:build_with_uv
ECHO [uv 1/2] Creating environment with uv (Python 3.12) ...
uv sync --python 3.12 --extra gui --extra build
if errorlevel 1 goto :error

ECHO [uv 2/2] Building dist\croco_wx.exe with PyInstaller ...
uv run pyinstaller --noconfirm --clean croco_wx_single.spec
if errorlevel 1 goto :error
goto :success

:success
ECHO.
ECHO Success:  Windows binary written to:  %CD%\dist\croco_wx.exe
ECHO.
goto :end

:error
ECHO.
ECHO Build failed. See the output above for details.
popd
exit /b 1

:end
popd
endlocal