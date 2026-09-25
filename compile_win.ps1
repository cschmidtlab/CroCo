#Requires -Version 5.1
<#
.SYNOPSIS
    Stage the CroCo sources in C:\Temp and build the Windows executable.

.DESCRIPTION
    Copies every file required by build_windows.bat (spec, project metadata,
    package sources and artwork) into a staging directory under C:\Temp, then
    runs build_windows.bat there. The finished binary ends up in
    <Destination>\dist\croco_wx.exe.

    build_windows.bat itself is location independent (it uses %~dp0), so the
    build can run from the staging copy instead of the repository checkout.

.PARAMETER Destination
    Staging directory. Defaults to C:\Temp\CroCo.

.PARAMETER Clean
    Remove an existing staging directory before copying. Without this switch
    the copy is done in place, which preserves a previously created .venv and
    speeds up subsequent builds.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\compile_win.ps1

.EXAMPLE
    .\compile_win.ps1 -Destination 'C:\Temp\CroCo' -Clean
#>
[CmdletBinding()]
param(
    [string]$Destination = 'C:\Temp\CroCo',
    [switch]$Clean
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message" -ForegroundColor Cyan
}

$SourceRoot = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($SourceRoot)) {
    $SourceRoot = (Get-Location).Path
}
Write-Step "Source root: $SourceRoot"

# Files and directories that build_windows.bat / the spec file need.
$RequiredItems = @(
    'build_windows.bat',
    'croco_wx_single.spec',
    'pyproject.toml',
    'README.md',
    'src',
    'artwork'
)

# Present in the repository but not strictly mandatory for the build.
$OptionalItems = @(
    'uv.lock',
    'environment.yaml',
    'xTable_definition.xlsx'
)

# Fail early if a required item is missing.
foreach ($item in $RequiredItems) {
    $path = Join-Path $SourceRoot $item
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required item not found: $path"
    }
}

if ($Clean -and (Test-Path -LiteralPath $Destination)) {
    Write-Step "Removing existing staging directory: $Destination"
    Remove-Item -LiteralPath $Destination -Recurse -Force
}
New-Item -ItemType Directory -Path $Destination -Force | Out-Null

function Copy-IntoStaging {
    param(
        [string]$RelativePath,
        [string[]]$ExcludeDirs = @()
    )

    $source = Join-Path $SourceRoot $RelativePath
    if (-not (Test-Path -LiteralPath $source)) {
        return
    }

    Write-Host "    copying $RelativePath"
    if (Test-Path -LiteralPath $source -PathType Container) {
        # robocopy copies the *contents* of the source, so aim it at the
        # matching subfolder of the destination to keep the tree layout.
        $target = Join-Path $Destination $RelativePath
        $robocopyArgs = @(
            $source, $target,
            '/E', '/NFL', '/NDL', '/NJH', '/NJS', '/NP', '/R:2', '/W:1'
        )
        if ($null -ne $ExcludeDirs -and $ExcludeDirs.Count -gt 0) {
            $robocopyArgs += '/XD'
            $robocopyArgs += $ExcludeDirs
        }
        & robocopy @robocopyArgs | Out-Null
        # robocopy exit codes 0-7 mean success, >= 8 indicates a failure.
        if ($LASTEXITCODE -ge 8) {
            throw "robocopy failed for '$RelativePath' (exit code $LASTEXITCODE)"
        }
    }
    else {
        Copy-Item -LiteralPath $source -Destination $Destination -Force
    }
}

# The build runs in a staging directory without git history, so extract the
# version from the source checkout here and hand it to the hatch-vcs build hook
# through the environment. If the version cannot be extracted, fail loudly.
$uv = Get-Command uv -ErrorAction SilentlyContinue
if (-not $uv) {
    throw "uv is required to extract the project version. " +
        "Install it from https://docs.astral.sh/uv/."
}

Write-Step "Extracting version from git"
$getVersionPy = "import sys; from setuptools_scm import get_version; " +
    "print(get_version(root=sys.argv[1], local_scheme='no-local-version'))"
$version = & $uv.Source run --quiet --no-project --with setuptools-scm `
    python -c $getVersionPy $SourceRoot | Select-Object -Last 1

if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($version)) {
    throw "Could not extract a version from '$SourceRoot'. " +
        "Make sure git is installed and the repository has tags."
}

$version = ([string]$version).Trim()
$env:SETUPTOOLS_SCM_PRETEND_VERSION = $version
$env:VCS_VERSIONING_PRETEND_VERSION = $version
Write-Step "Version: $version"

Write-Step "Copying build inputs to $Destination"
foreach ($item in ($RequiredItems + $OptionalItems)) {
    $excludes = @()
    if ($item -eq 'src') { $excludes = @('__pycache__') }
    Copy-IntoStaging -RelativePath $item -ExcludeDirs $excludes
}

$bat = Join-Path $Destination 'build_windows.bat'
Write-Step "Running build_windows.bat in $Destination"
$process = Start-Process -FilePath $bat `
    -WorkingDirectory $Destination `
    -NoNewWindow -Wait -PassThru
$exitCode = $process.ExitCode

if ($exitCode -ne 0) {
    throw "build_windows.bat failed with exit code $exitCode"
}

$exe = Join-Path $Destination 'dist\croco_wx.exe'
if (Test-Path -LiteralPath $exe) {
    Write-Step "Build succeeded: $exe"
    Copy-Item -LiteralPath $exe -Destination (New-Item -ItemType Directory -Force -Path (Join-Path $SourceRoot 'dist')).FullName -Force
}
else {
    Write-Warning "Build reported success but $exe was not found."
}
