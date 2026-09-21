<#
.SYNOPSIS
    Installs the build tools and builds the cccrpc~ Max/MSP external on Windows.

.DESCRIPTION
    On a machine with no developer tools this installs (via winget, which ships
    with Windows 10/11):
      - Git
      - CMake
      - Visual Studio 2022 Build Tools with the "Desktop development with C++"
        workload (MSVC + Windows SDK; several GB, takes a while)
    then clones max-sdk-base next to this repository, and builds cccrpc~.mxe64.

    Can be run from inside a libcccrt checkout (plugins\max\build-windows.ps1)
    or on its own, in which case it clones libcccrt into %USERPROFILE%\src.

.USAGE
    Open PowerShell *as Administrator* (the Build Tools installer needs it) and run:

        powershell -ExecutionPolicy Bypass -File build-windows.ps1

    or, without a checkout:

        Invoke-WebRequest https://raw.githubusercontent.com/chriskiefer/libcccrt/eigen-free-core/plugins/max/build-windows.ps1 -OutFile build-windows.ps1
        powershell -ExecutionPolicy Bypass -File build-windows.ps1

.PARAMETER Install
    After building, copy the external and help file into Max's Library folder
    (Documents\Max 9\Library or Documents\Max 8\Library, whichever exists).

.PARAMETER SkipToolInstall
    Assume git, cmake and MSVC are already installed and on PATH.

.PARAMETER MaxSdkBase
    Path to an existing max-sdk-base checkout. Default: a sibling of the repo.

.PARAMETER Branch
    libcccrt branch to clone when not run from inside a checkout.
#>
[CmdletBinding()]
param(
    [switch]$Install,
    [switch]$SkipToolInstall,
    [string]$MaxSdkBase = "",
    [string]$RepoUrl = "https://github.com/chriskiefer/libcccrt.git",
    [string]$Branch = "eigen-free-core"
)

$ErrorActionPreference = "Stop"

function Write-Step($msg) { Write-Host "`n==> $msg" -ForegroundColor Cyan }

# winget installs update the registry PATH, not the PATH of this process.
function Refresh-Path {
    $env:Path = [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" +
                [Environment]::GetEnvironmentVariable("Path", "User")
}

function Test-Command($name) {
    return $null -ne (Get-Command $name -ErrorAction SilentlyContinue)
}

function Install-WithWinget($id, $overrideArgs) {
    $wingetArgs = @("install", "--id", $id, "-e", "--accept-source-agreements", "--accept-package-agreements")
    if ($overrideArgs) { $wingetArgs += @("--override", $overrideArgs) }
    & winget @wingetArgs
    # winget returns non-zero when the package is already installed; that's fine
    Refresh-Path
}

function Get-VsInstallPath {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
    if (-not (Test-Path $vswhere)) { return $null }
    $path = & $vswhere -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath -latest
    if ($path) { return $path.Trim() }
    return $null
}

# ---------------------------------------------------------------- tools

if (-not $SkipToolInstall) {
    if (-not (Test-Command winget)) {
        throw "winget not found. Install 'App Installer' from the Microsoft Store, or install Git, CMake and Visual Studio Build Tools manually and rerun with -SkipToolInstall."
    }

    $isAdmin = ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
    if (-not $isAdmin) {
        Write-Warning "Not running as Administrator. The Visual Studio Build Tools install will prompt for elevation (or fail)."
    }

    Write-Step "Git"
    if (Test-Command git) { Write-Host "already installed: $(git --version)" }
    else { Install-WithWinget "Git.Git" }

    Write-Step "CMake"
    if (Test-Command cmake) { Write-Host "already installed: $((cmake --version)[0])" }
    else { Install-WithWinget "Kitware.CMake" }

    Write-Step "Visual Studio 2022 Build Tools (C++ workload)"
    $vs = Get-VsInstallPath
    if ($vs) {
        Write-Host "already installed: $vs"
    } else {
        Write-Host "Installing - this downloads several GB and can take 10-30 minutes..."
        Install-WithWinget "Microsoft.VisualStudio.2022.BuildTools" `
            "--passive --wait --norestart --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
        $vs = Get-VsInstallPath
        if (-not $vs) {
            throw "Visual Studio Build Tools installed but the C++ tools were not found. Open 'Visual Studio Installer', choose Modify, and add 'Desktop development with C++'. Then rerun this script."
        }
    }

    Refresh-Path
    foreach ($tool in @("git", "cmake")) {
        if (-not (Test-Command $tool)) { throw "$tool is installed but not on PATH. Open a new PowerShell window and rerun the script." }
    }
}

# ---------------------------------------------------------------- sources

Write-Step "libcccrt sources"
$repoRoot = $null
if ($PSScriptRoot -and (Test-Path (Join-Path $PSScriptRoot "..\..\core\rpc.hpp"))) {
    $repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
    Write-Host "using checkout at $repoRoot"
} else {
    $srcDir = Join-Path $env:USERPROFILE "src"
    $repoRoot = Join-Path $srcDir "libcccrt"
    if (Test-Path (Join-Path $repoRoot "core\rpc.hpp")) {
        Write-Host "using existing clone at $repoRoot"
    } else {
        New-Item -ItemType Directory -Force $srcDir | Out-Null
        Write-Host "cloning $RepoUrl ($Branch) into $repoRoot"
        git clone --branch $Branch $RepoUrl $repoRoot
    }
}

Write-Step "max-sdk-base"
if (-not $MaxSdkBase) { $MaxSdkBase = Join-Path (Split-Path $repoRoot -Parent) "max-sdk-base" }
if (Test-Path (Join-Path $MaxSdkBase "script\max-pretarget.cmake")) {
    Write-Host "using $MaxSdkBase"
} else {
    Write-Host "cloning into $MaxSdkBase"
    git clone --depth 1 https://github.com/Cycling74/max-sdk-base.git $MaxSdkBase
}

# ---------------------------------------------------------------- build

Write-Step "Building cccrpc~"
$externalDir = Join-Path $repoRoot "plugins\max\cccrpc~"
$buildDir = Join-Path $externalDir "build"
New-Item -ItemType Directory -Force $buildDir | Out-Null

cmake -S $externalDir -B $buildDir -G "Visual Studio 17 2022" -A x64 "-DMAX_SDK_BASE_PATH=$MaxSdkBase"
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }
cmake --build $buildDir --config Release
if ($LASTEXITCODE -ne 0) { throw "Build failed" }

$outDir = Join-Path $repoRoot "plugins\max\externals"
$external = Join-Path $outDir "cccrpc~.mxe64"
if (-not (Test-Path $external)) { throw "Build reported success but $external was not produced" }
Write-Host "`nBuilt: $external" -ForegroundColor Green

# ---------------------------------------------------------------- install

if ($Install) {
    Write-Step "Installing into Max Library"
    $docs = [Environment]::GetFolderPath("MyDocuments")
    $lib = $null
    foreach ($v in @("Max 9", "Max 8")) {
        $candidate = Join-Path $docs "$v\Library"
        if (Test-Path (Split-Path $candidate -Parent)) { $lib = $candidate; break }
    }
    if (-not $lib) {
        Write-Warning "No 'Documents\Max 9' or 'Documents\Max 8' folder found - is Max installed? Copy $external into a folder in Max's search path manually."
    } else {
        New-Item -ItemType Directory -Force $lib | Out-Null
        Copy-Item $external $lib -Force
        Copy-Item (Join-Path $externalDir "cccrpc~.maxhelp") $lib -Force
        Write-Host "copied cccrpc~.mxe64 and cccrpc~.maxhelp to $lib" -ForegroundColor Green
        Write-Host "Restart Max, then create a [cccrpc~] object or open its help (alt+click it)."
    }
} else {
    Write-Host "Copy it (and cccrpc~\cccrpc~.maxhelp) into Documents\Max 9\Library, or rerun with -Install."
}
