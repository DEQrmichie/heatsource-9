# This script installs hs.exe from this folder into the versioned folder at
# C:\Users\username\AppData\Local\Programs\HeatSource
# then copies the executable into C:\Users\username\AppData\Local\Programs\HeatSource\current
# then adds the 'current' folder to the user's PATH

param(
    [string]$Version = "9.0.1",
    [string]$SourceExePath = "",
    [string]$InstallRoot = ""
)

$ErrorActionPreference = "Stop"

function Resolve-NormalizedPath {
    param([string]$PathValue)

    if ([string]::IsNullOrWhiteSpace($PathValue)) {
        return ""
    }

    return [System.IO.Path]::GetFullPath($PathValue).TrimEnd("\")
}

if ([string]::IsNullOrWhiteSpace($SourceExePath)) {
    $localExe = Join-Path $PSScriptRoot "hs.exe"
    if (Test-Path -LiteralPath $localExe -PathType Leaf) {
        $SourceExePath = $localExe
    }
    else {
        $repoRoot = Split-Path -Parent $PSScriptRoot
        $SourceExePath = Join-Path $repoRoot "dist\_dev\hs.exe"
    }
}

if ([string]::IsNullOrWhiteSpace($InstallRoot)) {
    $InstallRoot = Join-Path $env:LOCALAPPDATA "Programs\HeatSource"
}

$SourceExePath = Resolve-NormalizedPath $SourceExePath
$InstallRoot = Resolve-NormalizedPath $InstallRoot

if (-not (Test-Path -LiteralPath $SourceExePath -PathType Leaf)) {
    throw "Source executable not found: $SourceExePath"
}

$versionDir = Join-Path $InstallRoot ("heatsource" + $Version)
$currentDir = Join-Path $InstallRoot "current"
$versionExePath = Join-Path $versionDir "hs.exe"
$currentExePath = Join-Path $currentDir "hs.exe"

New-Item -ItemType Directory -Force -Path $versionDir | Out-Null
New-Item -ItemType Directory -Force -Path $currentDir | Out-Null

Copy-Item -LiteralPath $SourceExePath -Destination $versionExePath -Force
Copy-Item -LiteralPath $versionExePath -Destination $currentExePath -Force

$userPath = [Environment]::GetEnvironmentVariable("Path", "User")
$pathEntries = @()

if (-not [string]::IsNullOrWhiteSpace($userPath)) {
    $pathEntries = @(
        $userPath.Split(";") |
        ForEach-Object { $_.Trim() } |
        Where-Object { $_ -ne "" }
    )
}

$normalizedCurrentDir = Resolve-NormalizedPath $currentDir
$pathExists = $false

foreach ($entry in $pathEntries) {
    if ((Resolve-NormalizedPath $entry) -ieq $normalizedCurrentDir) {
        $pathExists = $true
        break
    }
}

if (-not $pathExists) {
    $pathEntries += $currentDir
    $newUserPath = ($pathEntries -join ";")
    [Environment]::SetEnvironmentVariable("Path", $newUserPath, "User")

    if ([string]::IsNullOrWhiteSpace($env:Path)) {
        $env:Path = $currentDir
    }
    else {
        $env:Path = $currentDir + ";" + $env:Path
    }
}

Write-Host ""
Write-Host "Installed hs.exe to:"
Write-Host "  $versionExePath"
Write-Host ""
Write-Host "Added to user PATH:"
Write-Host "  $currentDir"
Write-Host ""
Write-Host "Open a new Command Prompt, PowerShell, or Windows Terminal window, then run:"
Write-Host "  hs -v"
Write-Host ""
Write-Host "If it prints the model version number, it worked."