function Assert-WindowsCppToolchain {
    if (-not [Environment]::Is64BitOperatingSystem) {
        throw "CPU-IPC's oneMKL build requires 64-bit Windows."
    }

    if ($null -ne (Get-Command cl.exe -ErrorAction SilentlyContinue)) {
        return
    }

    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio/Installer/vswhere.exe"
    if (Test-Path -LiteralPath $vswhere -PathType Leaf) {
        $installation = & $vswhere `
            -latest `
            -products * `
            -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
            -property installationPath
        if (-not [string]::IsNullOrWhiteSpace($installation)) {
            return
        }
    }

    throw @"
The Microsoft C++ build tools were not found. Install the "Desktop development
with C++" workload once, then rerun build.cmd. With winget:
  winget install --id Microsoft.VisualStudio.2022.BuildTools --override "--wait --passive --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
"@
}

function Resolve-CMakeExecutable {
    param([string]$VcpkgRoot = "")

    function Test-CMakeCandidate {
        param([string]$Path)

        if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
            return $false
        }
        $versionLine = & $Path --version | Select-Object -First 1
        if ($versionLine -notmatch "cmake version ([0-9]+\.[0-9]+(?:\.[0-9]+)?)") {
            return $false
        }
        [version]$version = $Matches[1]
        return $version -ge [version]"3.23"
    }

    $command = Get-Command cmake.exe -ErrorAction SilentlyContinue
    if ($null -ne $command -and (Test-CMakeCandidate -Path $command.Source)) {
        return $command.Source
    }

    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio/Installer/vswhere.exe"
    if (Test-Path -LiteralPath $vswhere -PathType Leaf) {
        $installations = & $vswhere -products * -property installationPath
        foreach ($installation in $installations) {
            $candidate = Join-Path $installation "Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
            if (Test-CMakeCandidate -Path $candidate) {
                return $candidate
            }
        }
    }

    if (-not [string]::IsNullOrWhiteSpace($VcpkgRoot)) {
        $toolsRoot = Join-Path $VcpkgRoot "downloads/tools"
        if (Test-Path -LiteralPath $toolsRoot -PathType Container) {
            $bundledCMake = Get-ChildItem -LiteralPath $toolsRoot -Filter cmake.exe -File -Recurse |
                Sort-Object LastWriteTime -Descending |
                Select-Object -First 1
            foreach ($candidate in $bundledCMake) {
                if (Test-CMakeCandidate -Path $candidate.FullName) {
                    return $candidate.FullName
                }
            }
        }
    }

    throw @"
CMake 3.23 or newer was not found on PATH, in Visual Studio, or in vcpkg's downloaded tools.
Install it once with `winget install --id Kitware.CMake`, then rerun build.cmd.
"@
}
