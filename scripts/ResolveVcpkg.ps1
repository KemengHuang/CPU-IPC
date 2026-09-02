function Get-CpuIpcVcpkgBootstrapRoot {
    param(
        [string]$RepositoryRoot,
        [string]$Revision
    )

    if (-not [string]::IsNullOrWhiteSpace($env:LOCALAPPDATA)) {
        $cacheRoot = Join-Path $env:LOCALAPPDATA "CPU-IPC"
    }
    else {
        $cacheRoot = Join-Path $RepositoryRoot "build/_deps"
    }
    $revisionTag = $Revision
    if ($revisionTag.Length -gt 12) {
        $revisionTag = $revisionTag.Substring(0, 12)
    }
    Join-Path $cacheRoot "vcpkg-$revisionTag"
}

function Resolve-VcpkgInstallation {
    param(
        [string]$RequestedRoot = "",
        [string]$BootstrapRoot = "",
        [string]$Revision = ""
    )

    $candidateRoot = $RequestedRoot
    $configuredSource = "-VcpkgRoot"
    if ([string]::IsNullOrWhiteSpace($candidateRoot)) {
        $candidateRoot = $env:VCPKG_ROOT
        $configuredSource = "VCPKG_ROOT"
    }

    $candidateExecutable = $null
    $hasConfiguredRoot = -not [string]::IsNullOrWhiteSpace($candidateRoot)
    if ($hasConfiguredRoot) {
        if (Test-Path -LiteralPath $candidateRoot -PathType Leaf) {
            $candidateExecutable = (Resolve-Path -LiteralPath $candidateRoot).Path
            $candidateRoot = Split-Path -Parent $candidateExecutable
        }
        else {
            foreach ($name in @("vcpkg.exe", "vcpkg")) {
                $path = [System.IO.Path]::Combine($candidateRoot, $name)
                if (Test-Path -LiteralPath $path -PathType Leaf) {
                    $candidateExecutable = (Resolve-Path -LiteralPath $path).Path
                    $candidateRoot = (Resolve-Path -LiteralPath $candidateRoot).Path
                    break
                }
            }
        }
        if ($null -eq $candidateExecutable) {
            throw "$configuredSource does not point to a valid vcpkg checkout or executable: $candidateRoot"
        }
    }

    if ($null -eq $candidateExecutable) {
        if ([string]::IsNullOrWhiteSpace($BootstrapRoot)) {
            throw "vcpkg was not found. Pass -VcpkgRoot or set VCPKG_ROOT."
        }
        $candidateRoot = [System.IO.Path]::GetFullPath($BootstrapRoot)
        $candidateExecutable = Join-Path $candidateRoot "vcpkg.exe"
        $git = $null
        $needsBootstrap = -not (Test-Path -LiteralPath $candidateExecutable -PathType Leaf)
        $gitDirectory = Join-Path $candidateRoot ".git"
        if (-not [string]::IsNullOrWhiteSpace($Revision)) {
            if (-not (Test-Path -LiteralPath $gitDirectory -PathType Container)) {
                $git = Get-Command git.exe -ErrorAction SilentlyContinue
                if ($null -eq $git) {
                    throw @"
Git is required to download the private vcpkg checkout. Install it once with
`winget install --id Git.Git`, then rerun build.cmd.
"@
                }
                New-Item -ItemType Directory -Force -Path (Split-Path -Parent $candidateRoot) | Out-Null
                New-Item -ItemType Directory -Force -Path $candidateRoot | Out-Null
                & $git.Source -C $candidateRoot init
                if ($LASTEXITCODE -ne 0) {
                    throw "Failed to initialize the private vcpkg checkout: $candidateRoot"
                }
                & $git.Source -C $candidateRoot remote add origin https://github.com/microsoft/vcpkg.git
                if ($LASTEXITCODE -ne 0) {
                    throw "Failed to configure the vcpkg remote: $candidateRoot"
                }
            }

            $currentRevision = ""
            $headFile = Join-Path $gitDirectory "HEAD"
            if (Test-Path -LiteralPath $headFile -PathType Leaf) {
                $headValue = (Get-Content -LiteralPath $headFile -Raw).Trim()
                if ($headValue -match "^[0-9a-fA-F]{40}$") {
                    $currentRevision = $headValue
                }
            }
            if ($currentRevision -ne $Revision) {
                if ($null -eq $git) {
                    $git = Get-Command git.exe -ErrorAction SilentlyContinue
                }
                if ($null -eq $git) {
                    throw @"
Git is required to update the private vcpkg checkout. Install it once with
`winget install --id Git.Git`, then rerun build.cmd.
"@
                }
                & $git.Source -C $candidateRoot fetch --depth 1 origin $Revision
                if ($LASTEXITCODE -ne 0) {
                    throw "Failed to download vcpkg revision $Revision"
                }
                & $git.Source -C $candidateRoot checkout --detach FETCH_HEAD
                if ($LASTEXITCODE -ne 0) {
                    throw "Failed to activate vcpkg revision $Revision"
                }
                $needsBootstrap = $true
            }
        }
        elseif (-not (Test-Path -LiteralPath $gitDirectory -PathType Container)) {
            $git = Get-Command git.exe -ErrorAction SilentlyContinue
            if ($null -eq $git) {
                throw "Git is required to download vcpkg. Install it with: winget install --id Git.Git"
            }
            New-Item -ItemType Directory -Force -Path (Split-Path -Parent $candidateRoot) | Out-Null
            & $git.Source clone --depth 1 https://github.com/microsoft/vcpkg.git $candidateRoot
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to clone vcpkg into $candidateRoot"
            }
            $needsBootstrap = $true
        }

        if ($needsBootstrap) {
            $bootstrap = Join-Path $candidateRoot "bootstrap-vcpkg.bat"
            if (-not (Test-Path -LiteralPath $bootstrap -PathType Leaf)) {
                throw "The local vcpkg checkout is incomplete: $candidateRoot"
            }
            & $bootstrap -disableMetrics
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to bootstrap vcpkg under $candidateRoot"
            }
        }
    }

    $toolchain = Join-Path $candidateRoot "scripts/buildsystems/vcpkg.cmake"
    $hasExecutable = Test-Path -LiteralPath $candidateExecutable -PathType Leaf
    $hasToolchain = Test-Path -LiteralPath $toolchain -PathType Leaf
    if (-not $hasExecutable -or -not $hasToolchain) {
        throw "Invalid vcpkg installation: $candidateRoot"
    }

    [PSCustomObject]@{
        Root = $candidateRoot
        Executable = $candidateExecutable
        Toolchain = $toolchain
    }
}
