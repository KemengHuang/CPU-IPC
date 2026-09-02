function Resolve-VcpkgInstallation {
    param(
        [string]$RequestedRoot = "",
        [string]$BootstrapRoot = ""
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
        foreach ($name in @("vcpkg.exe", "vcpkg")) {
            $command = Get-Command $name -ErrorAction SilentlyContinue
            if ($null -ne $command) {
                $candidateExecutable = $command.Source
                $candidateRoot = Split-Path -Parent $candidateExecutable
                break
            }
        }
    }

    if ($null -eq $candidateExecutable) {
        if ([string]::IsNullOrWhiteSpace($BootstrapRoot)) {
            throw "vcpkg was not found. Pass -VcpkgRoot, set VCPKG_ROOT, or add vcpkg to PATH."
        }
        $candidateRoot = [System.IO.Path]::GetFullPath($BootstrapRoot)
        $candidateExecutable = Join-Path $candidateRoot "vcpkg.exe"
        if (-not (Test-Path -LiteralPath $candidateExecutable -PathType Leaf)) {
            if (-not (Test-Path -LiteralPath (Join-Path $candidateRoot ".git"))) {
                New-Item -ItemType Directory -Force -Path (Split-Path -Parent $candidateRoot) | Out-Null
                & git clone --depth 1 https://github.com/microsoft/vcpkg.git $candidateRoot
                if ($LASTEXITCODE -ne 0) {
                    throw "Failed to clone vcpkg into $candidateRoot"
                }
            }
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
