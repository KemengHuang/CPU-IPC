param(
    [string]$VcpkgRoot = "",
    [string]$BuildDirectory = "",
    [ValidateSet("Release", "RelWithDebInfo", "Debug")]
    [string]$Configuration = "Release",
    [switch]$HeadlessOnly,
    [switch]$NonQuadraticBending,
    [switch]$DependenciesOnly,
    [switch]$Help
)

if ($Help) {
    Write-Host @"
Usage: .\build.cmd [options]
  -HeadlessOnly       Skip FreeGLUT and the viewer.
  -NonQuadraticBending
                      Build the complete dihedral-hinge bending model.
  -DependenciesOnly   Install/build dependencies without building CPU-IPC.
  -Configuration TYPE Release, RelWithDebInfo, or Debug.
  -BuildDirectory DIR Override the default .\build\cpu-ipc directory.
  -VcpkgRoot PATH     Use an explicitly managed vcpkg checkout.
"@
    exit 0
}

$ErrorActionPreference = "Stop"
$repositoryRoot = $PSScriptRoot
. (Join-Path $repositoryRoot "scripts/ResolveBuildTools.ps1")
. (Join-Path $repositoryRoot "scripts/ResolveVcpkg.ps1")
$vcpkgRevision = (Get-Content -LiteralPath (Join-Path $repositoryRoot "scripts/vcpkg-revision.txt") -Raw).Trim()
if ($vcpkgRevision -notmatch "^[0-9a-fA-F]{40}$") {
    throw "scripts/vcpkg-revision.txt must contain one 40-character Git revision."
}
$vcpkgBootstrapRoot = Get-CpuIpcVcpkgBootstrapRoot `
    -RepositoryRoot $repositoryRoot `
    -Revision $vcpkgRevision
Write-Host "[1/4] Checking the Windows C++ build environment"
Assert-WindowsCppToolchain

Write-Host "[2/4] Preparing vcpkg and numerical dependencies"
$vcpkg = Resolve-VcpkgInstallation `
    -RequestedRoot $VcpkgRoot `
    -BootstrapRoot $vcpkgBootstrapRoot `
    -Revision $vcpkgRevision
$VcpkgRoot = $vcpkg.Root
if ([string]::IsNullOrWhiteSpace($BuildDirectory)) {
    $BuildDirectory = Join-Path $repositoryRoot "build/cpu-ipc"
}

$vcpkgExecutable = $vcpkg.Executable
$triplet = "x64-windows"
$dependencies = @(
    "eigen3:$triplet",
    "tbb[core]:$triplet",
    "metis:$triplet",
    "intel-mkl:$triplet"
)
if (-not $HeadlessOnly) {
    $dependencies += "freeglut:$triplet"
}
& $vcpkgExecutable install @dependencies --recurse --no-print-usage
if ($LASTEXITCODE -ne 0) {
    throw "Dependency installation failed."
}

$cholmodInstall = Join-Path $repositoryRoot "build/cholmod-mkl-install"
& (Join-Path $repositoryRoot "scripts/build_cholmod_mkl.ps1") `
    -VcpkgRoot $VcpkgRoot `
    -InstallPrefix $cholmodInstall `
    -SkipDependencyInstall
if ($LASTEXITCODE -ne 0) {
    throw "Optimized CHOLMOD build failed."
}
if ($DependenciesOnly) {
    Write-Host "CPU-IPC dependencies are ready."
    Write-Host "vcpkg: $VcpkgRoot"
    Write-Host "Optimized CHOLMOD: $cholmodInstall"
    exit 0
}

$cmakeExecutable = Resolve-CMakeExecutable -VcpkgRoot $VcpkgRoot
Write-Host "[4/4] Configuring and building CPU-IPC"
$configureArguments = @(
    "-S", $repositoryRoot,
    "-B", $BuildDirectory,
    "-DCIPC_BUILD_VIEWER=$(if ($HeadlessOnly) { 'OFF' } else { 'ON' })",
    "-DCIPC_ENABLE_FRICTION=ON",
    "-DCIPC_ENABLE_QUADRATIC_BENDING=$(if ($NonQuadraticBending) { 'OFF' } else { 'ON' })",
    "-DCIPC_ENABLE_METIS_ORDERING=ON",
    "-DCIPC_ENABLE_PARDISO=ON",
    "-DCIPC_CHOLMOD_ROOT=$cholmodInstall",
    "-DCIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON"
)
$cmakeCache = Join-Path $BuildDirectory "CMakeCache.txt"
if (Test-Path -LiteralPath $cmakeCache -PathType Leaf) {
    $toolchainEntry = Select-String `
        -LiteralPath $cmakeCache `
        -Pattern "^CMAKE_TOOLCHAIN_FILE:[^=]+=" |
        Select-Object -First 1
    if ($null -eq $toolchainEntry) {
        throw "The existing build directory was not configured with vcpkg. Choose another -BuildDirectory."
    }
    $cachedToolchain = ($toolchainEntry.Line -split "=", 2)[1].Replace("/", "\")
    $requestedToolchain = $vcpkg.Toolchain.Replace("/", "\")
    if (-not $cachedToolchain.Equals($requestedToolchain, [StringComparison]::OrdinalIgnoreCase)) {
        throw "The existing build directory uses another vcpkg checkout. Choose another -BuildDirectory."
    }
}
else {
    $configureArguments += "-DCMAKE_TOOLCHAIN_FILE=$($vcpkg.Toolchain)"
}
& $cmakeExecutable @configureArguments
if ($LASTEXITCODE -ne 0) {
    throw "CPU-IPC configuration failed."
}
& $cmakeExecutable --build $BuildDirectory --config $Configuration --parallel
if ($LASTEXITCODE -ne 0) {
    throw "CPU-IPC build failed."
}

Write-Host "CPU-IPC $Configuration build completed: $BuildDirectory/$Configuration"
