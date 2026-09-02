param(
    [string]$VcpkgRoot = "",
    [string]$BuildDirectory = "",
    [string]$InstallPrefix = ""
)

$ErrorActionPreference = "Stop"
$repositoryRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot "ResolveVcpkg.ps1")
$vcpkg = Resolve-VcpkgInstallation `
    -RequestedRoot $VcpkgRoot `
    -BootstrapRoot (Join-Path $repositoryRoot "build/_deps/vcpkg")
$VcpkgRoot = $vcpkg.Root
if ([string]::IsNullOrWhiteSpace($BuildDirectory)) {
    $BuildDirectory = Join-Path $repositoryRoot "build/cholmod-mkl-lib"
}
if ([string]::IsNullOrWhiteSpace($InstallPrefix)) {
    $InstallPrefix = Join-Path $repositoryRoot "build/cholmod-mkl-install"
}

$vcpkgExecutable = $vcpkg.Executable

& $vcpkgExecutable install `
    "intel-mkl:x64-windows" `
    "suitesparse-cholmod[partition]:x64-windows" `
    --recurse
if ($LASTEXITCODE -ne 0) {
    throw "vcpkg dependency installation failed."
}

$sourceParent = Join-Path $VcpkgRoot "buildtrees/suitesparse-cholmod/src"
$sourceRoot = Get-ChildItem -LiteralPath $sourceParent -Directory |
    Where-Object { $_.Name.EndsWith(".clean") } |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1
if ($null -eq $sourceRoot) {
    throw "The extracted SuiteSparse/CHOLMOD source tree was not found."
}
$cholmodSource = Join-Path $sourceRoot.FullName "CHOLMOD"
$providerFile = Join-Path $repositoryRoot "cmake/CholmodMKLProvider.cmake"
$installedPrefix = Join-Path $VcpkgRoot "installed/x64-windows"

$configureArguments = @(
    "-S", $cholmodSource,
    "-B", $BuildDirectory,
    "-DCMAKE_PREFIX_PATH=$installedPrefix",
    "-DCMAKE_PROJECT_INCLUDE=$providerFile",
    "-DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON",
    "-DSuiteSparse_config_DIR=$installedPrefix/share/SuiteSparse_config",
    "-DAMD_DIR=$installedPrefix/share/AMD",
    "-DCOLAMD_DIR=$installedPrefix/share/COLAMD",
    "-DCAMD_DIR=$installedPrefix/share/CAMD",
    "-DCCOLAMD_DIR=$installedPrefix/share/CCOLAMD",
    "-DCMAKE_INSTALL_PREFIX=$InstallPrefix",
    "-DBUILD_SHARED_LIBS=ON",
    "-DBUILD_STATIC_LIBS=OFF",
    "-DBUILD_TESTING=OFF",
    "-DCHOLMOD_GPL=ON",
    "-DCHOLMOD_SUPERNODAL=ON",
    "-DCHOLMOD_PARTITION=ON",
    "-DCHOLMOD_MATRIXOPS=OFF",
    "-DCHOLMOD_MODIFY=OFF",
    "-DCHOLMOD_USE_OPENMP=OFF",
    "-DCHOLMOD_USE_CUDA=OFF",
    "-DSUITESPARSE_USE_CUDA=OFF",
    "-DSUITESPARSE_DEMOS=OFF"
)

& cmake @configureArguments
if ($LASTEXITCODE -ne 0) {
    throw "CHOLMOD/oneMKL configuration failed."
}
& cmake --build $BuildDirectory --config Release --parallel
if ($LASTEXITCODE -ne 0) {
    throw "CHOLMOD/oneMKL build failed."
}
& cmake --install $BuildDirectory --config Release
if ($LASTEXITCODE -ne 0) {
    throw "CHOLMOD/oneMKL installation failed."
}

Write-Host "Optimized CHOLMOD installed to: $InstallPrefix"
Write-Host "Configure CPU-IPC with: -DCIPC_CHOLMOD_ROOT=$InstallPrefix"
