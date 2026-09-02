param(
    [string]$VcpkgRoot = "",
    [string]$BuildDirectory = "",
    [string]$InstallPrefix = "",
    [switch]$SkipDependencyInstall
)

$ErrorActionPreference = "Stop"
$repositoryRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot "ResolveBuildTools.ps1")
. (Join-Path $PSScriptRoot "ResolveVcpkg.ps1")
$vcpkgRevision = (Get-Content -LiteralPath (Join-Path $PSScriptRoot "vcpkg-revision.txt") -Raw).Trim()
if ($vcpkgRevision -notmatch "^[0-9a-fA-F]{40}$") {
    throw "vcpkg-revision.txt must contain one 40-character Git revision."
}
$vcpkgBootstrapRoot = Get-CpuIpcVcpkgBootstrapRoot `
    -RepositoryRoot $repositoryRoot `
    -Revision $vcpkgRevision
Assert-WindowsCppToolchain
$vcpkg = Resolve-VcpkgInstallation `
    -RequestedRoot $VcpkgRoot `
    -BootstrapRoot $vcpkgBootstrapRoot `
    -Revision $vcpkgRevision
$VcpkgRoot = $vcpkg.Root
if ([string]::IsNullOrWhiteSpace($BuildDirectory)) {
    $BuildDirectory = Join-Path $repositoryRoot "build/suitesparse-mkl-bundle"
}
if ([string]::IsNullOrWhiteSpace($InstallPrefix)) {
    $InstallPrefix = Join-Path $repositoryRoot "build/cholmod-mkl-install"
}

$vcpkgExecutable = $vcpkg.Executable

if (-not $SkipDependencyInstall) {
    & $vcpkgExecutable install `
        "intel-mkl:x64-windows" `
        "tbb[core]:x64-windows" `
        "metis:x64-windows" `
        --recurse `
        --no-print-usage
    if ($LASTEXITCODE -ne 0) {
        throw "vcpkg dependency installation failed."
    }
}
$cmakeExecutable = Resolve-CMakeExecutable -VcpkgRoot $VcpkgRoot

$suiteSparseVersion = (Get-Content `
    -LiteralPath (Join-Path $PSScriptRoot "suitesparse-version.txt") `
    -Raw).Trim()
if ($suiteSparseVersion -notmatch "^[0-9]+\.[0-9]+\.[0-9]+$") {
    throw "suitesparse-version.txt must contain a semantic version."
}
$suiteSparseSourceParent = Join-Path $repositoryRoot "build/_deps"
$suiteSparseSource = Join-Path $suiteSparseSourceParent "SuiteSparse-$suiteSparseVersion"
$downloadDirectory = Join-Path $VcpkgRoot "downloads"
$fetchScript = Join-Path $repositoryRoot "cmake/FetchSuiteSparse.cmake"
$providerFile = Join-Path $repositoryRoot "cmake/CholmodMKLProvider.cmake"
$blasShim = Join-Path $repositoryRoot "cmake/SuiteSparseBLAS.cmake"
$installedPrefix = Join-Path $VcpkgRoot "installed/x64-windows"
# Increment recipe when the SuiteSparse configure arguments below change.
$bundleSignature = @(
    "schema=1",
    "recipe=1",
    "suitesparse=$suiteSparseVersion",
    "vcpkg_revision=$vcpkgRevision",
    "vcpkg=$([System.IO.Path]::GetFullPath($VcpkgRoot))",
    "vcpkg_executable=$((Get-FileHash -Algorithm SHA256 -LiteralPath $vcpkgExecutable).Hash)",
    "mkl_port=$((Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $VcpkgRoot 'ports/intel-mkl/vcpkg.json')).Hash)",
    "tbb_port=$((Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $VcpkgRoot 'ports/tbb/vcpkg.json')).Hash)",
    "metis_port=$((Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $VcpkgRoot 'ports/metis/vcpkg.json')).Hash)",
    "provider=$((Get-FileHash -Algorithm SHA256 -LiteralPath $providerFile).Hash)",
    "blas_shim=$((Get-FileHash -Algorithm SHA256 -LiteralPath $blasShim).Hash)",
    "fetch=$((Get-FileHash -Algorithm SHA256 -LiteralPath $fetchScript).Hash)"
) -join "`n"
$bundleStamp = Join-Path $InstallPrefix ".cpu-ipc-suite-sparse-bundle"
$requiredBundleFiles = @(
    "bin/suitesparseconfig.dll",
    "bin/amd.dll",
    "bin/camd.dll",
    "bin/ccolamd.dll",
    "bin/colamd.dll",
    "bin/cholmod.dll"
)
$bundleReady = Test-Path -LiteralPath $bundleStamp -PathType Leaf
foreach ($relativePath in $requiredBundleFiles) {
    if (-not (Test-Path -LiteralPath (Join-Path $InstallPrefix $relativePath) -PathType Leaf)) {
        $bundleReady = $false
    }
}
if ($bundleReady) {
    $storedSignature = (Get-Content -LiteralPath $bundleStamp -Raw).TrimEnd()
    $bundleReady = $storedSignature -eq $bundleSignature
}

if ($bundleReady) {
    Write-Host "[3/4] Reusing the optimized SuiteSparse/CHOLMOD bundle"
}
else {
    & $cmakeExecutable `
        "-DCPU_IPC_SUITESPARSE_SOURCE_PARENT=$suiteSparseSourceParent" `
        "-DCPU_IPC_DOWNLOAD_DIR=$downloadDirectory" `
        -P $fetchScript
    if ($LASTEXITCODE -ne 0) {
        throw "SuiteSparse source preparation failed."
    }

    $configureArguments = @(
        "-S", $suiteSparseSource,
        "-B", $BuildDirectory,
        "-DCMAKE_TOOLCHAIN_FILE=$($vcpkg.Toolchain)",
        "-DVCPKG_TARGET_TRIPLET=x64-windows",
        "-DCMAKE_PREFIX_PATH=$installedPrefix",
        "-DCMAKE_PROJECT_INCLUDE=$providerFile",
        "-DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON",
        "-DCMAKE_INSTALL_PREFIX=$InstallPrefix",
        "-DSUITESPARSE_ENABLE_PROJECTS=cholmod",
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
        "-DSUITESPARSE_USE_OPENMP=OFF",
        "-DSUITESPARSE_USE_FORTRAN=OFF",
        "-DSUITESPARSE_DEMOS=OFF",
        "-DSUITESPARSE_USE_STRICT=ON"
    )

    Write-Host "[3/4] Building the optimized SuiteSparse/CHOLMOD bundle with oneMKL"
    & $cmakeExecutable @configureArguments
    if ($LASTEXITCODE -ne 0) {
        throw "SuiteSparse/oneMKL configuration failed."
    }
    & $cmakeExecutable --build $BuildDirectory --config Release --parallel
    if ($LASTEXITCODE -ne 0) {
        throw "SuiteSparse/oneMKL build failed."
    }
    & $cmakeExecutable --install $BuildDirectory --config Release
    if ($LASTEXITCODE -ne 0) {
        throw "SuiteSparse/oneMKL installation failed."
    }
    [System.IO.File]::WriteAllText(
        $bundleStamp,
        $bundleSignature + [Environment]::NewLine,
        [System.Text.UTF8Encoding]::new($false))
}

Write-Host "Optimized CHOLMOD installed to: $InstallPrefix"
Write-Host "Configure CPU-IPC with: -DCIPC_CHOLMOD_ROOT=$InstallPrefix"
