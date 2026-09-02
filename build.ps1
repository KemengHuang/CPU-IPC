param(
    [string]$VcpkgRoot = "",
    [string]$BuildDirectory = "",
    [ValidateSet("Release", "RelWithDebInfo", "Debug")]
    [string]$Configuration = "Release",
    [switch]$HeadlessOnly
)

$ErrorActionPreference = "Stop"
$repositoryRoot = $PSScriptRoot
. (Join-Path $repositoryRoot "scripts/ResolveVcpkg.ps1")
$vcpkg = Resolve-VcpkgInstallation `
    -RequestedRoot $VcpkgRoot `
    -BootstrapRoot (Join-Path $repositoryRoot "build/_deps/vcpkg")
$VcpkgRoot = $vcpkg.Root
if ([string]::IsNullOrWhiteSpace($BuildDirectory)) {
    $BuildDirectory = Join-Path $repositoryRoot "build"
}

$vcpkgExecutable = $vcpkg.Executable
$triplet = "x64-windows"
if (-not $HeadlessOnly) {
    & $vcpkgExecutable install "freeglut:$triplet"
    if ($LASTEXITCODE -ne 0) {
        throw "Viewer dependency installation failed."
    }
}
& $vcpkgExecutable install `
    "eigen3:$triplet" `
    "tbb:$triplet" `
    "metis:$triplet"
if ($LASTEXITCODE -ne 0) {
    throw "Base dependency installation failed."
}

$cholmodInstall = Join-Path $repositoryRoot "build/cholmod-mkl-install"
& (Join-Path $repositoryRoot "scripts/build_cholmod_mkl.ps1") `
    -VcpkgRoot $VcpkgRoot `
    -InstallPrefix $cholmodInstall
if ($LASTEXITCODE -ne 0) {
    throw "Optimized CHOLMOD build failed."
}

$configureArguments = @(
    "-S", $repositoryRoot,
    "-B", $BuildDirectory,
    "-DCMAKE_TOOLCHAIN_FILE=$($vcpkg.Toolchain)",
    "-DCIPC_CHOLMOD_ROOT=$cholmodInstall",
    "-DCIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON"
)
if ($HeadlessOnly) {
    $configureArguments += "-DCIPC_BUILD_VIEWER=OFF"
}

& cmake @configureArguments
if ($LASTEXITCODE -ne 0) {
    throw "CPU-IPC configuration failed."
}
& cmake --build $BuildDirectory --config $Configuration --parallel
if ($LASTEXITCODE -ne 0) {
    throw "CPU-IPC build failed."
}

Write-Host "CPU-IPC $Configuration build completed: $BuildDirectory/$Configuration"
