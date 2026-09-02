param(
    [string]$VcpkgRoot = $env:VCPKG_ROOT,
    [string]$BuildDirectory = "",
    [ValidateSet("Release", "RelWithDebInfo", "Debug")]
    [string]$Configuration = "Release",
    [switch]$HeadlessOnly
)

$ErrorActionPreference = "Stop"
$repositoryRoot = $PSScriptRoot

if ([string]::IsNullOrWhiteSpace($VcpkgRoot)) {
    $defaultVcpkgRoot = "D:/VCPKG/vcpkg"
    if (Test-Path -LiteralPath "$defaultVcpkgRoot/vcpkg.exe") {
        $VcpkgRoot = $defaultVcpkgRoot
    }
    else {
        throw "Set VCPKG_ROOT or pass -VcpkgRoot."
    }
}
if ([string]::IsNullOrWhiteSpace($BuildDirectory)) {
    $BuildDirectory = Join-Path $repositoryRoot "build"
}

$vcpkgExecutable = Join-Path $VcpkgRoot "vcpkg.exe"
$triplet = "x64-windows"
& $vcpkgExecutable install `
    "eigen3:$triplet" `
    "freeglut:$triplet" `
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
    "-DCMAKE_TOOLCHAIN_FILE=$VcpkgRoot/scripts/buildsystems/vcpkg.cmake",
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
