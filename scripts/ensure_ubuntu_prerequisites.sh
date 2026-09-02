#!/usr/bin/env bash

# Sourced by build.sh. Keep unavoidable OS-level tooling here; all numerical
# libraries are installed separately through vcpkg.
ensure_ubuntu_prerequisites() {
    local cache_root="$1"
    local build_viewer="$2"
    local install_system_packages="$3"

    if [[ "$(uname -m)" != "x86_64" ]]; then
        echo "CPU-IPC's oneMKL build currently supports x86_64 Ubuntu/WSL only." >&2
        return 1
    fi

    # Reuse the no-sudo zip/unzip cache created by an earlier invocation.
    local user_tools="${cache_root}/tools"
    if [[ -d "${user_tools}/usr/bin" ]]; then
        export PATH="${user_tools}/usr/bin:${PATH}"
    fi

    # vcpkg needs zip/unzip very early. Ubuntu can download .deb archives as a
    # normal user, so these two small tools do not justify a sudo prompt.
    if ! command -v zip >/dev/null 2>&1 || ! command -v unzip >/dev/null 2>&1; then
        if command -v apt-get >/dev/null 2>&1 && command -v dpkg-deb >/dev/null 2>&1; then
            local package_cache="${cache_root}/apt-packages"
            mkdir -p "$package_cache" "$user_tools"
            (
                cd "$package_cache"
                apt-get download zip unzip
            )
            local package
            for package in "${package_cache}"/*.deb; do
                dpkg-deb -x "$package" "$user_tools"
            done
            export PATH="${user_tools}/usr/bin:${PATH}"
        fi
    fi

    local -a missing_packages=()
    command -v g++ >/dev/null 2>&1 || missing_packages+=(build-essential)
    command -v make >/dev/null 2>&1 || missing_packages+=(build-essential)
    command -v ninja >/dev/null 2>&1 || missing_packages+=(ninja-build)
    command -v git >/dev/null 2>&1 || missing_packages+=(git)
    command -v curl >/dev/null 2>&1 || missing_packages+=(curl)
    command -v tar >/dev/null 2>&1 || missing_packages+=(tar)
    command -v pkg-config >/dev/null 2>&1 || missing_packages+=(pkg-config)
    command -v zip >/dev/null 2>&1 || missing_packages+=(zip)
    command -v unzip >/dev/null 2>&1 || missing_packages+=(unzip)

    # FreeGLUT needs the core native X11/OpenGL development interfaces. Its
    # optional Xi/Xrandr/Xxf86vm paths can be disabled by the port, so checking
    # package names directly would request sudo on systems that already build
    # and run the viewer successfully.
    if [[ "$build_viewer" == "ON" ]] \
        && ! pkg-config --exists x11 gl glu 2>/dev/null; then
        missing_packages+=(libx11-dev libgl1-mesa-dev libglu1-mesa-dev)
    fi

    if [[ ${#missing_packages[@]} -ne 0 ]]; then
        # Preserve first occurrence while removing duplicates.
        local -a unique_packages=()
        local candidate existing seen
        for candidate in "${missing_packages[@]}"; do
            seen="OFF"
            for existing in "${unique_packages[@]}"; do
                if [[ "$candidate" == "$existing" ]]; then
                    seen="ON"
                    break
                fi
            done
            [[ "$seen" == "ON" ]] || unique_packages+=("$candidate")
        done

        if [[ "$install_system_packages" != "ON" ]]; then
            echo "Missing Ubuntu system packages: ${unique_packages[*]}" >&2
            echo "Install them with:" >&2
            echo "  sudo apt-get update && sudo apt-get install -y ${unique_packages[*]}" >&2
            return 1
        fi
        if ! command -v apt-get >/dev/null 2>&1; then
            echo "Missing build tools and apt-get is unavailable: ${unique_packages[*]}" >&2
            return 1
        fi

        local -a apt_command=(apt-get)
        if [[ "${EUID}" -ne 0 ]]; then
            if ! command -v sudo >/dev/null 2>&1; then
                echo "sudo is required once to install: ${unique_packages[*]}" >&2
                return 1
            fi
            apt_command=(sudo apt-get)
        fi

        echo "Installing the missing Ubuntu build prerequisites: ${unique_packages[*]}"
        "${apt_command[@]}" update
        "${apt_command[@]}" install -y "${unique_packages[@]}"
    fi

    local -a required_commands=(g++ make ninja git curl tar pkg-config zip unzip)
    local -a unresolved_commands=()
    local command_name
    for command_name in "${required_commands[@]}"; do
        command -v "$command_name" >/dev/null 2>&1 || unresolved_commands+=("$command_name")
    done
    if [[ ${#unresolved_commands[@]} -ne 0 ]]; then
        echo "Ubuntu prerequisite setup did not provide: ${unresolved_commands[*]}" >&2
        return 1
    fi
}

resolve_cmake_executable() {
    local vcpkg_root="$1"
    local -a candidates=()
    if command -v cmake >/dev/null 2>&1; then
        candidates+=("$(command -v cmake)")
    fi

    local candidate
    while IFS= read -r candidate; do
        candidates+=("$candidate")
    done < <(find "${vcpkg_root}/downloads/tools" -type f -name cmake \
        -path '*/bin/cmake' -perm -u+x 2>/dev/null | sort -r)

    local version minimum_version="3.23"
    for candidate in "${candidates[@]}"; do
        version="$("$candidate" --version 2>/dev/null | awk 'NR == 1 { print $3 }')"
        if [[ -n "$version" \
            && "$(printf '%s\n' "$minimum_version" "$version" | sort -V | head -n 1)" == "$minimum_version" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done

    echo "CMake 3.23 or newer was not found after vcpkg dependency setup." >&2
    return 1
}
