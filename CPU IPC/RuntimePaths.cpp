#include "RuntimePaths.h"

#include <algorithm>
#include <cerrno>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#ifndef CIPC_OUTPUT_DIR
#define CIPC_OUTPUT_DIR "Output"
#endif

namespace {

std::string joinPath(const std::string& directory, const std::string& filename)
{
    if (directory.empty()) {
        return filename;
    }

    const char last = directory.back();
    if (last == '/' || last == '\\') {
        return directory + filename;
    }
    return directory + "/" + filename;
}

bool createSingleDirectory(const std::string& path)
{
#ifdef _WIN32
    const int result = _mkdir(path.c_str());
#else
    const int result = mkdir(path.c_str(), 0755);
#endif
    return result == 0 || errno == EEXIST;
}

std::string& mutableOutputDirectory()
{
    static std::string path = CIPC_OUTPUT_DIR;
    return path;
}

bool createDirectories(const std::string& inputPath)
{
    if (inputPath.empty()) {
        return false;
    }

    std::string path = inputPath;
    std::replace(path.begin(), path.end(), '\\', '/');
    while (path.size() > 1 && path.back() == '/') {
        path.pop_back();
    }

    size_t start = 0;
    if (path.size() >= 3 && path[1] == ':' && path[2] == '/') {
        start = 3;
    }
    else if (!path.empty() && path[0] == '/') {
        start = 1;
    }

    for (size_t slash = path.find('/', start); slash != std::string::npos;
         slash = path.find('/', slash + 1)) {
        const std::string prefix = path.substr(0, slash);
        if (!prefix.empty() && !createSingleDirectory(prefix)) {
            return false;
        }
    }
    return createSingleDirectory(path);
}

std::string tempDirectory()
{
    return joinPath(mutableOutputDirectory(), "tempData");
}

std::string surfaceDirectory()
{
    return joinPath(mutableOutputDirectory(), "saveSurface");
}

std::string screenshotDirectory()
{
    return joinPath(mutableOutputDirectory(), "saveScreen");
}

} // namespace

namespace RuntimePaths {

bool setOutputDirectory(const std::string& directory)
{
    if (directory.empty()) {
        return false;
    }
    mutableOutputDirectory() = directory;
    return true;
}

const std::string& outputDirectory()
{
    return mutableOutputDirectory();
}

bool initialize()
{
    return createDirectories(outputDirectory())
        && createDirectories(tempDirectory())
        && createDirectories(surfaceDirectory())
        && createDirectories(screenshotDirectory());
}

std::string outputFile(const std::string& filename)
{
    return joinPath(outputDirectory(), filename);
}

std::string tempFile(const std::string& filename)
{
    return joinPath(tempDirectory(), filename);
}

std::string surfaceFile(const std::string& filename)
{
    return joinPath(surfaceDirectory(), filename);
}

std::string screenshotFile(const std::string& filename)
{
    return joinPath(screenshotDirectory(), filename);
}

} // namespace RuntimePaths
