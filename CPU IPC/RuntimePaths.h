#pragma once

#include <string>

namespace RuntimePaths {

// Must be called before initialize() or any path accessor is used.
// An empty directory is rejected.
bool setOutputDirectory(const std::string& directory);
const std::string& outputDirectory();

// Creates all runtime output directories. Safe to call repeatedly.
bool initialize();

std::string outputFile(const std::string& filename);
std::string tempFile(const std::string& filename);
std::string surfaceFile(const std::string& filename);
std::string screenshotFile(const std::string& filename);

} // namespace RuntimePaths
