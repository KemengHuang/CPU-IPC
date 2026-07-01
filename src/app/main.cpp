#include "viewer/GLViewer.h"
#include <cstdlib>
#include <iostream>
#include <string>

int runHeadlessSimulation(GLViewer& viewer, int frameCount) {
    viewer.simulator().buildModels(0, 3);
    for (int step = 0; step < frameCount; ++step) {
        int k = viewer.simulator().simulateStick(step);
        if (viewer.saveSurfaceFlag()) {
            viewer.saveSurfaceMesh("saveSurface/surf_");
        }
        std::printf("current step:  %d\n", step);
    }
    return 0;
}

int main(int argc, char** argv) {
    int headlessFrames = 0;
    bool headless = false;
    for (int i = 1; i < argc; ++i) {
        if ((std::string(argv[i]) == "--headless" || std::string(argv[i]) == "-headless") && i + 1 < argc) {
            headless = true;
            headlessFrames = std::atoi(argv[i + 1]);
            ++i;
        }
    }

    GLViewer viewer;
    if (headless) {
        viewer.saveSurfaceFlag() = true;
        return runHeadlessSimulation(viewer, headlessFrames);
    }

    viewer.init(argc, argv);
    viewer.runMainLoop();
    return 0;
}
