#include "viewer/GLViewer.h"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

int runHeadlessSimulation(int frameCount) {
    initScene();
    for (int step = 0; step < frameCount; ++step) {
        int k = simulator.simulateStick(step);
        if (saveSurface) {
            saveSurfaceMesh("saveSurface/surf_");
        }
        printf("current step:  %d\n", step);
    }
    return 0;
}

int main(int argc, char** argv)
{
    int headlessFrames = 0;
    bool headless = false;
    for (int i = 1; i < argc; ++i) {
        if ((string(argv[i]) == "--headless" || string(argv[i]) == "-headless") && i + 1 < argc) {
            headless = true;
            headlessFrames = std::atoi(argv[i + 1]);
            ++i;
        }
    }

    if (headless) {
        // In headless regression mode default to saving surfaces so the run
        // produces verifiable output. saveSurface can still be toggled if a
        // future flag is added.
        saveSurface = true;
        return runHeadlessSimulation(headlessFrames);
    }

    runGLUTMainLoop(argc, argv);
}
