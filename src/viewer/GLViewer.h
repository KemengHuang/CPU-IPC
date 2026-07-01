#pragma once
#ifndef GL_VIEWER_H
#define GL_VIEWER_H

#include "core/Simulator.h"
#include <string>

// Global simulator instance used by the GLUT viewer and headless runner.
extern FEMSimulator simulator;

// Global toggle for saving surface meshes each step.
extern bool saveSurface;

void initGL();
void initScene();

// Run the GLUT main loop. This function initializes GLUT, creates the window,
// wires callbacks, and never returns.
void runGLUTMainLoop(int argc, char** argv);

// Write the current tetrahedral surface mesh to an OBJ file under path.
void saveSurfaceMesh(const std::string& path);

#endif // GL_VIEWER_H
