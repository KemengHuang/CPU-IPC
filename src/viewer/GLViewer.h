#pragma once

#include "mesh/Mesh.h"
#include "core/Simulator.h"
#include "GL/glew.h"
#include "GL/freeglut.h"
#include <string>
#include <vector>

class GLViewer {
public:
    GLViewer();

    void init(int argc, char** argv);
    void runMainLoop();

    // Camera / interaction
    void onKeyboard(unsigned char key, int x, int y);
    void onMouse(int button, int state, int x, int y);
    void onMotion(int x, int y);
    void onReshape(int width, int height);
    void onDisplay();
    void onIdle();

    // Helpers
    void saveScreenshot(const std::string& file_name);
    void loadShaders();
    void saveSurfaceMesh(const std::string& path);

    // Accessors
    FEMSimulator& simulator() { return simulator_; }
    bool& stop() { return stop_; }
    bool& saveSurfaceFlag() { return saveSurface_; }

private:
    // Window / camera state
    float xRot_ = 0.0f, yRot_ = 0.0f;
    float xTrans_ = 0.0f, yTrans_ = 0.0f, zTrans_ = 0.0f;
    int ox_ = 0, oy_ = 0;
    int buttonState_ = 0;
    float xRotLength_ = 0.0f, yRotLength_ = 0.0f;
    float windowWidth_ = 1000.0f, windowHeight_ = 1000.0f;
    int sDimension_ = 3;

    // Runtime flags
    bool stop_ = true;
    bool screenshot_ = true;
    bool saveSurface_ = false;
    bool mouthOnly_ = false;
    bool isSetShader_ = false;
    bool isStop_ = false;

    // Unused legacy toggles (keep for compatibility with existing keyboard handler)
    bool drawFiber_ = false;
    bool drawSurface_ = false;
    bool drawMuscle_ = true;
    bool drawFat_ = true;
    bool drawSkull_ = true;
    bool drawJaw_ = true;
    bool drawMouth_ = true;

    // GL / shader resources
    GLuint pnVbo_ = 0;
    GLuint vao_ = 0;
    GLuint colorVbo_ = 0;
    GLuint normalVbo_ = 0;
    GLuint vertexShader_ = 0;
    GLuint fragmentShader_ = 0;
    GLuint shaderProgram_ = 0;

    // Simulation
    FEMSimulator simulator_;
    int surfNumId_ = 0;

    // Static helpers
    bool writeBitmapFile(int width, int height, const std::string& file_name, std::vector<unsigned char>& bitmapData);
};
