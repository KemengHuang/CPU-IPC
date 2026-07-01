#include "viewer/GLViewer.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

namespace {

static GLViewer* g_viewer = nullptr;

#pragma pack(push, 1)
typedef struct
{
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} mBITMAPFILEHEADER;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct
{
    uint32_t biSize;
    int32_t  biWidth;
    int32_t  biHeight;
    uint16_t biPlanes;
    uint16_t biBitCount;
    uint32_t biCompression;
    uint32_t biSizeImage;
    int32_t  biXPelsPerMeter;
    int32_t  biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
} mBITMAPINFOHEADER;
#pragma pack(pop)

static void displayThunk() { if (g_viewer) g_viewer->onDisplay(); }
static void reshapeThunk(int w, int h) { if (g_viewer) g_viewer->onReshape(w, h); }
static void keyboardThunk(unsigned char key, int x, int y) { if (g_viewer) g_viewer->onKeyboard(key, x, y); }
static void mouseThunk(int button, int state, int x, int y) { if (g_viewer) g_viewer->onMouse(button, state, x, y); }
static void motionThunk(int x, int y) { if (g_viewer) g_viewer->onMotion(x, y); }
static void idleThunk() { if (g_viewer) g_viewer->onIdle(); }

void drawBox3D(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(2.5f);
    glColor3f(0.8f, 0.8f, 0.8f);

    glBegin(GL_LINES);

    glVertex3f(ox, oy, oz);
    glVertex3f(ox + width, oy, oz);

    glVertex3f(ox, oy, oz);
    glVertex3f(ox, oy + height, oz);

    glVertex3f(ox, oy, oz);
    glVertex3f(ox, oy, oz + length);

    glVertex3f(ox + width, oy, oz);
    glVertex3f(ox + width, oy + height, oz);

    glVertex3f(ox + width, oy + height, oz);
    glVertex3f(ox, oy + height, oz);

    glVertex3f(ox, oy + height, oz + length);
    glVertex3f(ox, oy, oz + length);

    glVertex3f(ox, oy + height, oz + length);
    glVertex3f(ox, oy + height, oz);

    glVertex3f(ox + width, oy, oz);
    glVertex3f(ox + width, oy, oz + length);

    glVertex3f(ox, oy, oz + length);
    glVertex3f(ox + width, oy, oz + length);

    glVertex3f(ox + width, oy + height, oz);
    glVertex3f(ox + width, oy + height, oz + length);

    glVertex3f(ox + width, oy + height, oz + length);
    glVertex3f(ox + width, oy, oz + length);

    glVertex3f(ox, oy + height, oz + length);
    glVertex3f(ox + width, oy + height, oz + length);

    glEnd();
}

void drawTetMesh3D(FEMSimulator& simulator)
{
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.5f);
    glColor3f(0.9f, 0.1f, 0.1f);

    const model_tet& tetrahedra_meshes = simulator.getTetrahedraMeshes();
    const mesh3D& tetMesh = tetrahedra_meshes.mesh3Ds[0];
    const vector<Vector4i>& surf = tetMesh.surface;

    glBegin(GL_TRIANGLES);
    for (size_t j = 0; j < tetMesh.surface.size(); j++)
    {
        glVertex3f((tetMesh.vertexes[surf[j][0]][0]), (tetMesh.vertexes[surf[j][0]][1]), (tetMesh.vertexes[surf[j][0]][2]));
        glVertex3f((tetMesh.vertexes[surf[j][1]][0]), (tetMesh.vertexes[surf[j][1]][1]), (tetMesh.vertexes[surf[j][1]][2]));
        glVertex3f((tetMesh.vertexes[surf[j][2]][0]), (tetMesh.vertexes[surf[j][2]][1]), (tetMesh.vertexes[surf[j][2]][2]));
    }
    glEnd();

    glColor3f(0.9f, 0.9f, 0.9f);
    glLineWidth(0.1f);
    glBegin(GL_LINES);
    for (size_t j = 0; j < tetMesh.surfEdges.size(); j++)
    {
        glVertex3f((tetMesh.vertexes[tetMesh.surfEdges[j].first][0]), (tetMesh.vertexes[tetMesh.surfEdges[j].first][1]), (tetMesh.vertexes[tetMesh.surfEdges[j].first][2]));
        glVertex3f((tetMesh.vertexes[tetMesh.surfEdges[j].second][0]), (tetMesh.vertexes[tetMesh.surfEdges[j].second][1]), (tetMesh.vertexes[tetMesh.surfEdges[j].second][2]));

        glColor3f(0.9f, 0.9f, 0.9f);
        glLineWidth(0.1f);
    }
    glEnd();
}

void drawScene3D(float xTrans, float yTrans, float zTrans, float xRot, float yRot, FEMSimulator& simulator)
{
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(xTrans, yTrans, zTrans);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);

    drawBox3D(-1, -1, -1, 2, 2, 2);
    drawTetMesh3D(simulator);
    glPopMatrix();

    glutSwapBuffers();
}

} // namespace

GLViewer::GLViewer() = default;

void GLViewer::init(int argc, char** argv)
{
    g_viewer = this;

    glutInit(&argc, argv);
    glutSetOption(GLUT_MULTISAMPLE, 16);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);

    glutInitWindowSize(static_cast<int>(windowWidth_), static_cast<int>(windowHeight_));
    glutInitWindowPosition(0, 0);
    glutCreateWindow("FEM");

    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
    }
    else
    {
        std::cerr << "Status: Using GLEW " << glewGetString(GLEW_VERSION) << std::endl;
    }
    glClearColor(0.0, 0.0, 0.0, 1.0);

    simulator_.buildModels(0, 3);

    if (!isSetShader_)
    {
        glViewport(0, 0, static_cast<int>(windowWidth_), static_cast<int>(windowHeight_));
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, static_cast<float>(windowWidth_) / windowHeight_, 10.1f, 500.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -3.0f);
    }
    else
    {
        glGenBuffers(1, &pnVbo_);
        glGenVertexArrays(1, &vao_);
        loadShaders();
    }

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);

    glutDisplayFunc(displayThunk);
    glutReshapeFunc(reshapeThunk);
    glutKeyboardFunc(keyboardThunk);
    glutMouseFunc(mouseThunk);
    glutMotionFunc(motionThunk);
    glutIdleFunc(idleThunk);
}

void GLViewer::runMainLoop()
{
    glutMainLoop();
}

bool GLViewer::writeBitmapFile(int width, int height, const std::string& file_name, std::vector<unsigned char>& bitmapData)
{
    mBITMAPFILEHEADER bitmapFileHeader;
    memset(&bitmapFileHeader, 0, sizeof(mBITMAPFILEHEADER));
    bitmapFileHeader.bfSize = sizeof(mBITMAPFILEHEADER);
    bitmapFileHeader.bfType = 0x4d42;  //BM
    bitmapFileHeader.bfOffBits = sizeof(mBITMAPFILEHEADER) + sizeof(mBITMAPINFOHEADER);

    mBITMAPINFOHEADER bitmapInfoHeader;
    memset(&bitmapInfoHeader, 0, sizeof(mBITMAPINFOHEADER));
    bitmapInfoHeader.biSize = sizeof(mBITMAPINFOHEADER);
    bitmapInfoHeader.biWidth = width;
    bitmapInfoHeader.biHeight = height;
    bitmapInfoHeader.biPlanes = 1;
    bitmapInfoHeader.biBitCount = 24;
    bitmapInfoHeader.biCompression = 0L;
    bitmapInfoHeader.biSizeImage = width * std::abs(height) * 3;

    unsigned char tempRGB;
    for (int imageIdx = 0; imageIdx < static_cast<int>(bitmapInfoHeader.biSizeImage); imageIdx += 3)
    {
        tempRGB = bitmapData[imageIdx];
        bitmapData[imageIdx] = bitmapData[imageIdx + 2];
        bitmapData[imageIdx + 2] = tempRGB;
    }

    std::ofstream file(file_name, std::ios::binary);
    if (!file)
    {
        return false;
    }

    file.write(reinterpret_cast<const char*>(&bitmapFileHeader), sizeof(mBITMAPFILEHEADER));
    file.write(reinterpret_cast<const char*>(&bitmapInfoHeader), sizeof(mBITMAPINFOHEADER));
    file.write(reinterpret_cast<const char*>(bitmapData.data()), bitmapInfoHeader.biSizeImage);

    return file.good();
}

void GLViewer::saveScreenshot(const std::string& file_name)
{
    int data_len = static_cast<int>(windowHeight_ * windowWidth_ * 3); // bytes
    std::vector<unsigned char> screenData(data_len, 0);
    glReadPixels(0, 0, static_cast<int>(windowWidth_), static_cast<int>(windowHeight_), GL_RGB, GL_UNSIGNED_BYTE, screenData.data());

    writeBitmapFile(static_cast<int>(windowWidth_), static_cast<int>(windowHeight_), file_name + ".bmp", screenData);
}

void GLViewer::loadShaders()
{
    std::string vsPath = std::string(CIPC_ASSETS_DIR) + "shader/shader.vs";
    std::string fsPath = std::string(CIPC_ASSETS_DIR) + "shader/shader.fs";

    std::ifstream vsFile(vsPath);
    if (!vsFile)
    {
        std::cerr << "Failed to open vertex shader: " << vsPath << std::endl;
        return;
    }

    std::ifstream fsFile(fsPath);
    if (!fsFile)
    {
        std::cerr << "Failed to open fragment shader: " << fsPath << std::endl;
        return;
    }

    std::string vsSource((std::istreambuf_iterator<char>(vsFile)), std::istreambuf_iterator<char>());
    std::string fsSource((std::istreambuf_iterator<char>(fsFile)), std::istreambuf_iterator<char>());

    vertexShader_ = glCreateShader(GL_VERTEX_SHADER);
    fragmentShader_ = glCreateShader(GL_FRAGMENT_SHADER);

    const char* vsSourceCStr = vsSource.c_str();
    const char* fsSourceCStr = fsSource.c_str();

    glShaderSource(vertexShader_, 1, &vsSourceCStr, nullptr);
    glShaderSource(fragmentShader_, 1, &fsSourceCStr, nullptr);

    int success = 0;

    glCompileShader(vertexShader_);
    glGetShaderiv(vertexShader_, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char info_log[5000];
        glGetShaderInfoLog(vertexShader_, 5000, nullptr, info_log);
        std::cerr << "Error in vertex shader compilation!" << std::endl;
        std::cerr << "Info Log: " << info_log << std::endl;
    }

    glCompileShader(fragmentShader_);
    glGetShaderiv(fragmentShader_, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char info_log[5000];
        glGetShaderInfoLog(fragmentShader_, 5000, nullptr, info_log);
        std::cerr << "Error in fragment shader compilation!" << std::endl;
        std::cerr << "Info Log: " << info_log << std::endl;
    }

    shaderProgram_ = glCreateProgram();
    glAttachShader(shaderProgram_, vertexShader_);
    glAttachShader(shaderProgram_, fragmentShader_);
    glLinkProgram(shaderProgram_);
    glUseProgram(shaderProgram_);
}

void GLViewer::saveSurfaceMesh(const std::string& path)
{
    std::stringstream ss;
    ss << path;
    ss.fill('0');
    ss.width(5);
    ss << surfNumId_;
    surfNumId_++;
    ss << ".obj";
    std::string file_path = ss.str();

    ofstream outSurf(file_path);
    model_tet tetrahedra_meshes = simulator_.getTetrahedraMeshes();
    mesh3D& meshTemp = tetrahedra_meshes.mesh3Ds[0];
    map<int, int> meshToSurf;
    for (size_t i = 0; i < meshTemp.surfVerts.size(); i++)
    {
        const auto& pos = meshTemp.vertexes[meshTemp.surfVerts[i]];
        outSurf << "v " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
        meshToSurf[meshTemp.surfVerts[i]] = static_cast<int>(i);
    }

    for (size_t i = 0; i < meshTemp.surface.size(); i++)
    {
        const auto& tri = meshTemp.surface[i];
        outSurf << "f " << meshToSurf[tri[0]] + 1 << " " << meshToSurf[tri[1]] + 1 << " " << meshToSurf[tri[2]] + 1 << endl;
    }
    outSurf.close();
}

void GLViewer::onDisplay()
{
    drawScene3D(xTrans_, yTrans_, zTrans_, xRot_, yRot_, simulator_);

    int step = 0;
    int k = simulator_.simulateStick(step);

    if (false)
    {
        std::stringstream ss;
        ss << "saveScreen/step_";
        ss.fill('0');
        ss.width(5);
        ss << step / 10;
        if (step % 10 == 0)
        {
            saveScreenshot(ss.str());
        }
    }

    if (true)
    {
        saveSurfaceMesh("saveSurface/surf_");
    }

    printf("current step:  %d\n", step);

    if (step >= 119)
    {
    }
}

void GLViewer::onIdle()
{
    glutPostRedisplay();
}

void GLViewer::onReshape(int width, int height)
{
    windowWidth_ = static_cast<float>(width);
    windowHeight_ = static_cast<float>(height);

    glViewport(0, 0, width, height);
    if (!isSetShader_)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, static_cast<float>(width) / height, 0.1, 500.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -3.0f);
    }
}

void GLViewer::onKeyboard(unsigned char key, int /*x*/, int /*y*/)
{
    if (key == 'w')
    {
        zTrans_ += .01f;
    }

    if (key == 's')
    {
        zTrans_ -= .01f;
    }

    if (key == 'a')
    {
        xTrans_ += .01f;
    }

    if (key == 'd')
    {
        xTrans_ -= .01f;
    }

    if (key == 'q')
    {
        yTrans_ -= .01f;
    }

    if (key == 'e')
    {
        yTrans_ += .01f;
    }

    if (key == ' ')
    {
        stop_ = !stop_;
    }

    if (key == '/')
    {
        screenshot_ = !screenshot_;
    }

    if (key == 'k')
    {
        drawSurface_ = !drawSurface_;
    }

    if (key == 'f')
    {
        drawFiber_ = !drawFiber_;
    }

    if (key == '1')
    {
        drawMuscle_ = !drawMuscle_;
    }

    if (key == '2')
    {
        drawFat_ = !drawFat_;
    }

    if (key == '3')
    {
        drawSkull_ = !drawSkull_;
    }

    if (key == '4')
    {
        drawJaw_ = !drawJaw_;
    }

    if (key == '5')
    {
        drawMouth_ = !drawMouth_;
    }

    if (key == '9')
    {
        saveSurface_ = !saveSurface_;
    }

    if (key == 'm')
    {
        mouthOnly_ = !mouthOnly_;
    }

    if (key == '0')
    {
        isStop_ = !isStop_;
    }
    glutPostRedisplay();
}

void GLViewer::onMouse(int /*button*/, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        buttonState_ = 1;
    }
    else if (state == GLUT_UP)
    {
        buttonState_ = 0;
    }

    ox_ = x; oy_ = y;

    glutPostRedisplay();
}

void GLViewer::onMotion(int x, int y)
{
    float dx = static_cast<float>(x - ox_);
    float dy = static_cast<float>(y - oy_);

    if (buttonState_ == 1)
    {
        xRot_ += dy / 5.0f;
        yRot_ += dx / 5.0f;
    }

    ox_ = x; oy_ = y;

    glutPostRedisplay();
}
