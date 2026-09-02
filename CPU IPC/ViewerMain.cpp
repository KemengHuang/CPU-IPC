#include <GL/freeglut.h>

#include "RuntimePaths.h"
#include "Simulator.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int initialWindowWidth = 1000;
constexpr int initialWindowHeight = 1000;

FEMSimulator simulator;

float rotationX = 0.0f;
float rotationY = 0.0f;
float translationX = 0.0f;
float translationY = 0.0f;
float translationZ = 0.0f;
int previousMouseX = 0;
int previousMouseY = 0;
int currentWindowWidth = initialWindowWidth;
int currentWindowHeight = initialWindowHeight;
int surfaceFrame = 0;
bool mouseDragging = false;
bool simulationPaused = true;
bool screenshotsEnabled = false;
bool surfaceExportEnabled = false;

SimulationScene extractViewerSceneArgument(int& argc, char** argv)
{
    SimulationScene scene = SimulationScene::ClothOverBunny;
    int outputArgument = 1;
    for (int inputArgument = 1; inputArgument < argc; ++inputArgument) {
        const std::string argument = argv[inputArgument];
        if (argument == "--scene") {
            if (++inputArgument >= argc) {
                throw std::invalid_argument("--scene requires a value");
            }
            scene = parseSimulationScene(argv[inputArgument]);
        }
        else {
            argv[outputArgument++] = argv[inputArgument];
        }
    }
    argc = outputArgument;
    argv[outputArgument] = nullptr;
    return scene;
}

#pragma pack(push, 1)
struct BitmapFileHeader {
    std::uint16_t type;
    std::uint32_t size;
    std::uint16_t reserved1;
    std::uint16_t reserved2;
    std::uint32_t pixelOffset;
};

struct BitmapInfoHeader {
    std::uint32_t size;
    std::int32_t width;
    std::int32_t height;
    std::uint16_t planes;
    std::uint16_t bitsPerPixel;
    std::uint32_t compression;
    std::uint32_t imageSize;
    std::int32_t horizontalResolution;
    std::int32_t verticalResolution;
    std::uint32_t colorsUsed;
    std::uint32_t importantColors;
};
#pragma pack(pop)

bool writeBitmap(
    int width,
    int height,
    const std::string& fileName,
    const std::vector<unsigned char>& pixels)
{
    const std::uint32_t imageSize = static_cast<std::uint32_t>(pixels.size());
    BitmapFileHeader fileHeader{};
    fileHeader.type = 0x4d42;
    fileHeader.pixelOffset = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);
    fileHeader.size = fileHeader.pixelOffset + imageSize;

    BitmapInfoHeader infoHeader{};
    infoHeader.size = sizeof(BitmapInfoHeader);
    infoHeader.width = width;
    infoHeader.height = height;
    infoHeader.planes = 1;
    infoHeader.bitsPerPixel = 24;
    infoHeader.imageSize = imageSize;

    std::FILE* file = std::fopen(fileName.c_str(), "wb");
    if (file == nullptr) {
        return false;
    }

    const bool written =
        std::fwrite(&fileHeader, sizeof(fileHeader), 1, file) == 1
        && std::fwrite(&infoHeader, sizeof(infoHeader), 1, file) == 1
        && std::fwrite(pixels.data(), pixels.size(), 1, file) == 1;
    std::fclose(file);
    return written;
}

bool captureScreenshot(int width, int height, const std::string& fileName)
{
    if (width <= 0 || height <= 0) {
        return false;
    }

    const std::size_t sourceRowBytes = static_cast<std::size_t>(width) * 3;
    const std::size_t bitmapRowBytes = (sourceRowBytes + 3) & ~std::size_t(3);
    std::vector<unsigned char> rgb(sourceRowBytes * static_cast<std::size_t>(height));
    std::vector<unsigned char> bgr(bitmapRowBytes * static_cast<std::size_t>(height), 0);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, rgb.data());
    for (int row = 0; row < height; ++row) {
        for (int column = 0; column < width; ++column) {
            const std::size_t source = static_cast<std::size_t>(row) * sourceRowBytes
                + static_cast<std::size_t>(column) * 3;
            const std::size_t destination = static_cast<std::size_t>(row) * bitmapRowBytes
                + static_cast<std::size_t>(column) * 3;
            bgr[destination] = rgb[source + 2];
            bgr[destination + 1] = rgb[source + 1];
            bgr[destination + 2] = rgb[source];
        }
    }

    return writeBitmap(width, height, fileName + ".bmp", bgr);
}

void drawBoundingBox(
    float originX,
    float originY,
    float originZ,
    float width,
    float height,
    float depth)
{
    const float x1 = originX + width;
    const float y1 = originY + height;
    const float z1 = originZ + depth;
    const float corners[8][3] = {
        { originX, originY, originZ }, { x1, originY, originZ },
        { originX, y1, originZ }, { x1, y1, originZ },
        { originX, originY, z1 }, { x1, originY, z1 },
        { originX, y1, z1 }, { x1, y1, z1 }
    };
    constexpr int edges[12][2] = {
        { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 3 },
        { 1, 5 }, { 2, 3 }, { 2, 6 }, { 3, 7 },
        { 4, 5 }, { 4, 6 }, { 5, 7 }, { 6, 7 }
    };

    glLineWidth(2.5f);
    glColor3f(0.8f, 0.8f, 0.8f);
    glBegin(GL_LINES);
    for (const auto& edge : edges) {
        glVertex3fv(corners[edge[0]]);
        glVertex3fv(corners[edge[1]]);
    }
    glEnd();
}

void drawSimulationSurface()
{
    const SimulationModel& model = simulator.getModel();
    if (model.meshes.empty()) {
        return;
    }

    const mesh3D& mesh = model.meshes.front();
    glColor3f(0.9f, 0.1f, 0.1f);
    glBegin(GL_TRIANGLES);
    for (const Eigen::Vector4i& face : mesh.surface) {
        for (int corner = 0; corner < 3; ++corner) {
            const Eigen::Vector3d& position = mesh.vertexes[face[corner]];
            glVertex3d(position.x(), position.y(), position.z());
        }
    }
    glEnd();

    glColor3f(0.9f, 0.9f, 0.9f);
    glLineWidth(0.1f);
    glBegin(GL_LINES);
    for (const auto& edge : mesh.surfEdges) {
        const Eigen::Vector3d& first = mesh.vertexes[edge.first];
        const Eigen::Vector3d& second = mesh.vertexes[edge.second];
        glVertex3d(first.x(), first.y(), first.z());
        glVertex3d(second.x(), second.y(), second.z());
    }
    glEnd();
}

void renderScene()
{
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(translationX, translationY, translationZ);
    glRotatef(rotationX, 1.0f, 0.0f, 0.0f);
    glRotatef(rotationY, 0.0f, 1.0f, 0.0f);
    drawBoundingBox(-1.0f, -1.0f, -1.0f, 2.0f, 2.0f, 2.0f);
    drawSimulationSurface();
    glPopMatrix();
    glutSwapBuffers();
}

void saveScreenshotIfDue(int step)
{
    if (step % 10 != 0) {
        return;
    }

    std::ostringstream name;
    name << RuntimePaths::screenshotFile("step_");
    name.fill('0');
    name.width(5);
    name << step / 10;
    if (!captureScreenshot(currentWindowWidth, currentWindowHeight, name.str())) {
        std::cerr << "Warning: failed to save screenshot " << name.str() << ".bmp\n";
    }
}

void saveSurfaceMesh()
{
    const SimulationModel& model = simulator.getModel();
    if (model.meshes.empty()) {
        return;
    }

    std::ostringstream name;
    name << RuntimePaths::surfaceFile("surf_");
    name.fill('0');
    name.width(5);
    name << surfaceFrame++ << ".obj";

    std::ofstream output(name.str());
    if (!output) {
        std::cerr << "Warning: failed to write surface mesh to " << name.str() << '\n';
        return;
    }

    const mesh3D& mesh = model.meshes.front();
    std::vector<int> meshToSurface(static_cast<std::size_t>(mesh.vertexNum), -1);
    for (int surfaceVertex = 0;
         surfaceVertex < static_cast<int>(mesh.surfVerts.size());
         ++surfaceVertex) {
        const int vertex = mesh.surfVerts[surfaceVertex];
        const Eigen::Vector3d& position = mesh.vertexes[vertex];
        output << "v " << position.x() << ' ' << position.y() << ' '
               << position.z() << '\n';
        meshToSurface[vertex] = surfaceVertex;
    }
    for (const Eigen::Vector4i& face : mesh.surface) {
        output << "f " << meshToSurface[face[0]] + 1 << ' '
               << meshToSurface[face[1]] + 1 << ' '
               << meshToSurface[face[2]] + 1 << '\n';
    }
}

void display()
{
    renderScene();
    if (simulationPaused) {
        return;
    }

    int step = 0;
    simulator.simulateStep(step);
    if (screenshotsEnabled) {
        saveScreenshotIfDue(step);
    }
    if (surfaceExportEnabled) {
        saveSurfaceMesh();
    }
    std::cout << "current step: " << step << '\n';
}

void configureProjection(int width, int height)
{
    currentWindowWidth = (std::max)(1, width);
    currentWindowHeight = (std::max)(1, height);
    glViewport(0, 0, currentWindowWidth, currentWindowHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(
        45.0,
        static_cast<double>(currentWindowWidth) / currentWindowHeight,
        0.1,
        500.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -3.0f);
}

void initializeViewer(SimulationScene scene)
{
    simulator.buildModels(scene);
    configureProjection(initialWindowWidth, initialWindowHeight);
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
}

void idle()
{
    glutPostRedisplay();
}

void reshape(int width, int height)
{
    configureProjection(width, height);
}

void keyboard(unsigned char key, int, int)
{
    switch (key) {
    case 'w': translationZ += 0.01f; break;
    case 's': translationZ -= 0.01f; break;
    case 'a': translationX += 0.01f; break;
    case 'd': translationX -= 0.01f; break;
    case 'q': translationY -= 0.01f; break;
    case 'e': translationY += 0.01f; break;
    case ' ': simulationPaused = !simulationPaused; break;
    case '/': screenshotsEnabled = !screenshotsEnabled; break;
    case '9': surfaceExportEnabled = !surfaceExportEnabled; break;
    default: break;
    }
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        mouseDragging = state == GLUT_DOWN;
    }
    previousMouseX = x;
    previousMouseY = y;
    glutPostRedisplay();
}

void motion(int x, int y)
{
    if (mouseDragging) {
        rotationX += static_cast<float>(y - previousMouseY) / 5.0f;
        rotationY += static_cast<float>(x - previousMouseX) / 5.0f;
    }
    previousMouseX = x;
    previousMouseY = y;
    glutPostRedisplay();
}

} // namespace

int main(int argc, char** argv)
{
    SimulationScene scene = SimulationScene::ClothOverBunny;
    try {
        scene = extractViewerSceneArgument(argc, argv);
    }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        std::cerr << "Usage: " << argv[0]
                  << " [--scene cloth-bunny|twisting-mat|twisting-mat-soft|bunny2]\n";
        return 1;
    }

    glutInit(&argc, argv);
    glutSetOption(GLUT_MULTISAMPLE, 16);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutInitWindowSize(initialWindowWidth, initialWindowHeight);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("CPU IPC");

    initializeViewer(scene);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);
    glutMainLoop();
    return 0;
}
