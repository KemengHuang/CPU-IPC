#include "viewer/GLViewer.h"

#include "GL/glew.h"
#include "GL/freeglut.h"
#include "mesh/Mesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

using namespace std;
using namespace Eigen;

int surfNumId = 0;
float xRot = 0.0f;
float yRot = 0.f;
float xTrans = 0;
float yTrans = 0;
float zTrans = 0;
int ox;
int oy;
int buttonState;
float xRotLength = 0.0f;
float yRotLength = 0.0f;
float window_width = 1000;
float window_height = 1000;
int s_dimention = 3;
bool stop = true;
bool screenshot = true;
bool saveSurface = false;
bool mouthOnly = false;
bool isSetShader = false;
bool isStop = false;

fiber_obj fiberObj;

bool drawFiber = false;
bool drawSurface = false;
bool drawMuscle = true;
bool drawFat = true;
bool drawSkull = true;
bool drawJaw = true;
bool drawMouth = true;

FEMSimulator simulator;

GLuint PN_vbo_;
GLuint VAO;
GLuint color_vbo_;
//GLuint color_vao_;
GLuint normal_vbo_;
//GLuint normal_vao_;
GLuint v;
GLuint f;
GLuint shaderProgram;


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

bool WriteBitmapFile(int width, int height, const std::string& file_name, unsigned char* bitmapData)
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
    bitmapInfoHeader.biSizeImage = width * abs(height) * 3;

    //////////////////////////////////////////////////////////////////////////
    FILE* filePtr;
    unsigned char tempRGB;
    int           imageIdx;

    for (imageIdx = 0; imageIdx < (int)bitmapInfoHeader.biSizeImage; imageIdx += 3)
    {
        tempRGB = bitmapData[imageIdx];
        bitmapData[imageIdx] = bitmapData[imageIdx + 2];
        bitmapData[imageIdx + 2] = tempRGB;
    }

    filePtr = fopen(file_name.c_str(), "wb");
    if (NULL == filePtr)
    {
        return false;
    }

    fwrite(&bitmapFileHeader, sizeof(mBITMAPFILEHEADER), 1, filePtr);

    fwrite(&bitmapInfoHeader, sizeof(mBITMAPINFOHEADER), 1, filePtr);

    fwrite(bitmapData, bitmapInfoHeader.biSizeImage, 1, filePtr);

    fclose(filePtr);
    return true;
}

void SaveScreenShot(int width, int height, const std::string& file_name)
{
    int data_len = height * width * 3;      // bytes
    void* screen_data = malloc(data_len);
    memset(screen_data, 0, data_len);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, screen_data);

    WriteBitmapFile(width, height, file_name + ".bmp", (unsigned char*)screen_data);

    free(screen_data);
}

void set_shaders()
{
    char* vs = NULL;
    char* fs = NULL;

    vs = (char*)malloc(sizeof(char) * 10000);
    fs = (char*)malloc(sizeof(char) * 10000);
    memset(vs, 0, sizeof(char) * 10000);
    memset(fs, 0, sizeof(char) * 10000);

    FILE* fp;
    char c;
    int count;

    fp = fopen((std::string(CIPC_ASSETS_DIR) + "shader/shader.vs").c_str(), "r");
    count = 0;
    while ((c = fgetc(fp)) != EOF)
    {
        vs[count] = c;
        count++;
    }
    fclose(fp);

    fp = fopen((std::string(CIPC_ASSETS_DIR) + "shader/shader.fs").c_str(), "r");
    count = 0;
    while ((c = fgetc(fp)) != EOF)
    {
        fs[count] = c;
        count++;
    }
    fclose(fp);

    v = glCreateShader(GL_VERTEX_SHADER);
    f = glCreateShader(GL_FRAGMENT_SHADER);

    const char* vv;
    const char* ff;
    vv = vs;
    ff = fs;

    glShaderSource(v, 1, &vv, NULL);
    glShaderSource(f, 1, &ff, NULL);

    int success;

    glCompileShader(v);
    glGetShaderiv(v, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char info_log[5000];
        glGetShaderInfoLog(v, 5000, NULL, info_log);
        printf("Error in vertex shader compilation!\n");
        printf("Info Log: %s\n", info_log);
    }

    glCompileShader(f);
    glGetShaderiv(f, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char info_log[5000];
        glGetShaderInfoLog(f, 5000, NULL, info_log);
        printf("Error in fragment shader compilation!\n");
        printf("Info Log: %s\n", info_log);
    }

    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, v);
    glAttachShader(shaderProgram, f);
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    free(vs);
    free(fs);
}

void draw_box2D(float ox, float oy, float width, float height)
{
    glLineWidth(2.5f);
    glColor3f(0.8f, 0.8f, 0.8f);

    glBegin(GL_LINES);

    glVertex3f(ox, oy, 0);
    glVertex3f(ox + width, oy, 0);

    glVertex3f(ox, oy, 0);
    glVertex3f(ox, oy + height, 0);

    glVertex3f(ox + width, oy, 0);
    glVertex3f(ox + width, oy + height, 0);

    glVertex3f(ox + width, oy + height, 0);
    glVertex3f(ox, oy + height, 0);

    glEnd();
}

void draw_box3D(float ox, float oy, float oz, float width, float height, float length)
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

void draw_lines(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(0.5f);
    glColor3f(0.8f, 0.8f, 0.8f);

    glBegin(GL_LINES);
    int numbers = 20;
    for (int i = 0; i <= numbers; i++) {
        //glVertex3f(ox, oy, oz);
        glVertex3f(ox + width * i / numbers, oy, 0);
        glVertex3f(ox + width * i / numbers, oy + height, 0);
    }

    for (int i = 0; i <= numbers; i++) {
        //glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy + height * i / numbers, 0);
        glVertex3f(ox + width, oy + height * i / numbers, 0);
    }

    glEnd();


    glLineWidth(1.5f);
    glColor3f(0.8f, 0.8f, 0.f);
    glBegin(GL_LINES);
    glVertex3f(ox + width / 2, oy, 0);
    glVertex3f(ox + width / 2, oy + height, 0);

    glVertex3f(ox, oy + height / 2, 0);
    glVertex3f(ox + width, oy + height / 2, 0);

    glEnd();
}

void draw_tet_mesh3D()
{
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.5f);
    glColor3f(0.9f, 0.1f, 0.1f);
    const model_tet& tetrahedra_meshes = simulator.getTetrahedraMeshes();
    const mesh3D& tetMesh = tetrahedra_meshes.mesh3Ds[0];
    const vector<Vector4i>& surf = tetMesh.surface;//obj.faces;
    glBegin(GL_TRIANGLES);

    for (int j = 0; j < tetMesh.surface.size(); j++) {
        glVertex3f((tetMesh.vertexes[surf[j][0]][0]), (tetMesh.vertexes[surf[j][0]][1]), (tetMesh.vertexes[surf[j][0]][2]));
        glVertex3f((tetMesh.vertexes[surf[j][1]][0]), (tetMesh.vertexes[surf[j][1]][1]), (tetMesh.vertexes[surf[j][1]][2]));
        glVertex3f((tetMesh.vertexes[surf[j][2]][0]), (tetMesh.vertexes[surf[j][2]][1]), (tetMesh.vertexes[surf[j][2]][2]));
    }
    glEnd();

    glColor3f(0.9f, 0.9f, 0.9f);
    glLineWidth(0.1f);
    glBegin(GL_LINES);

    for (int j = 0; j < tetMesh.surfEdges.size(); j++) {
        glVertex3f((tetMesh.vertexes[tetMesh.surfEdges[j].first][0]), (tetMesh.vertexes[tetMesh.surfEdges[j].first][1]), (tetMesh.vertexes[tetMesh.surfEdges[j].first][2]));
        glVertex3f((tetMesh.vertexes[tetMesh.surfEdges[j].second][0]), (tetMesh.vertexes[tetMesh.surfEdges[j].second][1]), (tetMesh.vertexes[tetMesh.surfEdges[j].second][2]));

        glColor3f(0.9f, 0.9f, 0.9f);
        glLineWidth(0.1f);
    }
    glEnd();

}


void draw_Scene3D() {
    //face.mesh3Ds[0] = mesh3d;
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(xTrans, yTrans, zTrans);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);

    draw_box3D(-1, -1, -1, 2, 2, 2);
    //draw_lines(-1, -1, -1, 2, 2, 2);
    draw_tet_mesh3D();
    glPopMatrix();


    glutSwapBuffers();
    //glFlush();
}
double mfsum = 0;
double total_time = 0;
int total_cg_iterations = 0;
int total_newton_iterations = 0;
int start = -1;

void saveScreenPic(const string& path, int step_index) {
    std::stringstream ss;
    ss << path;
    ss.fill('0');
    ss.width(5);
    ss << step_index / 10;
    if (step_index % 10 != 0) return;
    std::string file_path = ss.str();
    //if (step / 10 != start) {
    //start = step / 10;
    SaveScreenShot(window_width, window_height, file_path);
    //}
}

void saveSurfaceMesh(const string& path) {
    std::stringstream ss;
    ss << path;
    ss.fill('0');
    ss.width(5);
    //if(surfNumId%10!=0) return;
    ss << surfNumId;
    surfNumId++;
    //if(surfNumId%3!=0) return;
    ss << ".obj";
    std::string file_path = ss.str();
    ofstream outSurf(file_path);
    model_tet tetrahedra_meshes = simulator.getTetrahedraMeshes();
    mesh3D& meshTemp = tetrahedra_meshes.mesh3Ds[0];
    map<int, int> meshToSurf;//(meshTemp.surfVerts.size());
    for (int i = 0; i < meshTemp.surfVerts.size(); i++) {
        const auto& pos = meshTemp.vertexes[meshTemp.surfVerts[i]];
        outSurf << "v " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
        meshToSurf[meshTemp.surfVerts[i]] = i;
    }

    for (int i = 0; i < meshTemp.surface.size(); i++) {
        const auto& tri = meshTemp.surface[i];
        outSurf << "f " << meshToSurf[tri[0]] + 1 << " " << meshToSurf[tri[1]] + 1 << " " << meshToSurf[tri[2]] + 1 << endl;
    }
    outSurf.close();
}
vector<int> newtonIt;

vector<Vector3d> tpv[5];

void display(void) {

    draw_Scene3D();


    //if (stop) return;

    int step = 0;
    int k = simulator.simulateStick(step);

    if (false) {
        saveScreenPic("saveScreen/step_", step);
    }

    if (true) {
        saveSurfaceMesh("saveSurface/surf_");
    }


    printf("current step:  %d\n", step);
    newtonIt.push_back(k);

    if (step>=119) {
    }
}

void initScene() {
    simulator.buildModels(0, 3);
}

void initGL()
{
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
    }
    std::cerr << "Status: Using GLEW " << glewGetString(GLEW_VERSION) << std::endl;
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

static void init(void)
{
    initGL();
    initScene();

    if (!isSetShader) {
        glViewport(0, 0, window_width, window_height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, (float)window_width / window_height, 10.1f, 500.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -3.0f);
    }
    else {
        glGenBuffers(1, &PN_vbo_);
        glGenVertexArrays(1, &VAO);
    }
    //glEnable(GL_DEPTH_TEST);
}

void idle_func()
{
    glutPostRedisplay();
}

void reshape_func(GLint width, GLint height)
{
    glViewport(0, 0, width, height);
    if (!isSetShader) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        gluPerspective(45.0, (float)width / height, 0.1, 500.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -3.0f);
    }
    //glTranslatef(0.5f, 0.5f, -4.0f);
}

void keyboard_func(unsigned char key, int x, int y)
{
    if (key == 'w')
    {
        zTrans += .01f;
    }

    if (key == 's')
    {
        zTrans -= .01f;
    }

    if (key == 'a')
    {
        xTrans += .01f;
    }

    if (key == 'd')
    {
        xTrans -= .01f;
    }

    if (key == 'q')
    {
        yTrans -= .01f;
    }

    if (key == 'e')
    {
        yTrans += .01f;
    }

    if (key == ' ')
    {
        stop = !stop;
    }

    if (key == '/')
    {
        screenshot = !screenshot;
    }

    if (key == 'k')
    {
        drawSurface = !drawSurface;
    }

    if (key == 'f')
    {
        drawFiber = !drawFiber;
    }

    if (key == '1')
    {
        drawMuscle = !drawMuscle;
    }

    if (key == '2')
    {
        drawFat = !drawFat;
    }

    if (key == '3')
    {
        drawSkull = !drawSkull;
    }

    if (key == '4')
    {
        drawJaw = !drawJaw;
    }

    if (key == '5')
    {
        drawMouth = !drawMouth;
    }

    if (key == '9')
    {
        saveSurface = !saveSurface;
    }

    if (key == 'm')
    {
        mouthOnly = !mouthOnly;
    }

    if (key == '0')
    {
        isStop = !isStop;
    }
    glutPostRedisplay();
}

void special_keyboard_func(int key, int x, int y)
{
    glutPostRedisplay();
}

void mouse_func(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        buttonState = 1;
    }
    else if (state == GLUT_UP)
    {
        buttonState = 0;
    }

    ox = x; oy = y;

    glutPostRedisplay();
}

void motion_func(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

    if (buttonState == 1)
    {
        xRot += dy / 5.0f;
        yRot += dx / 5.0f;
    }

    ox = x; oy = y;

    glutPostRedisplay();
}

void runGLUTMainLoop(int argc, char** argv)
{
    glutInit(&argc, argv);

    glutSetOption(GLUT_MULTISAMPLE, 16);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);

    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("FEM");

    init();

    if (isSetShader) {
        set_shaders();
    }

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);


    glEnable(GL_MULTISAMPLE);
    glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);


    glutDisplayFunc(display);


    //glutDisplayFunc(display_func);
    glutReshapeFunc(reshape_func);
    glutKeyboardFunc(keyboard_func);
    //glutSpecialFunc(special_keyboard_func);
    glutMouseFunc(mouse_func);
    glutMotionFunc(motion_func);
    glutIdleFunc(idle_func);


    glutMainLoop();
}
