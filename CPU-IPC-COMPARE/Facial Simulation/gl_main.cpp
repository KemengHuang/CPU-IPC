#include "GL\glew.h"
#include "GL\freeglut.h"
#include <fstream>
#include <iostream>
#include<cuda_runtime.h>
#include "fem3D.h"
#include "fem_timer.h"
#include "Simulator.h"
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
//#include <unsupported/Eigen/SparseExtra>
using namespace std;
int step = 0;
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
bool screenshot = false;
bool saveSurface = false;
bool mouthOnly = false;
bool isSetShader = false;
//mesh2D mesh2d;
//mesh3D mesh3d;
//model_obj triangle_meshes;
//model_tet tetrahedra_meshes;
//model_tet face2;
//mesh_obj insideBound1;
//mesh_obj insideBound2;

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




bool WriteBitmapFile(int width, int height, const std::string& file_name, unsigned char* bitmapData)
{
    BITMAPFILEHEADER bitmapFileHeader;
    memset(&bitmapFileHeader, 0, sizeof(BITMAPFILEHEADER));
    bitmapFileHeader.bfSize = sizeof(BITMAPFILEHEADER);
    bitmapFileHeader.bfType = 0x4d42;   //BM  
    bitmapFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);

    BITMAPINFOHEADER bitmapInfoHeader;
    memset(&bitmapInfoHeader, 0, sizeof(BITMAPINFOHEADER));
    bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
    bitmapInfoHeader.biWidth = width;
    bitmapInfoHeader.biHeight = height;
    bitmapInfoHeader.biPlanes = 1;
    bitmapInfoHeader.biBitCount = 24;
    bitmapInfoHeader.biCompression = BI_RGB;
    bitmapInfoHeader.biSizeImage = width * abs(height) * 3;

    //////////////////////////////////////////////////////////////////////////  
    FILE* filePtr;
    unsigned char tempRGB;
    int imageIdx;

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

    fwrite(&bitmapFileHeader, sizeof(BITMAPFILEHEADER), 1, filePtr);

    fwrite(&bitmapInfoHeader, sizeof(BITMAPINFOHEADER), 1, filePtr);

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

    fp = fopen("shader/shader.vs", "r");
    count = 0;
    while ((c = fgetc(fp)) != EOF)
    {
        vs[count] = c;
        count++;
    }
    fclose(fp);

    fp = fopen("shader/shader.fs", "r");
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

//void draw_mesh2D()
//{
//    glEnable(GL_DEPTH_TEST);
//    glLineWidth(1.5f);
//    glColor3f(0.8f, 0.1f, 0.8f);
//    glBegin(GL_LINES);
//    for (int i = 0; i < mesh2d.triangleNum; i++) {
//        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][0]][0], mesh2d.vertexes[mesh2d.triangles[i][0]][1], 0.0f);
//        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][1]][0], mesh2d.vertexes[mesh2d.triangles[i][1]][1], 0.0f);
//
//        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][0]][0], mesh2d.vertexes[mesh2d.triangles[i][0]][1], 0.0f);
//        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][2]][0], mesh2d.vertexes[mesh2d.triangles[i][2]][1], 0.0f);
//
//        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][1]][0], mesh2d.vertexes[mesh2d.triangles[i][1]][1], 0.0f);
//        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][2]][0], mesh2d.vertexes[mesh2d.triangles[i][2]][1], 0.0f);
//    }
//    glEnd();
//}

void draw_tet_mesh3D()
{
    glEnable(GL_DEPTH_TEST);
    glLineWidth(0.5f);
    glColor3f(0.8f, 0.1f, 0.8f);
    glBegin(GL_TRIANGLES);
    int drawtype = 2;
    model_tet tetrahedra_meshes = simulator.getTetrahedraMeshes();
    for (int f = 0; f < tetrahedra_meshes.mesh3Ds.size(); f++) {
        mesh3D& meshTemp = tetrahedra_meshes.mesh3Ds[f];
        for (int i = 0; i < meshTemp.tetrahedraNum; i++) {
            //if (meshTemp.isInside[i] != drawtype) continue;
            //Vector3d center = meshTemp.vertexes[meshTemp.tetrahedras[i][0]] + meshTemp.vertexes[meshTemp.tetrahedras[i][1]] + meshTemp.vertexes[meshTemp.tetrahedras[i][2]] + meshTemp.vertexes[meshTemp.tetrahedras[i][3]];
            //center = center / 4;
            ////if (center[0] < 0) //continue;
            //if (true)
            //{
            //    if (meshTemp.isInside[i] == -1 && !drawFat) {
            //        continue;
            //    }

            //    if (meshTemp.isMuscle[i] && !drawMuscle) {
            //        continue;
            //    }
            //}
            //if (meshTemp.isJaw[i] && !drawJaw) {
            //    continue;
            //}
            //if (meshTemp.isSkull[i] && !drawSkull) {
            //    continue;
            //}
            //if (meshTemp.isMouth[i] && !drawMouth) {
            //    continue;
            //}

            float red = 0, green = 0, blue = 0;
                red = 0.3;
 
                blue = 0.6;

                green = 0.3;
            glColor3f(red, green, blue);
            //if (meshTemp.rehabilitate[i] == -1) {

                //glColor3f(red, green, blue);
                for (int j = 0; j < 4; j++) {

                    glColor3f(red, green, blue);
                    if ((j) % 4 == 2) glColor3f(red / 2, green / 2, blue / 2);
                    if (meshTemp.Constraints[meshTemp.tetrahedras[i][j]].determinant() < 0.5) glColor3f(1, 1, 1);
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][j]][0], meshTemp.vertexes[meshTemp.tetrahedras[i][j]][1], meshTemp.vertexes[meshTemp.tetrahedras[i][j]][2]);
                    glColor3f(red, green, blue);
                    if ((j + 1) % 4 == 2) glColor3f(red / 2, green / 2, blue / 2);
                    if (meshTemp.Constraints[meshTemp.tetrahedras[i][(j + 1) % 4]].determinant() < 0.5) glColor3f(1, 1, 1);
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][0], meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][1], meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][2]);
                    glColor3f(red, green, blue);
                    if ((j + 2) % 4 == 2) glColor3f(red / 2, green / 2, blue / 2);
                    if (meshTemp.Constraints[meshTemp.tetrahedras[i][(j + 2) % 4]].determinant() < 0.5) glColor3f(1, 1, 1);
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][0], meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][1], meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][2]);

                }
                //break;
            //}
        }
    }
    glEnd();
    glColor3f(0.8f, 0.8f, 0.1f);
    //glDisable(GL_DEPTH_TEST);
    glLineWidth(0.5f);
    glBegin(GL_LINES);
    double offset = 1.0;
    for (int f = 0; f < tetrahedra_meshes.mesh3Ds.size(); f++) {
        mesh3D& meshTemp = tetrahedra_meshes.mesh3Ds[f];
        for (int i = 0; i < meshTemp.tetrahedraNum; i++) {
                for (int j = 0; j < 4; j++) {
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][j]][0] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][j]][1] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][j]][2] * offset);
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][0] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][1] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][2] * offset);

                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][0] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][1] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 1) % 4]][2] * offset);
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][0] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][1] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][2] * offset);

                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][j]][0] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][j]][1] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][j]][2] * offset);
                    glVertex3f(meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][0] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][1] * offset, meshTemp.vertexes[meshTemp.tetrahedras[i][(j + 2) % 4]][2] * offset);
                }

        }
    }
    glEnd();
}




//bool drawSkin = false;
//void draw_mesh3D()
//{
//    glEnable(GL_DEPTH_TEST);
//    glLineWidth(1.5f);
//    glColor3f(0.5f, 0.5f, 0.5f);
//    glBegin(GL_TRIANGLES);
//    for (int i = 0; i < triangle_meshes.meshes.size(); i++) {
//        if (!drawSkin && i == 2) continue;
//        for (int j = 0; j < triangle_meshes.meshes[i].faceNum; j++) {
//            for (int t = 0; t < 3; t++) {
//                int p = triangle_meshes.meshes[i].faces[j][t];
//               
//                glVertex3f(triangle_meshes.meshes[i].vertexes[p][0], triangle_meshes.meshes[i].vertexes[p][1], triangle_meshes.meshes[i].vertexes[p][2]);
//            }
//        }
//    }
//    glEnd();
//
//    glColor3f(0.8f, 0.8f, 0.1f);
//    //glDisable(GL_DEPTH_TEST);
//    glLineWidth(0.1f);
//    glBegin(GL_LINES);
//    double offset = 1.0;
//    for (int i = 0; i < triangle_meshes.meshes.size(); i++) {
//        if (!drawSkin && i == 2) continue;
//        for (int j = 0; j < triangle_meshes.meshes[i].faceNum; j++) {
//            //for (int t = 0; t < 3; t++) {
//            int p0 = triangle_meshes.meshes[i].faces[j][0];
//            int p1 = triangle_meshes.meshes[i].faces[j][1];
//            int p2 = triangle_meshes.meshes[i].faces[j][2];
//
//            //glVertex3f(xx, yy, zz);
//            glVertex3f(triangle_meshes.meshes[i].vertexes[p0][0], triangle_meshes.meshes[i].vertexes[p0][1], triangle_meshes.meshes[i].vertexes[p0][2]);
//            glVertex3f(triangle_meshes.meshes[i].vertexes[p1][0], triangle_meshes.meshes[i].vertexes[p1][1], triangle_meshes.meshes[i].vertexes[p1][2]);
//
//            glVertex3f(triangle_meshes.meshes[i].vertexes[p1][0], triangle_meshes.meshes[i].vertexes[p1][1], triangle_meshes.meshes[i].vertexes[p1][2]);
//            glVertex3f(triangle_meshes.meshes[i].vertexes[p2][0], triangle_meshes.meshes[i].vertexes[p2][1], triangle_meshes.meshes[i].vertexes[p2][2]);
//
//            glVertex3f(triangle_meshes.meshes[i].vertexes[p2][0], triangle_meshes.meshes[i].vertexes[p2][1], triangle_meshes.meshes[i].vertexes[p2][2]);
//            glVertex3f(triangle_meshes.meshes[i].vertexes[p0][0], triangle_meshes.meshes[i].vertexes[p0][1], triangle_meshes.meshes[i].vertexes[p0][2]);
//            //}
//        }
//    }
//    glEnd();
//}

void draw_mesh3D()
{
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.5f);
    glColor3f(0.5f, 0.5f, 0.5f);
    glBegin(GL_TRIANGLES);
    model_tet tetrahedra_meshes = simulator.getTetrahedraMeshes();
    for (int f = 0; f < tetrahedra_meshes.mesh3Ds.size(); f++) {

        const vector<Vector4i>& surf = tetrahedra_meshes.mesh3Ds[f].surface;

        for (int j = 0; j < surf.size(); j++) {
            if (mouthOnly) {
                //Vector3d center = tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]] + tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]] + tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]];
                //center /= 3;
                //if (center.y() > 0.0 || center.y() < -0.35) continue;
                //if (center.z() < 0.0) continue;
                //if (center.x() > 0.23 || center.x() < -0.23) continue;
                //if (tetrahedra_meshes.mesh3Ds[f].isJaw[surf[j][3]] || tetrahedra_meshes.mesh3Ds[f].isSkull[surf[j][3]]) {
                //    glColor3f(0.8f, 0.1f, 0.8f);
                //}
                //else {
                //    continue;//glColor3f(0.8f, 0.8f, 0.1f);
                //}

                for (int t = 0; t < 3; t++) {
                    glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][t]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][t]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][t]][2]);
                }
            }
        }
    }
    glEnd();

    glColor3f(0.9f, 0.9f, 0.9f);
    //glDisable(GL_DEPTH_TEST);
    glLineWidth(0.1f);
    glBegin(GL_LINES);
    //tetrahedra_meshes = simulator.getTetrahedraMeshes();
    for (int f = 0; f < tetrahedra_meshes.mesh3Ds.size(); f++) {

        const vector<Vector4i>& surf = tetrahedra_meshes.mesh3Ds[f].surface;
        for (int j = 0; j < surf.size(); j++) {
            //Vector3d center = tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]] + tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]] + tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]];
            //center /= 3;
            //if (mouthOnly) {
            //    if (center.y() > 0.0 || center.y() < -0.35) continue;
            //    if (center.z() < 0.0) continue;
            //    if (center.x() > 0.23 || center.x() < -0.23) continue;
            //    if (tetrahedra_meshes.mesh3Ds[f].isJaw[surf[j][3]] || tetrahedra_meshes.mesh3Ds[f].isSkull[surf[j][3]]) {
            //        glColor3f(0.0f, 0.0f, 0.0f);
            //    }
            //    else {
            //        glColor3f(0.8f, 0.8f, 0.1f);
            //    }
            //}

            Vector3d v1 = tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]] - tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]];
            Vector3d v2 = tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]] - tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]];
            Vector3d normal = v1.cross(v2).normalized() * 0.01;

            //glVertex3f(center[0], center[1], center[2]);
            //glVertex3f(center[0] + normal[0], center[1] + normal[1], center[2] + normal[2]);
            //glVertex3f(xx, yy, zz);
            glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]][2]);
            glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]][2]);

            glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][1]][2]);
            glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]][2]);

            glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][0]][2]);
            glVertex3f(tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]][0], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]][1], tetrahedra_meshes.mesh3Ds[f].vertexes[surf[j][2]][2]);
            //}
        }
    }
    glEnd();
}


//void draw_Scene2D() {
//    glEnable(GL_DEPTH_TEST);
//    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   //sf ±³¾°ÑÕÉ«
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    draw_box2D(-1, -1, 2, 2);
//    draw_mesh2D();
//    glutSwapBuffers();
//}
int counttt = 0;
vector<float3> getRenderGeometry(int& number) {
    model_tet tetrahedra_meshes = simulator.getTetrahedraMeshes();
    mesh3D& meshTemp = tetrahedra_meshes.mesh3Ds[0];
    vector<Vector3d> meshNormal(meshTemp.vertexNum, Vector3d(0, 0, 0));
    number = meshTemp.surface.size();
    vector<float3> pos_normal_color(3 * number * 3);
#ifdef USE_TBB
    vector<tbb::spin_mutex> countMutex(meshTemp.vertexNum);
    tbb::parallel_for(0, number, 1, [&](int i)
#else
    for (int i = 0; i < number; i++)
#endif
    {
        int tetId = meshTemp.surface[i][3];
        int v0 = meshTemp.surface[i][0];
        int v1 = meshTemp.surface[i][1];
        int v2 = meshTemp.surface[i][2];
        Vector3d vt0 = meshTemp.vertexes[v0];// Vector3d(meshTemp.vertexes[v0][0], meshTemp.vertexes[v0][1], meshTemp.vertexes[v0][2]);
        Vector3d vt1 = meshTemp.vertexes[v1];// Vector3d(meshTemp.vertexes[v1][0], meshTemp.vertexes[v1][1], meshTemp.vertexes[v1][2]);
        Vector3d vt2 = meshTemp.vertexes[v2];// Vector3d(meshTemp.vertexes[v2][0], meshTemp.vertexes[v2][1], meshTemp.vertexes[v2][2]);
        Vector3d vec1 = vt1 - vt0;
        Vector3d vec2 = vt2 - vt0;
        Vector3d normal = vec1.cross(vec2).normalized();

        pos_normal_color[i * 9] = make_float3(vt0[0], vt0[1], vt0[2]);
        pos_normal_color[i * 9 + 3] = make_float3(vt1[0], vt1[1], vt1[2]);
        pos_normal_color[i * 9 + 6] = make_float3(vt2[0], vt2[1], vt2[2]);

            pos_normal_color[i * 9 + 2] = make_float3(0.6875f, 0.51953f, 0.38671f);
            pos_normal_color[i * 9 + 5] = make_float3(0.6875f, 0.51953f, 0.38671f);
            pos_normal_color[i * 9 + 8] = make_float3(0.6875f, 0.51953f, 0.38671f);

#ifdef USE_TBB
        countMutex[v0].lock();
        meshNormal[v0] += normal;
        countMutex[v0].unlock();
        countMutex[v1].lock();
        meshNormal[v1] += normal;
        countMutex[v1].unlock();
        countMutex[v2].lock();
        meshNormal[v2] += normal;
        countMutex[v2].unlock();
#else
        meshNormal[v0] += normal;
        meshNormal[v1] += normal;
        meshNormal[v2] += normal;
#endif
    }
#ifdef USE_TBB
    );
#endif

#ifdef USE_TBB
    tbb::parallel_for(0, number, 1, [&](int i)
#else
    for (int i = 0; i < number; i++)
#endif
    {
        int v0 = meshTemp.surface[i][0];
        int v1 = meshTemp.surface[i][1];
        int v2 = meshTemp.surface[i][2];
        //meshNormal[v0].normalize(); meshNormal[v1].normalize(); meshNormal[v2].normalize();
        pos_normal_color[i * 9 + 1] = make_float3(meshNormal[v0][0], meshNormal[v0][1], meshNormal[v0][2]);
        pos_normal_color[i * 9 + 4] = make_float3(meshNormal[v1][0], meshNormal[v1][1], meshNormal[v1][2]);
        pos_normal_color[i * 9 + 7] = make_float3(meshNormal[v2][0], meshNormal[v2][1], meshNormal[v2][2]);
    }
#ifdef USE_TBB
    );
#endif
    return pos_normal_color;
}


void draw_face_withShader() {


    int number;
    vector<float3> pos_normal_color = getRenderGeometry(number);

    glBindBuffer(GL_ARRAY_BUFFER, PN_vbo_);
    glBufferData(GL_ARRAY_BUFFER, 9 * number * sizeof(float3), &pos_normal_color[0], GL_DYNAMIC_DRAW);

    glBindVertexArray(VAO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float3), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float3), (GLvoid*)(sizeof(float3)));
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float3), (GLvoid*)(2 * sizeof(float3)));
    glEnableVertexAttribArray(2);
    //glBindVertexArray(0);

    glUseProgram(shaderProgram);
    //const glm::vec3 objectColor(0.9375f, 0.82031f, 0.78125f);
    const glm::vec3 objectColor(0.6875f, 0.51953f, 0.38671f);
    glUniform3fv(glGetUniformLocation(shaderProgram, "objectColor"), 1, &objectColor[0]);
    const glm::vec3 lightColor(1.0f, 1.0f, 1.0f);

    glUniform3fv(glGetUniformLocation(shaderProgram, "lightColor"), 1, &lightColor[0]);
    const glm::vec3 lightPos(0.0f, 0.5f, 3.5f);
    glUniform3fv(glGetUniformLocation(shaderProgram, "lightPos"), 1, &lightPos[0]);


    // create transformations
    glm::mat4 model = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 view = glm::mat4(1.0f);
    glm::mat4 projection = glm::mat4(1.0f);
    //model = glm::rotate(model, 45.f, glm::vec3(0.0f, 1.0f, 0.0f));
    view = glm::translate(view, glm::vec3(xTrans, yTrans, zTrans - 3));
    view = glm::rotate(view, xRot * 0.2f, glm::vec3(1.0f, 0.0f, 0.0f));
    view = glm::rotate(view, yRot * 0.2f, glm::vec3(0.0f, 1.0f, 0.0f));
    projection = glm::perspective(glm::radians(45.0f), (float)window_width / (float)window_height, 0.1f, 500.0f);
    // retrieve the matrix uniform locations
    unsigned int modelLoc = glGetUniformLocation(shaderProgram, "model");
    unsigned int viewLoc = glGetUniformLocation(shaderProgram, "view");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, &projection[0][0]);

    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, 3 * number);
}

void draw_Scene3D() {
    //face.mesh3Ds[0] = mesh3d;
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (!isSetShader) {
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glTranslatef(xTrans, yTrans, zTrans);
        glRotatef(xRot, 1.0f, 0.0f, 0.0f);
        glRotatef(yRot, 0.0f, 1.0f, 0.0f);

        draw_box3D(-1, -1, -1, 2, 2, 2);
        //draw_lines(-1, -1, -1, 2, 2, 2);
        if (true) {
            if (!drawSurface) {
                draw_tet_mesh3D();
            }
            else {
                draw_mesh3D();
                //draw_landMarks();
            }
        }
        glPopMatrix();
    }
    else {
        draw_face_withShader();
    }
    glutSwapBuffers();
    //glFlush();
}
double mfsum = 0;
double total_time = 0;
int total_cg_iterations = 0;
int total_newton_iterations = 0;
int start = -1;

void saveScreenPic(const string& path) {
    std::stringstream ss;
    ss << path;
    ss.fill('0');
    ss.width(5);
    ss << step;// / 10;
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
    ss << surfNumId++;// / 10;
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

void display(void)
{

    draw_Scene3D();

    if (saveSurface)
    {
        saveSurfaceMesh("OBJ/saveSurface/surf_");
    }


    if (stop) return;

    if (screenshot)
    {
        saveScreenPic("saveScreen/step_");
        step++;

    }
    int k = simulator.simulateStick();

}

void initScene0() {
    simulator.buildModels(0, 3);
}


void init(void)
{
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
    }
    std::cerr << "Status: Using GLEW " << glewGetString(GLEW_VERSION) << std::endl;
    glClearColor(0.0, 0.0, 0.0, 1.0);

    initScene0();

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
    //window_width = width;
    //window_height = height;

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



int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    //glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);

    glutSetOption(GLUT_MULTISAMPLE, 16);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);

    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("FEM");

    init();

    if (isSetShader) {
        set_shaders();
    }
    //glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
    //glEnable(GL_POINT_SPRITE_ARB);
    //glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
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
    //return 0;
}

