#ifndef VIEWER4_H
#define VIEWER4_H

#include "HW.h"
#include "Helpers/mesh_importer.h"
#include "Helpers/vis.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <QOpenGLFunctions>
#include <QtOpenGL>
#include <vector>
#include <Eigen>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <math.h>
#include <QRadioButton>
#include <QLineEdit>
#include <Helpers/stb_image.h>
#include <QOpenGLTexture>

class Viewer4: public HW
{
    Q_OBJECT
public:
    Viewer4		(const QGLFormat &glf, QWidget *parent = 0);
    QGroupBox*	controlPanel	();		// create control panel
    void		reset		();		// reset parameters
    void		initVertexBuffer();		// init vertices
    void        initBuffer_original();
    void        initBuffer_16();
    void        initBuffer_32();
    void        initBuffer_64();
    void        initBuffer_128();
    void        initBuffer1_original();
    void        initBuffer1_16();
    void        initBuffer1_32();
    void        initBuffer1_64();
    void        initBuffer1_128();
    void        initBuffer2_original();
    void        initBuffer2_16();
    void        initBuffer2_32();
    void        initBuffer2_64();
    void        initBuffer2_128();
    void        initBuffer_WALL();
    void        initBuffer_FLOOR();
    void		initShaders	();		// init shaders
    void        showFramerate(QGroupBox *groupBox);
    void        loadPositions(int type);
    void        draw_Armadillo(Eigen::Matrix4f model,Eigen::Matrix4f view);
    void        draw_Bunny(Eigen::Matrix4f model,Eigen::Matrix4f view);
    void        draw_Anno(Eigen::Matrix4f model,Eigen::Matrix4f view);
    void        LoadDic();
    void        implementColored(QGroupBox*groupBox);
protected:
    void		initializeGL	();		// init GL state
    void		resizeGL	(int, int);	// resize GL widget
    void		paintGL		();		// render GL scene
    void		divideTriangle(vec2, vec2, vec2, int);	// subdivide triangle
    void		triangle(vec3, vec3, vec3);	// process single triangle

private:

    bool colored;
    //Model 0 -> Armadillo
    GLuint        vertexBufferVBO_ID;
    GLuint        normalBufferVBO_ID;
    GLuint        facesBufferVBO_ID;
    GLuint        VAO_ID;
    GLuint        VBO_ID;

    GLuint        vertexBufferVBO_ID_16;
    GLuint        normalBufferVBO_ID_16;
    GLuint        facesBufferVBO_ID_16;
    GLuint        VAO_ID_16;
    GLuint        VBO_ID_16;

    GLuint        vertexBufferVBO_ID_32;
    GLuint        normalBufferVBO_ID_32;
    GLuint        facesBufferVBO_ID_32;
    GLuint        VAO_ID_32;
    GLuint        VBO_ID_32;

    GLuint        vertexBufferVBO_ID_64;
    GLuint        normalBufferVBO_ID_64;
    GLuint        facesBufferVBO_ID_64;
    GLuint        VAO_ID_64;
    GLuint        VBO_ID_64;

    GLuint        vertexBufferVBO_ID_128;
    GLuint        normalBufferVBO_ID_128;
    GLuint        facesBufferVBO_ID_128;
    GLuint        VAO_ID_128;
    GLuint        VBO_ID_128;

    float scaleArm;
    float transYArm;


    //Model 1 -> Bunny
    GLuint        vertexBufferVBO_ID1;
    GLuint        normalBufferVBO_ID1;
    GLuint        facesBufferVBO_ID1;
    GLuint        VAO_ID1;
    GLuint        VBO_ID1;

    GLuint        vertexBufferVBO_ID1_16;
    GLuint        normalBufferVBO_ID1_16;
    GLuint        facesBufferVBO_ID1_16;
    GLuint        VAO_ID1_16;
    GLuint        VBO_ID1_16;

    GLuint        vertexBufferVBO_ID1_32;
    GLuint        normalBufferVBO_ID1_32;
    GLuint        facesBufferVBO_ID1_32;
    GLuint        VAO_ID1_32;
    GLuint        VBO_ID1_32;

    GLuint        vertexBufferVBO_ID1_64;
    GLuint        normalBufferVBO_ID1_64;
    GLuint        facesBufferVBO_ID1_64;
    GLuint        VAO_ID1_64;
    GLuint        VBO_ID1_64;

    GLuint        vertexBufferVBO_ID1_128;
    GLuint        normalBufferVBO_ID1_128;
    GLuint        facesBufferVBO_ID1_128;
    GLuint        VAO_ID1_128;
    GLuint        VBO_ID1_128;

    float scaleBun;
    float transYBun;

    //Model 2 -> Anno
    GLuint        vertexBufferVBO_ID2;
    GLuint        normalBufferVBO_ID2;
    GLuint        facesBufferVBO_ID2;
    GLuint        VAO_ID2;
    GLuint        VBO_ID2;

    GLuint        vertexBufferVBO_ID2_16;
    GLuint        normalBufferVBO_ID2_16;
    GLuint        facesBufferVBO_ID2_16;
    GLuint        VAO_ID2_16;
    GLuint        VBO_ID2_16;

    GLuint        vertexBufferVBO_ID2_32;
    GLuint        normalBufferVBO_ID2_32;
    GLuint        facesBufferVBO_ID2_32;
    GLuint        VAO_ID2_32;
    GLuint        VBO_ID2_32;

    GLuint        vertexBufferVBO_ID2_64;
    GLuint        normalBufferVBO_ID2_64;
    GLuint        facesBufferVBO_ID2_64;
    GLuint        VAO_ID2_64;
    GLuint        VBO_ID2_64;

    GLuint        vertexBufferVBO_ID2_128;
    GLuint        normalBufferVBO_ID2_128;
    GLuint        facesBufferVBO_ID2_128;
    GLuint        VAO_ID2_128;
    GLuint        VBO_ID2_128;

    float scaleAnn;
    float transYAnn;

    // Walls

    GLuint        vertexBufferVBO_WALL;
    GLuint        normalBufferVBO_WALL;
    GLuint        facesBufferVBO_WALL;
    GLuint        VAO_WALL;
    GLuint        VBO_WALL;
    QOpenGLTexture        *wall_tex;
    int         size_faces_wall;

    // Floor
    GLuint        vertexBufferVBO_FLOOR;
    GLuint        normalBufferVBO_FLOOR;
    GLuint        facesBufferVBO_FLOOR;
    GLuint        VAO_FLOOR;
    GLuint        VBO_FLOOR;
    QOpenGLTexture        *floor_tex;
    int         size_faces_floor;


    typedef void (APIENTRY *_glGenVertexArrays) (GLsizei, GLuint*);
    typedef void (APIENTRY *_glBindVertexArray) (GLuint);
    _glGenVertexArrays glGenVertexArrays;
    _glBindVertexArray glBindVertexArray;
    MI            mesh_importer;
    QLabel *m_fpsCount;
    QString fpsString;

    QOpenGLTimerQuery timer;
    GLuint64 resTime;
    float fps;

    int LoD;

    //Grid
    int threshold;
    float N=61;
    long int scene_triangles;
    long int scene_vertices;
    std::vector <float> diag;
    bool max_throughput=false;
    std::vector<int> hysteresis;

    struct models{
        int cellID;
        int x;
        int z;
        float rotation;
        int type;
        int lod;
    };

    std::vector<models> objects;
    std::vector<int> dic;




public slots:
    void importModel ();
    void changeColored(int c);
};

#endif // VIEWER4_H
