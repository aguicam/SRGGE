#ifndef VIEWER3_H
#define VIEWER3_H



#include "HW.h"
#include "Helpers/mesh_importer.h"
#include <iostream>
#include <iomanip>
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
#include <QCheckBox>

class Viewer3 : public HW
{
    Q_OBJECT
public:
    Viewer3		(const QGLFormat &glf, QWidget *parent = 0);
    QGroupBox*	controlPanel	();		// create control panel
    void		reset		();		// reset parameters
    void		initVertexBuffer();		// init vertices
    void        initBuffer_original();
    void        initBuffer_16();
    void        initBuffer_32();
    void        initBuffer_64();
    void        initBuffer_128();
    void		initShaders	();		// init shaders
    void        showFramerate(QGroupBox *groupBox);
    void        implementNxN(QGroupBox*groupBox);
    void        implementColored(QGroupBox*groupBox);


protected:
    void		initializeGL	();		// init GL state
    void		resizeGL	(int, int);	// resize GL widget
    void		paintGL		();		// render GL scene
    void		divideTriangle(vec2, vec2, vec2, int);	// subdivide triangle
    void		triangle(vec3, vec3, vec3);	// process single triangle

private:
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
    std::vector<int> num_face;
    std::vector<int> num_vertices;
    bool colored;

    //Grid
    std::vector<float> grid;
    int threshold;
    float N=9;
    long int scene_triangles;
    long int scene_vertices;
    float diag;
    bool max_throughput=false;
    std::vector<int> mod_Lod;
    std::vector<int> hysteresis;
    QLineEdit *echoLineEdit = new QLineEdit;



public slots:
    void importModel ();
    void exportModel();
    void getLineText();
    void changeColored(int c);
};

#endif // VIEWER3_H

