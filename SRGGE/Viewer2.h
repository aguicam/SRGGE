#ifndef VIEWER2_H
#define VIEWER2_H

#include "HW.h"
#include "Helpers/mesh_importer.h"
#include <iostream>
#include <iomanip>

class Viewer2 : public HW
{
    Q_OBJECT
public:
    Viewer2		(const QGLFormat &glf, QWidget *parent = 0);
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
    void        chooseLOD(QGroupBox *groupBox);

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
    QLineEdit *echoLineEdit = new QLineEdit;
    float N=1;
    QLabel *m_fpsCount;
    QString fpsString;

    QOpenGLTimerQuery timer;
    GLuint64 resTime;
    float fps;

    int LoD;


public slots:
    void importModel ();
    void exportModel();
    void setLOD(int b);
};

#endif // VIEWER1_H
