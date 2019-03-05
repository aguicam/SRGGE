#include "Viewer2.h"
#include "Helpers/mesh_importer.h"
#include <QOpenGLFunctions>
#include <QtOpenGL>
#include <vector>
#include <Eigen>
#include <QOpenGLFunctions>
#include <QtOpenGL>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <math.h>
#include <QRadioButton>
#include <QLineEdit>



//shader ID
enum {VIEW2};

Viewer2::Viewer2(const QGLFormat &glf, QWidget *parent) : HW(glf, parent) , mesh_importer(this)
{
    // init vars
    setFocusPolicy(Qt::StrongFocus);
}

void Viewer2::initializeGL()
{
    // initialize GL function resolution for current context
    initializeGLFunctions();

    // init vertex and fragment shaders
    initShaders();

    glEnable(GL_NORMALIZE);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);


    // initialize vertex buffer and write positions to vertex shader
    initVertexBuffer();

    // init state variables
    glClearColor(1.0, 1.0, 1.0, 1.0);	// set background color
    glColor3f   (1.0, 1.0, 0.0);		// set foreground color

    if(!timer.isCreated()){
        timer.create();
    }

    LoD=4;

}

void Viewer2::resizeGL(int w, int h)
{

    if (h == 0) h = 1;
    width_ = w;
    height_ = h;

    camera_.SetViewport(0, 0, w, h);
    camera_.SetProjection(kFieldOfView, kZNear, kZFar);

    m_winW = w;
    m_winH = h;

    return;
}

void Viewer2::paintGL()
{

     if(timer.isCreated()) timer.begin();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // use glsl program
    glUseProgram(m_program[VIEW2].programId());

    // pass the following parameters to vertex the shader:
    // projection matrix, modelview matrix, theta, and twist
    camera_.SetViewport();

    Eigen::Matrix4f projection = camera_.SetProjection();
    Eigen::Matrix4f view = camera_.SetView();
    Eigen::Matrix4f model = camera_.SetModel();

    Eigen::Matrix4f t = view * model;
    Eigen::Matrix3f normal;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

    normal = normal.inverse().transpose();

    glUniformMatrix4fv(const_uniform[VIEW2][CONST_UNIFORM_PROJ], 1, GL_FALSE, projection.data());
    glUniformMatrix4fv(const_uniform[VIEW2][CONST_UNIFORM_MODEL], 1, GL_FALSE, model.data());
    glUniformMatrix4fv(const_uniform[VIEW2][CONST_UNIFORM_VIEW], 1, GL_FALSE, view.data());
    glUniformMatrix3fv(const_uniform[VIEW2][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());

    const Eigen::Affine3f piRotation(Eigen::AngleAxisf(M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
    model = piRotation* model ;
    // draw triangles
    if (mesh_ == nullptr) return;

    //////////
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");
    switch(LoD){

    case 0:
        glBindVertexArray(VAO_ID_16);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_16);
        glDrawElements(GL_TRIANGLES,mesh_->faces_16.size(), GL_UNSIGNED_INT, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);
        break;
    case 1:
        glBindVertexArray(VAO_ID_32);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_32);
        glDrawElements(GL_TRIANGLES,mesh_->faces_32.size(), GL_UNSIGNED_INT, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);
        break;
    case 2:
        glBindVertexArray(VAO_ID_64);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_64);
        glDrawElements(GL_TRIANGLES,mesh_->faces_64.size(), GL_UNSIGNED_INT, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);
        break;
    case 3:
        glBindVertexArray(VAO_ID_128);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_128);
        glDrawElements(GL_TRIANGLES,mesh_->faces_128.size(), GL_UNSIGNED_INT, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);
        break;
    case 4:
        glBindVertexArray(VAO_ID);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID);
        glDrawElements(GL_TRIANGLES,mesh_->faces_.size(), GL_UNSIGNED_INT, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);
        break;

    }





    // terminate program; rendering is done
    glUseProgram(0);


    if(timer.isCreated()){
         timer.end();
         resTime = timer.waitForResult();
         if(resTime==0){return;}

         fps = 1.0/((float)resTime/1000000000.0);

         std::stringstream s_fps;
         s_fps << std::fixed << std::setprecision(1);
         s_fps << fps << " FPS" << std::ends;
         s_fps << std::resetiosflags(std::ios_base::fixed | std::ios_base::floatfield);
         fpsString = QString::fromStdString(s_fps.str());
         m_fpsCount->setText(fpsString);
     }




}

QGroupBox* Viewer2::controlPanel()
{
    // init group box
    QGroupBox *groupBox = HW::controlPanel();
    groupBox->setStyleSheet(GroupBoxStyle);
    mesh_importer.AddImportExportPLY(groupBox);
    chooseLOD(groupBox);
    showFramerate(groupBox);
    //HW::AddImportExportPLY(groupBox);

    return(groupBox);
}

void Viewer2::reset()
{
    mesh_ = nullptr;	// recompute geometry
    initVertexBuffer();
    // draw
    updateGL();
}

void Viewer2::initShaders()
{
    // init uniforms hash table based on uniform variable names and location IDs
    UniformMap uniforms;

    // compile shader, bind attribute vars, link shader, and initialize uniform var table
    initShader(VIEW2, QString(":/SRGGE/Shaders/s1s2s3.vert"), QString(":/SRGGE/Shaders/s1s2s3.frag"), uniforms);
}

void Viewer2::initVertexBuffer()
{

    if (mesh_ == nullptr) return;

    initBuffer_original();
    initBuffer_16();
    initBuffer_32();
    initBuffer_64();
    initBuffer_128();

   updateGL();


}


void Viewer2::importModel(){
    QString filename;

      filename = QFileDialog::getOpenFileName(this, tr("Load model"), "./",
                                              tr("PLY Files ( *.ply )"));
      if (!filename.isNull()) {
          if(!mesh_importer.LoadModel(filename))
            QMessageBox::warning(this, tr("Error"),
                               tr("The file could not be opened"));
      }

}
void Viewer2::exportModel(){
    QString filename;

      filename = QFileDialog::getSaveFileName(this, tr("Export model"), "untitled.ply",
                                              tr("PLY Files ( *.ply )"));
      if (!filename.isNull()) {
          if(!mesh_importer.SaveModel(filename,LoD,*mesh_))
            QMessageBox::warning(this, tr("Error"),
                               tr("The file could not be saved"));
      }

}


void Viewer2::showFramerate(QGroupBox *groupBox){
    m_fpsCount = new QLabel("0");
    QFrame* line = new QFrame();

    line->setFrameShape(QFrame::HLine);
    QLabel *title = new QLabel(" FPS:");

    auto layout = dynamic_cast<QGridLayout*>(groupBox->layout());

    int row = layout->rowCount() + 1;
    int headRow = row;
    row++;
    layout->addWidget(title,row,0);
 //   row++;
    layout->addWidget(m_fpsCount,row, 1);

    layout->addWidget(line,headRow, 0, 1 ,layout->columnCount());

    groupBox->setLayout(layout);
}

void Viewer2::chooseLOD(QGroupBox *groupBox){

    QLabel *levelDetail= new QLabel(tr("Level of detail:"));

    QSpinBox *levelDetailBox = new QSpinBox;
    levelDetailBox->setRange(0, 4);
    levelDetailBox->setSingleStep(1);
    levelDetailBox->setValue(4);

    QFrame* line = new QFrame();
    line->setFrameShape(QFrame::HLine);
    QLabel *title = new QLabel("Vertex Clustering:");


    connect(levelDetailBox, SIGNAL(valueChanged(int)),this, SLOT(setLOD(int)));
    auto layout = dynamic_cast<QGridLayout*>(groupBox->layout());

    int row = layout->rowCount() + 1;
    int headRow = row;
    row++;
    layout->addWidget(title,row,0);
    row++;
    layout->addWidget(levelDetail,row, 0);
    layout->addWidget(levelDetailBox,row, 1);

    layout->addWidget(line,headRow, 0, 1 ,layout->columnCount());

    groupBox->setLayout(layout);

}
void Viewer2::setLOD(int b){
    LoD=b;

    switch(LoD){
    case 0:
        mesh_importer.updateNumVert(mesh_->vertices_16.size()/3,(int)mesh_->faces_16.size() / 3);
        break;
    case 1:
        mesh_importer.updateNumVert(mesh_->vertices_32.size()/3,(int)mesh_->faces_32.size() / 3);
        break;
    case 2:
        mesh_importer.updateNumVert(mesh_->vertices_64.size()/3,(int)mesh_->faces_64.size() / 3);
        break;
    case 3:
        mesh_importer.updateNumVert(mesh_->vertices_128.size()/3,(int)mesh_->faces_128.size() / 3);
        break;
    case 4:
        mesh_importer.updateNumVert((int)mesh_->vertices_.size() / 3,(int)mesh_->faces_.size() / 3);
        break;
    }

    //mesh_importer.updateNumVert((int) N*N*mesh_->vertices_.size() / 3,(int)N*N*mesh_->faces_.size() / 3);
    initVertexBuffer();
}


//***********INIT VERTEX FOR ALL LEVEL OF DETAIL***************//
void Viewer2::initBuffer_original(){

    glGenBuffers(1, &vertexBufferVBO_ID);
    glGenBuffers(1, &normalBufferVBO_ID);
    glGenBuffers(1, &facesBufferVBO_ID);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_->vertices_.size()/3;i++){

        data.push_back(mesh_->vertices_[3*i]);
        data.push_back(mesh_->vertices_[3*i+1]);
        data.push_back(mesh_->vertices_[3*i+2]);
        data.push_back(mesh_->normals_[3*i]);
        data.push_back(mesh_->normals_[3*i+1]);
        data.push_back(mesh_->normals_[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID);

    glBindVertexArray(VAO_ID);
    glGenBuffers(1, &VBO_ID);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_->vertices_.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_->faces_.size() * sizeof(int), &mesh_->faces_[0], GL_STATIC_DRAW);
     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
     glBindVertexArray(0);
}
void Viewer2::initBuffer_16(){
    glGenBuffers(1, &vertexBufferVBO_ID_16);
    glGenBuffers(1, &normalBufferVBO_ID_16);
    glGenBuffers(1, &facesBufferVBO_ID_16);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_->vertices_16.size()/3;i++){

        data.push_back(mesh_->vertices_16[3*i]);
        data.push_back(mesh_->vertices_16[3*i+1]);
        data.push_back(mesh_->vertices_16[3*i+2]);
        data.push_back(mesh_->normals_16[3*i]);
        data.push_back(mesh_->normals_16[3*i+1]);
        data.push_back(mesh_->normals_16[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID_16);

    glBindVertexArray(VAO_ID_16);
    glGenBuffers(1, &VBO_ID_16);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID_16);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_->vertices_16.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_16);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_->faces_16.size() * sizeof(int), &mesh_->faces_16[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer2::initBuffer_32(){
    glGenBuffers(1, &vertexBufferVBO_ID_32);
    glGenBuffers(1, &normalBufferVBO_ID_32);
    glGenBuffers(1, &facesBufferVBO_ID_32);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_->vertices_32.size()/3;i++){

        data.push_back(mesh_->vertices_32[3*i]);
        data.push_back(mesh_->vertices_32[3*i+1]);
        data.push_back(mesh_->vertices_32[3*i+2]);
        data.push_back(mesh_->normals_32[3*i]);
        data.push_back(mesh_->normals_32[3*i+1]);
        data.push_back(mesh_->normals_32[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID_32);

    glBindVertexArray(VAO_ID_32);
    glGenBuffers(1, &VBO_ID_32);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID_32);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_->vertices_32.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_32);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_->faces_32.size() * sizeof(int), &mesh_->faces_32[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer2::initBuffer_64(){
    glGenBuffers(1, &vertexBufferVBO_ID_64);
    glGenBuffers(1, &normalBufferVBO_ID_64);
    glGenBuffers(1, &facesBufferVBO_ID_64);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_->vertices_64.size()/3;i++){

        data.push_back(mesh_->vertices_64[3*i]);
        data.push_back(mesh_->vertices_64[3*i+1]);
        data.push_back(mesh_->vertices_64[3*i+2]);
        data.push_back(mesh_->normals_64[3*i]);
        data.push_back(mesh_->normals_64[3*i+1]);
        data.push_back(mesh_->normals_64[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID_64);

    glBindVertexArray(VAO_ID_64);
    glGenBuffers(1, &VBO_ID_64);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID_64);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_->vertices_64.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_64);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_->faces_64.size() * sizeof(int), &mesh_->faces_64[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer2::initBuffer_128(){
    glGenBuffers(1, &vertexBufferVBO_ID_128);
    glGenBuffers(1, &normalBufferVBO_ID_128);
    glGenBuffers(1, &facesBufferVBO_ID_128);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_->vertices_128.size()/3;i++){

        data.push_back(mesh_->vertices_128[3*i]);
        data.push_back(mesh_->vertices_128[3*i+1]);
        data.push_back(mesh_->vertices_128[3*i+2]);
        data.push_back(mesh_->normals_128[3*i]);
        data.push_back(mesh_->normals_128[3*i+1]);
        data.push_back(mesh_->normals_128[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID_128);

    glBindVertexArray(VAO_ID_128);
    glGenBuffers(1, &VBO_ID_128);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID_128);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_->vertices_128.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_128);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_->faces_128.size() * sizeof(int), &mesh_->faces_128[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}

