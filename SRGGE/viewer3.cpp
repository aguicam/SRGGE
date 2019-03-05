#include "viewer3.h"
#include "Helpers/mesh_importer.h"



//shader ID
enum {VIEW3};

Viewer3::Viewer3(const QGLFormat &glf, QWidget *parent) : HW(glf, parent) , mesh_importer(this)
{
    // init vars
    setFocusPolicy(Qt::StrongFocus);
}

void Viewer3::initializeGL()
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
    // Create the grid
    for(int k =0;k<N;k++){
        for(int i =0;i<N;i++){
            if(fmod(N,2)==1){
                grid.push_back((k+1)/2*cos(M_PI*k));
                grid.push_back(0);
                grid.push_back((i+1)/2*cos(M_PI*i));
            }
            else{
                grid.push_back(k/2*cos(M_PI*k)+0.5*cos(M_PI*k));
                grid.push_back(0);
                grid.push_back(i/2*cos(M_PI*i)+0.5*cos(M_PI*i));
            }

        }
    }
    //Initialize all positions to the lowest LOD
    for(int i=0;i<(int)grid.size()/3;i++){
        mod_Lod.push_back(0);
    }

//We load the Armadillo
    if(!mesh_importer.LoadModel(QString("models/Armadillo.ply")))
      QMessageBox::warning(this, tr("Error"),
                         tr("The file could not be opened"));
//We initialize  the hysteresis vector
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);

    //Save the number of faces and vertices
    num_face.push_back(mesh_->faces_16.size()/3);
    num_face.push_back(mesh_->faces_32.size()/3);
    num_face.push_back(mesh_->faces_64.size()/3);
    num_face.push_back(mesh_->faces_128.size()/3);
    num_face.push_back(mesh_->faces_.size()/3);

    num_vertices.push_back(mesh_->vertices_16.size()/3);
    num_vertices.push_back(mesh_->vertices_32.size()/3);
    num_vertices.push_back(mesh_->vertices_64.size()/3);
    num_vertices.push_back(mesh_->vertices_128.size()/3);
    num_vertices.push_back(mesh_->vertices_.size()/3);

    //Compute initial number of faces and vertices
    scene_triangles=N*N*mesh_->faces_[0];
    scene_vertices=N*N*num_vertices[0];
    mesh_importer.updateNumVert(scene_vertices,scene_triangles);
    //Compute the diagonal for the cost calculation
    diag=sqrt(pow(mesh_->max_[0]-mesh_->min_[0],2)+pow(mesh_->max_[1]-mesh_->min_[1],2)+pow(mesh_->max_[2]-mesh_->min_[2],2));



}

void Viewer3::resizeGL(int w, int h)
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

void Viewer3::paintGL()
{

     if(timer.isCreated()) timer.begin();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // use glsl program
    glUseProgram(m_program[VIEW3].programId());

    // pass the following parameters to vertex the shader:
    // projection matrix, modelview matrix, theta, and twist
    camera_.SetViewport();

    Eigen::Matrix4f projection = camera_.SetProjection();
    Eigen::Matrix4f view = camera_.SetView();
    Eigen::Matrix4f model = camera_.SetModel();


    const Eigen::Affine3f piRotation(Eigen::AngleAxisf(M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
    model = piRotation* model ;

    Eigen::Matrix4f t = view * model;
    Eigen::Matrix3f normal;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

    normal = normal.inverse().transpose();

    glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_PROJ], 1, GL_FALSE, projection.data());
    glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_MODEL], 1, GL_FALSE, model.data());
    glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_VIEW], 1, GL_FALSE, view.data());
    glUniformMatrix3fv(const_uniform[VIEW3][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
    glUniform1i(glGetUniformLocation(m_program[VIEW3].programId(), "colored"), (int)colored);

    // draw triangles
    if (mesh_ == nullptr) return;
    /*******************************************************/

    if(!max_throughput){  // If  we are under the maximum througput
        float low_cost=INT_MAX;
        int l_pos=0;
        //We look for the model with the lower cost to increase the LOD
        for(int i=0;i<(int)mod_Lod.size();i++){
            //We take into account only the models that can increase one LoD, and are not in the hysteresis vector
            if(mod_Lod[i]!=4&&(i!=hysteresis[0]&&i!=hysteresis[1]&&i!=hysteresis[2]&&i!=hysteresis[3]&&i!=hysteresis[4])){
                //We compute the cost as (2^(L+1)-2^L)*D/d
                Eigen::Vector4f view_c = view*(Eigen::Vector4f (grid[3*i],grid[3*i+1],grid[3*i+2],1));
                float D = sqrt(pow(view_c[0],2)+pow(view_c[1],2)+pow(view_c[2],2));
                float cost_diff= pow(2,mod_Lod[i]+1)*D/diag-pow(2,mod_Lod[i])*D/diag;
                // We save the position of he model that has the lowest cost
                if(std::min(low_cost,cost_diff)!=low_cost){
                    low_cost=cost_diff;
                    l_pos=i;
                }
             }

        }
        if(mod_Lod[l_pos]!=4){ // We increase the LoD of the model
         mod_Lod[l_pos]+=1;}
         hysteresis.push_back(l_pos); // We add the model to the hysteresis vector
         hysteresis.erase(hysteresis.begin());// We take out the first position of the hysteresis vector
        // We update count oftriangles and vertices of the scene
        scene_triangles= scene_triangles- num_face[mod_Lod[l_pos]-1]+num_face[mod_Lod[l_pos]];
        scene_vertices= scene_vertices- num_vertices[mod_Lod[l_pos]-1]+num_vertices[mod_Lod[l_pos]];
    mesh_importer.updateNumVert(scene_vertices,scene_triangles);
    }else{
        // If we had maximum throuput, we perform the same opperations but looking for the maximum cost to decrease one LoD
        float high_cost=INT_MIN;
        int h_pos=0;
        for(int i=0;i<(int)mod_Lod.size();i++){
             //We take into account only the models that can decrease one LoD, and are not in the hysteresis vector
            if(mod_Lod[i]!=0&&(i!=hysteresis[0]&&i!=hysteresis[1]&&i!=hysteresis[2]&&i!=hysteresis[3]&&i!=hysteresis[4])){
                Eigen::Vector4f view_c = view*(Eigen::Vector4f (grid[3*i],grid[3*i+1],grid[3*i+2],1));//model*
                float D = sqrt(pow(view_c[0],2)+pow(view_c[1],2)+pow(view_c[2],2));
                float cost_diff= pow(2,mod_Lod[i]+1)*D/diag-pow(2,mod_Lod[i])*D/diag;
                //
                if(std::max(high_cost,cost_diff)!=high_cost){
                    high_cost=cost_diff;
                    h_pos=i;
                }
            }

        }

        mod_Lod[h_pos]-=1;// We decrease the LoD of the model
        hysteresis.push_back(h_pos);// We add the model to the hysteresis vector
        hysteresis.erase(hysteresis.begin());// We take out the first position of the hysteresis vector
        // We update count oftriangles and vertices of the scene
        scene_triangles= scene_triangles- num_face[mod_Lod[h_pos]+1]+num_face[mod_Lod[h_pos]];
        scene_vertices= scene_vertices- num_vertices[mod_Lod[h_pos]+1]+num_vertices[mod_Lod[h_pos]];
    mesh_importer.updateNumVert(scene_vertices,scene_triangles);


    }
    /*******************************************************/
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");
    Eigen::Matrix4f translation;
        /*** We draw LOD 0 ***/
        glBindVertexArray(VAO_ID_16);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_16);
        for(int l=0;l<(int)mod_Lod.size();l++){
            if(mod_Lod[l]==0){
                Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(grid[3*l],grid[3*l+1],grid[3*l+2])));
                translation = (kTran.matrix()*model);
                glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
                glUniform1i(glGetUniformLocation(m_program[VIEW3].programId(), "lod"), 0);
                glDrawElements(GL_TRIANGLES,mesh_->faces_16.size(), GL_UNSIGNED_INT, 0);
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);

        /*** We draw LOD 1 ***/
        glBindVertexArray(VAO_ID_32);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_32);
        for(int l=0;l<(int)mod_Lod.size();l++){
            if(mod_Lod[l]==1){
                Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(grid[3*l],grid[3*l+1],grid[3*l+2])));
                translation = (kTran.matrix()*model);
                glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
                glUniform1i(glGetUniformLocation(m_program[VIEW3].programId(), "lod"), 1);
                glDrawElements(GL_TRIANGLES,mesh_->faces_32.size(), GL_UNSIGNED_INT, 0);
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);


        /*** We draw LOD 2 ***/
        glBindVertexArray(VAO_ID_64);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_64);
        for(int l=0;l<(int)mod_Lod.size();l++){
            if(mod_Lod[l]==2){
                Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(grid[3*l],grid[3*l+1],grid[3*l+2])));
                translation = (kTran.matrix()*model);
                glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
                glUniform1i(glGetUniformLocation(m_program[VIEW3].programId(), "lod"), 2);
                glDrawElements(GL_TRIANGLES,mesh_->faces_64.size(), GL_UNSIGNED_INT, 0);
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);


        /*** We draw LOD 3 ***/
        glBindVertexArray(VAO_ID_128);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_128);
        for(int l=0;l<(int)mod_Lod.size();l++){
            if(mod_Lod[l]==3){
                Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(grid[3*l],grid[3*l+1],grid[3*l+2])));
                translation = (kTran.matrix()*model);
                glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
                glUniform1i(glGetUniformLocation(m_program[VIEW3].programId(), "lod"), 3);
                glDrawElements(GL_TRIANGLES,mesh_->faces_128.size(), GL_UNSIGNED_INT, 0);
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);


        /*** We draw LOD 4 ***/
        glBindVertexArray(VAO_ID);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID);
        for(int l=0;l<(int)mod_Lod.size();l++){
            if(mod_Lod[l]==4){
                Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(grid[3*l],grid[3*l+1],grid[3*l+2])));
                translation = (kTran.matrix()*model);
                glUniformMatrix4fv(const_uniform[VIEW3][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
                glUniform1i(glGetUniformLocation(m_program[VIEW3].programId(), "lod"), 4);
                glDrawElements(GL_TRIANGLES,mesh_->faces_.size(), GL_UNSIGNED_INT, 0);
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        glBindVertexArray(0);


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
    //We check the state of the maximmum throughput
    if(fps<32){max_throughput=true;}
    if(fps>=32){max_throughput=false;}
//    if(scene_triangles>16000000){max_throughput=true;}
//    if(scene_triangles<=16000000){max_throughput=false;}

}

QGroupBox* Viewer3::controlPanel()
{
    // init group box
    QGroupBox *groupBox = HW::controlPanel();
    groupBox->setStyleSheet(GroupBoxStyle);
    mesh_importer.AddImportExportPLY(groupBox);
    implementNxN(groupBox);
    implementColored(groupBox);
    showFramerate(groupBox);
    //HW::AddImportExportPLY(groupBox);

    return(groupBox);
}

void Viewer3::reset()
{
    mesh_ = nullptr;	// recompute geometry
    initVertexBuffer();
    // draw
    updateGL();
}

void Viewer3::initShaders()
{
    // init uniforms hash table based on uniform variable names and location IDs
    UniformMap uniforms;

    // compile shader, bind attribute vars, link shader, and initialize uniform var table
    initShader(VIEW3, QString(":/SRGGE/Shaders/s4s5s6.vert"), QString(":/SRGGE/Shaders/s4s5s6.frag"), uniforms);
}

void Viewer3::initVertexBuffer()
{

    if (mesh_ == nullptr) return;

    initBuffer_original();
    initBuffer_16();
    initBuffer_32();
    initBuffer_64();
    initBuffer_128();

   updateGL();


}


void Viewer3::importModel(){
    QString filename;

      filename = QFileDialog::getOpenFileName(this, tr("Load model"), "./",
                                              tr("PLY Files ( *.ply )"));
      if (!filename.isNull()) {
          if(!mesh_importer.LoadModel(filename))
            QMessageBox::warning(this, tr("Error"),
                               tr("The file could not be opened"));
      }

}

void Viewer3::exportModel(){
    QString filename;

      filename = QFileDialog::getSaveFileName(this, tr("Export model"), "untitled.ply",
                                              tr("PLY Files ( *.ply )"));
      if (!filename.isNull()) {
          if(!mesh_importer.SaveModel(filename,LoD,*mesh_))
            QMessageBox::warning(this, tr("Error"),
                               tr("The file could not be saved"));
      }

}

void Viewer3::showFramerate(QGroupBox *groupBox){
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

void Viewer3::implementNxN(QGroupBox*groupBox){
    QPushButton *buttonOK  = new QPushButton("OK");
    echoLineEdit->setPlaceholderText("Introduce N");
    echoLineEdit->setFocus();

    QFrame* line = new QFrame();
    line->setFrameShape(QFrame::HLine);
    auto layout = dynamic_cast<QGridLayout*>(groupBox->layout());

    int row = layout->rowCount() + 1;

    QLabel *titleB = new QLabel(" NxN ");

    layout->addWidget(line,row,0, 1 ,layout->columnCount());
    row++;
    layout->addWidget(titleB,row,0);
    row++;
    layout->addWidget(echoLineEdit,row,0);
    layout->addWidget(buttonOK,row,1);
    HW::connect(buttonOK, SIGNAL(clicked()), this,SLOT(getLineText()));
}
void Viewer3::implementColored(QGroupBox*groupBox){
    QCheckBox *boxColor  = new QCheckBox("Colored LoDs:");
    boxColor->setChecked(true);
    colored =true;

    QFrame* line = new QFrame();
    line->setFrameShape(QFrame::HLine);
    auto layout = dynamic_cast<QGridLayout*>(groupBox->layout());

    int row = layout->rowCount() + 1;

    layout->addWidget(line,row,0, 1 ,layout->columnCount());
    row++;
    layout->addWidget(boxColor ,row,0);
    HW::connect(boxColor, SIGNAL(stateChanged(int)), this,SLOT(changeColored(int)));
}
void Viewer3::changeColored(int c){
    colored = c;
    update();
}

void Viewer3::getLineText(){
    // We get the line
    N=(echoLineEdit->text()).toFloat();
    grid.clear();

    //We build the new grid
    for(int k =0;k<N;k++){
        for(int i =0;i<N;i++){
            if(fmod(N,2)==1){
                grid.push_back((k+1)/2*cos(M_PI*k));
                grid.push_back(0);
                grid.push_back((i+1)/2*cos(M_PI*i));
            }
            else{
                grid.push_back(k/2*cos(M_PI*k)+0.5*cos(M_PI*k));
                grid.push_back(0);
                grid.push_back(i/2*cos(M_PI*i)+0.5*cos(M_PI*i));
            }

        }
    }
    //Initiallize all models LODs
    mod_Lod.clear();
    for(int i=0;i<(int)grid.size()/3;i++){
        mod_Lod.push_back(0);
    }
    //Initiallize hysteresis vector
    hysteresis.clear();
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
//Update num of faces and vertices
    scene_triangles=N*N*num_face[0];
    scene_vertices=N*N*num_vertices[0];
    mesh_importer.updateNumVert(scene_vertices,scene_triangles);
    initVertexBuffer();
}


//***********INIT VERTEX FOR ALL LEVEL OF DETAIL***************//
void Viewer3::initBuffer_original(){

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
void Viewer3::initBuffer_16(){
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
void Viewer3::initBuffer_32(){
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
void Viewer3::initBuffer_64(){
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
void Viewer3::initBuffer_128(){
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

