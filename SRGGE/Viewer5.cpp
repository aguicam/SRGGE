#include "Viewer5.h"
#include "Helpers/mesh_importer.h"



//shader ID
enum {VIEW5,WALL,FLOOR};

Viewer5::Viewer5(const QGLFormat &glf, QWidget *parent) : HW(glf, parent) , mesh_importer(this)
{
    // init vars
    setFocusPolicy(Qt::StrongFocus);
}

void Viewer5::initializeGL()
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
    //initVertexBuffer();

    // init state variables
    glClearColor(1.0, 1.0, 1.0, 1.0);	// set background color
    glColor3f   (1.0, 1.0, 0.0);		// set foreground color

    if(!timer.isCreated()){
        timer.create();
    }

    /************VISIBILITY*****************/
//    Vis visib = Vis();
//    //visib.CalculateVisibility(); //BASIC
//    visib.CalculateBresenham();//ADVANCED
//    visib.makeUnique();
//    visib.writeToFile();
    /****************************************/

    LoD=4;

    /*******We load the 3 models*********/
    std::vector<QString> strM;
    strM.push_back(QString("models/Simplified/Armadillo0.ply"));
    strM.push_back(QString("models/Simplified/Armadillo1.ply"));
    strM.push_back(QString("models/Simplified/Armadillo2.ply"));
    strM.push_back(QString("models/Simplified/Armadillo3.ply"));
    strM.push_back(QString("models/Armadillo.ply"));
    if(!mesh_importer.Load3Models(strM,0))
      QMessageBox::warning(this, tr("Error"),
                         tr("The file could not be opened"));
    strM.clear();
    strM.push_back(QString("models/Simplified/bunny0.ply"));
    strM.push_back(QString("models/Simplified/bunny1.ply"));
    strM.push_back(QString("models/Simplified/bunny2.ply"));
    strM.push_back(QString("models/Simplified/bunny3.ply"));
    strM.push_back(QString("models/bunny.ply"));

    if(!mesh_importer.Load3Models(strM,1))
      QMessageBox::warning(this, tr("Error"),
                         tr("The file could not be opened"));
    strM.clear();
    strM.push_back(QString("models/Simplified/Anno0.ply"));
    strM.push_back(QString("models/Simplified/Anno1.ply"));
    strM.push_back(QString("models/Simplified/Anno2.ply"));
    strM.push_back(QString("models/Simplified/Anno3.ply"));
    strM.push_back(QString("models/Anno.ply"));

    if(!mesh_importer.Load3Models(strM,2))
      QMessageBox::warning(this, tr("Error"),
                         tr("The file could not be opened"));


    /*****We are going to read the positions of each model*********/
    objects.clear();
//    //ARMADILLO
    loadPositions(0);
//    //BUNNY
    loadPositions(1);
//    //ANNO
    loadPositions(2);
    //WALL
    loadPositions(3);

    loadVisibility(); //We load the visibility information
    LoadDic();

    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);
    hysteresis.push_back(-1);


      diag.clear();
      diag.push_back(sqrt(pow(mesh_->max_[0]-mesh_->min_[0],2)+pow(mesh_->max_[1]-mesh_->min_[1],2)+pow(mesh_->max_[2]-mesh_->min_[2],2)));
      diag.push_back(sqrt(pow(mesh_1->max_[0]-mesh_1->min_[0],2)+pow(mesh_1->max_[1]-mesh_1->min_[1],2)+pow(mesh_1->max_[2]-mesh_1->min_[2],2)));
      diag.push_back(sqrt(pow(mesh_2->max_[0]-mesh_2->min_[0],2)+pow(mesh_2->max_[1]-mesh_2->min_[1],2)+pow(mesh_2->max_[2]-mesh_2->min_[2],2)));
// initialize vertex buffer and write positions to vertex shader
initVertexBuffer();

}
void Viewer5::resizeGL(int w, int h)
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

void Viewer5::paintGL()
{

     if(timer.isCreated()) timer.begin();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // use glsl program
    glUseProgram(m_program[VIEW5].programId());

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

//    cam.Position;
    Eigen::Matrix4f view_inverse= view.inverse();
    // We compute the current cell from the position of the camara
    Eigen::Vector3f cam_pos;
    cam_pos<<view_inverse(12),view_inverse(13),view_inverse(14);
    int i_c = 2*floor(cam_pos[0]);
    if(i_c<0) i_c = abs(i_c)+1;
    int k_c =2*floor(cam_pos[2]);
    if(k_c<0) k_c = abs(k_c)+1;
    int curr_cell=61*i_c+k_c;
    std::vector<int> curr_vis;
    for(int v=0;v<(int)visi.size();v++){// We take the visibility vector from the current cell
    if(visi[v].cellID==curr_cell) curr_vis=visi[v].visibility;
    }

  //  std::cout<<"Cam pos: "<<std::to_string(cam_pos[0])<<" , "<<std::to_string(cam_pos[2])<<" cell: "<<std::to_string(curr_cell)<<std::endl;

    /**DRAW WALLS**/
    glUseProgram(m_program[WALL].programId());
    glUniformMatrix4fv(const_uniform[WALL][CONST_UNIFORM_PROJ], 1, GL_FALSE, projection.data());
    glUniformMatrix4fv(const_uniform[WALL][CONST_UNIFORM_MODEL], 1, GL_FALSE, model.data());
    glUniformMatrix4fv(const_uniform[WALL][CONST_UNIFORM_VIEW], 1, GL_FALSE, view.data());
    glUniformMatrix3fv(const_uniform[WALL][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());


    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");
    Eigen::Matrix4f translation;
    glBindVertexArray(VAO_WALL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_WALL);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].type==3){
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,0,objects[l].z)));
            translation = (kTran.matrix()*model);
            glUniformMatrix4fv(const_uniform[WALL][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());

            glDrawElements(GL_TRIANGLES,size_faces_wall, GL_UNSIGNED_INT, 0);
        }
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);


    /**DRAW FLOOR**/
    glUseProgram(m_program[FLOOR].programId());
    glActiveTexture(GL_TEXTURE1);
    floor_tex->bind(1);
    glUniform1i(glGetUniformLocation(m_program[FLOOR].programId(), "floor_tex"), 1);

    glUniformMatrix4fv(const_uniform[FLOOR][CONST_UNIFORM_PROJ], 1, GL_FALSE, projection.data());
    glUniformMatrix4fv(const_uniform[FLOOR][CONST_UNIFORM_MODEL], 1, GL_FALSE, model.data());
    glUniformMatrix4fv(const_uniform[FLOOR][CONST_UNIFORM_VIEW], 1, GL_FALSE, view.data());
    glUniformMatrix3fv(const_uniform[FLOOR][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());

    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");

    glBindVertexArray(VAO_FLOOR);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_FLOOR);
    glDrawElements(GL_TRIANGLES,size_faces_floor, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);

    if (mesh_ == nullptr || mesh_1==nullptr|| mesh_2==nullptr) return;
    glUseProgram(m_program[VIEW5].programId());


    glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_PROJ], 1, GL_FALSE, projection.data());
    glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE, model.data());
    glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_VIEW], 1, GL_FALSE, view.data());
    glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
    glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "colored"), (int)colored);
    /***Computation of Session 4 ***/
    if(!max_throughput){
        float low_cost=INT_MAX;
        int l_pos=0;
        for(int i=0;i<(int)objects.size();i++){
            std::vector<int>::iterator it;
            it = std::find (curr_vis.begin(), curr_vis.end(), objects[i].cellID);
            if(objects[i].lod!=4&&objects[i].type<3&&it!=curr_vis.end()&&(i!=hysteresis[0]&&i!=hysteresis[1]&&i!=hysteresis[2]&&i!=hysteresis[3]&&i!=hysteresis[4])){
                Eigen::Vector4f view_c = view*(Eigen::Vector4f (objects[i].x,0,objects[i].z,1));//model*
                float D = sqrt(pow(view_c[0],2)+pow(view_c[1],2)+pow(view_c[2],2));
                float cost_diff= pow(2,objects[i].lod+1)*D/diag[objects[i].type]-pow(2,objects[i].lod)*D/diag[objects[i].type];
                if(std::min(low_cost,cost_diff)!=low_cost){
                    low_cost=cost_diff;
                    l_pos=i;
                }
            }

        }
        objects[l_pos].lod+=1;
        hysteresis.push_back(l_pos);
        hysteresis.erase(hysteresis.begin());

    }else{
        float high_cost=INT_MIN;
        int h_pos=0;
        for(int i=0;i<(int)objects.size();i++){
            std::vector<int>::iterator it;
            it = std::find (curr_vis.begin(), curr_vis.end(), objects[i].cellID);
            if(objects[i].lod!=0&&objects[i].type<3&&it!=curr_vis.end()&&(i!=hysteresis[0]&&i!=hysteresis[1]&&i!=hysteresis[2]&&i!=hysteresis[3]&&i!=hysteresis[4])){
                Eigen::Vector4f view_c = view*(Eigen::Vector4f (objects[i].x,0,objects[i].z,1));//model*
                float D = sqrt(pow(view_c[0],2)+pow(view_c[1],2)+pow(view_c[2],2));
                float cost_diff= pow(2,objects[i].lod+1)*D/diag[objects[i].type]-pow(2,objects[i].lod)*D/diag[objects[i].type];
                if(std::max(high_cost,cost_diff)!=high_cost){
                    high_cost=cost_diff;
                    h_pos=i;
                }
            }
        }

        objects[h_pos].lod-=1;
        hysteresis.push_back(h_pos);
        hysteresis.erase(hysteresis.begin());

    }

    // We draw all three models taking into account the visibility
    draw_Armadillo(model,view,curr_vis);
    draw_Bunny(model,view,curr_vis);
    draw_Anno(model,view,curr_vis);

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

    if(fps<32){max_throughput=true;}
    if(fps>=32){max_throughput=false;}


}

QGroupBox* Viewer5::controlPanel()
{
    // init group box
    QGroupBox *groupBox = HW::controlPanel();
    groupBox->setStyleSheet(GroupBoxStyle);
    implementColored(groupBox);
    showFramerate(groupBox);
    return(groupBox);
}

void Viewer5::reset()
{
    mesh_ = nullptr;	// recompute geometry
    initVertexBuffer();
    // draw
    updateGL();
}

void Viewer5::initShaders()
{
    // init uniforms hash table based on uniform variable names and location IDs
    UniformMap uniforms;

    // compile shader, bind attribute vars, link shader, and initialize uniform var table
    initShader(VIEW5, QString(":/SRGGE/Shaders/s4s5s6.vert"), QString(":/SRGGE/Shaders/s4s5s6.frag"), uniforms);
    initShader(WALL, QString(":/SRGGE/Shaders/wall.vert"), QString(":/SRGGE/Shaders/wall.frag"), uniforms);
    initShader(FLOOR, QString(":/SRGGE/Shaders/floor.vert"), QString(":/SRGGE/Shaders/floor.frag"), uniforms);
}

void Viewer5::initVertexBuffer()
{

    //Wall
    initBuffer_WALL();
    //Floor
    initBuffer_FLOOR();

    if (mesh_ != nullptr)
    {
    //Armadillo
    initBuffer_original();
    initBuffer_16();
    initBuffer_32();
    initBuffer_64();
    initBuffer_128();
    }
    if(mesh_1 != nullptr)
    {
    //Bunny
    initBuffer1_original();
    initBuffer1_16();
    initBuffer1_32();
    initBuffer1_64();
    initBuffer1_128();
    }
    if(mesh_2 != nullptr)
    {
        //Anno
        initBuffer2_original();
        initBuffer2_16();
        initBuffer2_32();
        initBuffer2_64();
        initBuffer2_128();
    }
   updateGL();


}


void Viewer5::importModel(){
    QString filename;

      filename = QFileDialog::getOpenFileName(this, tr("Load model"), "./",
                                              tr("PLY Files ( *.ply )"));
      if (!filename.isNull()) {
          if(!mesh_importer.LoadModel(filename))
            QMessageBox::warning(this, tr("Error"),
                               tr("The file could not be opened"));
      }

}



void Viewer5::showFramerate(QGroupBox *groupBox){
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


/**DRAW FUNCTIONS***/
void Viewer5::draw_Armadillo(Eigen::Matrix4f model,Eigen::Matrix4f view,std::vector<int> visibility){
    Eigen::Matrix4f translation;
    /**DRAW ARMADILLO **/
    /*** We draw LOD 0 ***/
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");
    glBindVertexArray(VAO_ID_16);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_16);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==0&&objects[l].type==0){
            // If the object to be drawn is in the current visivility vector, we draw it, otherwise we don't
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
                Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleArm,scaleArm,scaleArm)));
                Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYArm,objects[l].z)));
                Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
                translation = (kTran.matrix()*kSc*piRotation* model);
                Eigen::Matrix4f t = view *piRotation* model;
                Eigen::Matrix3f normal;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

                normal = normal.inverse().transpose();
                glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
                glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
                glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 0);
                glDrawElements(GL_TRIANGLES,mesh_->faces_16.size(), GL_UNSIGNED_INT, 0);
            }//
        }
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);

    //THE SAME WAY FOR THE REST OF THE LODs

    /*** We draw LOD 1 ***/
    glBindVertexArray(VAO_ID_32);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_32);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==1&&objects[l].type==0){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleArm,scaleArm,scaleArm)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYArm,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 1);
            glDrawElements(GL_TRIANGLES,mesh_->faces_32.size(), GL_UNSIGNED_INT, 0);
            }
        }
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 2 ***/
    glBindVertexArray(VAO_ID_64);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_64);
    for(int l=0;l<(int)objects.size();l++){

        if(objects[l].lod==2&&objects[l].type==0){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleArm,scaleArm,scaleArm)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYArm,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 2);
            glDrawElements(GL_TRIANGLES,mesh_->faces_64.size(), GL_UNSIGNED_INT, 0);
        }
        }
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 3 ***/
    glBindVertexArray(VAO_ID_128);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID_128);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==3&&objects[l].type==0){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleArm,scaleArm,scaleArm)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYArm,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 3);
            glDrawElements(GL_TRIANGLES,mesh_->faces_128.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 4 ***/
    glBindVertexArray(VAO_ID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==4&&objects[l].type==0){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleArm,scaleArm,scaleArm)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYArm,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 4);
            glDrawElements(GL_TRIANGLES,mesh_->faces_.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);


}

void Viewer5::draw_Bunny(Eigen::Matrix4f model,Eigen::Matrix4f view,std::vector<int> visibility){
    Eigen::Matrix4f translation;
    /**DRAW Bunny 0 **/
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");
    /*** We draw LOD 0 ***/
    glBindVertexArray(VAO_ID1_16);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_16);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==0&&objects[l].type==1){
            // If the object to be drawn is in the current visivility vector, we draw it, otherwise we don't
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleBun,scaleBun,scaleBun)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYBun,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 0);
            glDrawElements(GL_TRIANGLES,mesh_1->faces_16.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);

    //THE SAME WAY FOR THE REST OF THE LODs

    /*** We draw LOD 1 ***/
    glBindVertexArray(VAO_ID1_32);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_32);
    for(int l=0;l<(int)objects.size();l++){

        if(objects[l].lod==1&&objects[l].type==1){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleBun,scaleBun,scaleBun)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYBun,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 1);
            glDrawElements(GL_TRIANGLES,mesh_1->faces_32.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 2 ***/
    glBindVertexArray(VAO_ID1_64);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_64);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==2&&objects[l].type==1){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleBun,scaleBun,scaleBun)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYBun,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 2);
            glDrawElements(GL_TRIANGLES,mesh_1->faces_64.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 3 ***/
    glBindVertexArray(VAO_ID1_128);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_128);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==3&&objects[l].type==1){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleBun,scaleBun,scaleBun)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYBun,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 3);
            glDrawElements(GL_TRIANGLES,mesh_1->faces_128.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 3 ***/
    glBindVertexArray(VAO_ID1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==4&&objects[l].type==1){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleBun,scaleBun,scaleBun)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYBun,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 4);
            glDrawElements(GL_TRIANGLES,mesh_1->faces_.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
}

void Viewer5::draw_Anno(Eigen::Matrix4f model,Eigen::Matrix4f view, std::vector<int> visibility){
    Eigen::Matrix4f translation;
    /**DRAW ANNO 0 **/
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");
    /*** We draw LOD 0 ***/
    glBindVertexArray(VAO_ID2_16);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_16);
    for(int l=0;l<(int)objects.size();l++){

        if(objects[l].lod==0&&objects[l].type==2){
            // If the object to be drawn is in the current visivility vector, we draw it, otherwise we don't
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleAnn,scaleAnn,scaleAnn)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYAnn,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 0);
            glDrawElements(GL_TRIANGLES,mesh_2->faces_16.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);

        //THE SAME WAY FOR THE REST OF THE LODs

    /*** We draw LOD 1 ***/
    glBindVertexArray(VAO_ID2_32);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_32);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==1&&objects[l].type==2){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleAnn,scaleAnn,scaleAnn)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYAnn,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 1);
            glDrawElements(GL_TRIANGLES,mesh_2->faces_32.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 2 ***/
    glBindVertexArray(VAO_ID2_64);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_64);
    for(int l=0;l<(int)objects.size();l++){

        if(objects[l].lod==2&&objects[l].type==2){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleAnn,scaleAnn,scaleAnn)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYAnn,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 2);
            glDrawElements(GL_TRIANGLES,mesh_2->faces_64.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 3 ***/
    glBindVertexArray(VAO_ID2_128);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_128);
    for(int l=0;l<(int)objects.size();l++){
        if(objects[l].lod==3&&objects[l].type==2){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleAnn,scaleAnn,scaleAnn)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYAnn,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 3);
            glDrawElements(GL_TRIANGLES,mesh_2->faces_128.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);
    /*** We draw LOD 4 ***/
    glBindVertexArray(VAO_ID2);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2);
    for(int l=0;l<(int)objects.size();l++){

        if(objects[l].lod==4&&objects[l].type==2){
            std::vector<int>::iterator it;
            it = std::find (visibility.begin(), visibility.end(), objects[l].cellID);
            if(it!= visibility.end()){
            Eigen::Affine3f kSc (Eigen::Scaling(Eigen::Vector3f(scaleAnn,scaleAnn,scaleAnn)));
            Eigen::Affine3f kTran(Eigen::Translation3f(Eigen::Vector3f(objects[l].x,transYAnn,objects[l].z)));
            Eigen::Affine3f piRotation(Eigen::AngleAxisf(objects[l].rotation*M_PI, Eigen::Vector3f(0.0, 1.0, 0.0)));
            translation = (kTran.matrix()*kSc*piRotation* model);
            Eigen::Matrix4f t = view *piRotation* model;
            Eigen::Matrix3f normal;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

            normal = normal.inverse().transpose();
            glUniformMatrix3fv(const_uniform[VIEW5][CONST_UNIFORM_NORMAL], 1, GL_FALSE, normal.data());
            glUniformMatrix4fv(const_uniform[VIEW5][CONST_UNIFORM_MODEL], 1, GL_FALSE,translation.data());
            glUniform1i(glGetUniformLocation(m_program[VIEW5].programId(), "lod"), 4);
            glDrawElements(GL_TRIANGLES,mesh_2->faces_.size(), GL_UNSIGNED_INT, 0);
        }}
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glBindVertexArray(0);


}

void Viewer5::implementColored(QGroupBox*groupBox){
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
void Viewer5::changeColored(int c){
    colored = c;
    update();
}




//********** GET MODELS POSITIONS**************//
void Viewer5::loadPositions(int type){
    if(type==0){
        std::ifstream infile("Map/Armadillo_pos.txt");
        std::string line;
        int L;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss>>L;

        for(int l=0;l<L;l++){
            std::getline(infile, line);
            std::istringstream iss(line);
            int cell;
            iss>>cell;

            int k =fmod(cell,N);
            int i = (cell-k)/N;

            models m;
            m.cellID=cell;
            m.x=(i+fmod(i,2))/2*cos(M_PI*i);
            m.z=(k+fmod(k,2))/2*cos(M_PI*k);
            m.type=type;
            m.lod=0;

            std::getline(infile, line);
            std::istringstream iss2(line);
            float rot;
            iss2>>rot;
            m.rotation=rot;
            objects.push_back(m);
        }

    }
    if(type==1){
       std::ifstream infile("Map/Bunny_pos.txt");
        std::string line;
        int L;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss>>L;

        for(int l=0;l<L;l++){
            std::getline(infile, line);
            std::istringstream iss(line);
            int cell;
            iss>>cell;

            int k =fmod(cell,N);
            int i = (cell-k)/N;

            models m;
            m.cellID=cell;
            m.x=(i+fmod(i,2))/2*cos(M_PI*i);
            m.z=(k+fmod(k,2))/2*cos(M_PI*k);
            m.type=type;
            m.lod=0;

            std::getline(infile, line);
            std::istringstream iss2(line);
            float rot;
            iss2>>rot;
            m.rotation=rot;
            objects.push_back(m);
        }

    }
    if(type==2){

        std::ifstream infile("Map/Anno_pos.txt");
        std::string line;
        int L;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss>>L;

        for(int l=0;l<L;l++){
            std::getline(infile, line);
            std::istringstream iss(line);
            int cell;
            iss>>cell;

            int k =fmod(cell,N);
            int i = (cell-k)/N;

            models m;
            m.cellID=cell;
            m.x=(i+fmod(i,2))/2*cos(M_PI*i);
            m.z=(k+fmod(k,2))/2*cos(M_PI*k);
            m.type=type;
            m.lod=0;

            std::getline(infile, line);
            std::istringstream iss2(line);
            float rot;
            iss2>>rot;
            m.rotation=rot;
            objects.push_back(m);
        }

    }
    if(type==3){
        std::ifstream infile("Map/Wall_pos.txt");
        std::string line;
        int L;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss>>L;

        for(int l=0;l<L;l++){
            std::getline(infile, line);
            std::istringstream iss(line);
            int cell;
            iss>>cell;

            int k =fmod(cell,N);
            int i = (cell-k)/N;

            models m;
            m.x=(i+fmod(i,2))/2*cos(M_PI*i);
            m.z=(k+fmod(k,2))/2*cos(M_PI*k);
            m.type=type;
            m.lod=0;
            m.rotation=0;
            objects.push_back(m);
        }


    }

}

//********** GET MODELS POSITIONS**************//
void Viewer5::loadVisibility(){

    std::ifstream infile("Map/Visibility.txt");
    std::string line;
    while(std::getline( infile, line )){
      std::istringstream is( line );
      vis v;
      int n;
      is >> n;
      v.cellID=n;
         while( is >> n ) {
        v.visibility.push_back(n);
        }
      visi.push_back(v);
    }
}
void Viewer5::LoadDic(){
    std::ifstream infile("Map/DIC.txt");
    std::string line;
    for(int l=0;l<61*61;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;
        dic.push_back(cellR);
    }
}

//***********INIT VERTEX FOR ALL LEVEL OF DETAIL*******ARMADILLO********//
void Viewer5::initBuffer_original(){

    scaleArm =2/std::max(mesh_->max_[0]-mesh_->min_[0],mesh_->max_[2]-mesh_->min_[2]);
    transYArm=scaleArm*(mesh_->max_[1]+mesh_->min_[1])*0.63;
    diag[0]=diag[0]*scaleArm;

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
void Viewer5::initBuffer_16(){
    std::cout<<"inicializo 16 Armadillo"<<std::endl;
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
void Viewer5::initBuffer_32(){
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
void Viewer5::initBuffer_64(){
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
void Viewer5::initBuffer_128(){
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


//***********INIT VERTEX FOR ALL LEVEL OF DETAIL*******BUNNY********//
void Viewer5::initBuffer1_original(){

    scaleBun=2/std::max(mesh_1->max_[0]-mesh_1->min_[0],mesh_1->max_[2]-mesh_1->min_[2]);
    transYBun=-0.3;//scaleBun*(mesh_1->max_[1]+mesh_1->min_[1])/2;
    diag[1]=diag[1]*scaleBun;
    glGenBuffers(1, &vertexBufferVBO_ID1);
    glGenBuffers(1, &normalBufferVBO_ID1);
    glGenBuffers(1, &facesBufferVBO_ID1);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_1->vertices_.size()/3;i++){

        data.push_back(mesh_1->vertices_[3*i]);
        data.push_back(mesh_1->vertices_[3*i+1]);
        data.push_back(mesh_1->vertices_[3*i+2]);
        data.push_back(mesh_1->normals_[3*i]);
        data.push_back(mesh_1->normals_[3*i+1]);
        data.push_back(mesh_1->normals_[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID1);

    glBindVertexArray(VAO_ID1);
    glGenBuffers(1, &VBO_ID1);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID1);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_1->vertices_.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_1->faces_.size() * sizeof(int), &mesh_1->faces_[0], GL_STATIC_DRAW);
     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
     glBindVertexArray(0);
}
void Viewer5::initBuffer1_16(){
    glGenBuffers(1, &vertexBufferVBO_ID1_16);
    glGenBuffers(1, &normalBufferVBO_ID1_16);
    glGenBuffers(1, &facesBufferVBO_ID1_16);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_1->vertices_16.size()/3;i++){

        data.push_back(mesh_1->vertices_16[3*i]);
        data.push_back(mesh_1->vertices_16[3*i+1]);
        data.push_back(mesh_1->vertices_16[3*i+2]);
        data.push_back(mesh_1->normals_16[3*i]);
        data.push_back(mesh_1->normals_16[3*i+1]);
        data.push_back(mesh_1->normals_16[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID1_16);

    glBindVertexArray(VAO_ID1_16);
    glGenBuffers(1, &VBO_ID1_16);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID1_16);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_1->vertices_16.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_16);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_1->faces_16.size() * sizeof(int), &mesh_1->faces_16[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer5::initBuffer1_32(){
    glGenBuffers(1, &vertexBufferVBO_ID1_32);
    glGenBuffers(1, &normalBufferVBO_ID1_32);
    glGenBuffers(1, &facesBufferVBO_ID1_32);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_1->vertices_32.size()/3;i++){

        data.push_back(mesh_1->vertices_32[3*i]);
        data.push_back(mesh_1->vertices_32[3*i+1]);
        data.push_back(mesh_1->vertices_32[3*i+2]);
        data.push_back(mesh_1->normals_32[3*i]);
        data.push_back(mesh_1->normals_32[3*i+1]);
        data.push_back(mesh_1->normals_32[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID1_32);

    glBindVertexArray(VAO_ID1_32);
    glGenBuffers(1, &VBO_ID1_32);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID1_32);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_1->vertices_32.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_32);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_1->faces_32.size() * sizeof(int), &mesh_1->faces_32[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer5::initBuffer1_64(){
    glGenBuffers(1, &vertexBufferVBO_ID1_64);
    glGenBuffers(1, &normalBufferVBO_ID1_64);
    glGenBuffers(1, &facesBufferVBO_ID1_64);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_1->vertices_64.size()/3;i++){

        data.push_back(mesh_1->vertices_64[3*i]);
        data.push_back(mesh_1->vertices_64[3*i+1]);
        data.push_back(mesh_1->vertices_64[3*i+2]);
        data.push_back(mesh_1->normals_64[3*i]);
        data.push_back(mesh_1->normals_64[3*i+1]);
        data.push_back(mesh_1->normals_64[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID1_64);

    glBindVertexArray(VAO_ID1_64);
    glGenBuffers(1, &VBO_ID1_64);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID1_64);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_1->vertices_64.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_64);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_1->faces_64.size() * sizeof(int), &mesh_1->faces_64[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer5::initBuffer1_128(){
    glGenBuffers(1, &vertexBufferVBO_ID1_128);
    glGenBuffers(1, &normalBufferVBO_ID1_128);
    glGenBuffers(1, &facesBufferVBO_ID1_128);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_1->vertices_128.size()/3;i++){

        data.push_back(mesh_1->vertices_128[3*i]);
        data.push_back(mesh_1->vertices_128[3*i+1]);
        data.push_back(mesh_1->vertices_128[3*i+2]);
        data.push_back(mesh_1->normals_128[3*i]);
        data.push_back(mesh_1->normals_128[3*i+1]);
        data.push_back(mesh_1->normals_128[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID1_128);

    glBindVertexArray(VAO_ID1_128);
    glGenBuffers(1, &VBO_ID1_128);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID1_128);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_1->vertices_128.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID1_128);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_1->faces_128.size() * sizeof(int), &mesh_1->faces_128[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}

//***********INIT VERTEX FOR ALL LEVEL OF DETAIL*******ANNO********//
void Viewer5::initBuffer2_original(){

    scaleAnn =2/std::max(mesh_2->max_[0]-mesh_2->min_[0],mesh_2->max_[2]-mesh_2->min_[2]);
    transYAnn=0;//scaleAnn*(mesh_2->max_[1]+mesh_2->min_[1])*0.63;
    diag[2]=diag[2]*scaleAnn;
    glGenBuffers(1, &vertexBufferVBO_ID2);
    glGenBuffers(1, &normalBufferVBO_ID2);
    glGenBuffers(1, &facesBufferVBO_ID2);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_2->vertices_.size()/3;i++){

        data.push_back(mesh_2->vertices_[3*i]);
        data.push_back(mesh_2->vertices_[3*i+1]);
        data.push_back(mesh_2->vertices_[3*i+2]);
        data.push_back(mesh_2->normals_[3*i]);
        data.push_back(mesh_2->normals_[3*i+1]);
        data.push_back(mesh_2->normals_[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID2);

    glBindVertexArray(VAO_ID2);
    glGenBuffers(1, &VBO_ID2);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID2);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_2->vertices_.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_2->faces_.size() * sizeof(int), &mesh_2->faces_[0], GL_STATIC_DRAW);
     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
     glBindVertexArray(0);
}
void Viewer5::initBuffer2_16(){
    glGenBuffers(1, &vertexBufferVBO_ID2_16);
    glGenBuffers(1, &normalBufferVBO_ID2_16);
    glGenBuffers(1, &facesBufferVBO_ID2_16);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_2->vertices_16.size()/3;i++){

        data.push_back(mesh_2->vertices_16[3*i]);
        data.push_back(mesh_2->vertices_16[3*i+1]);
        data.push_back(mesh_2->vertices_16[3*i+2]);
        data.push_back(mesh_2->normals_16[3*i]);
        data.push_back(mesh_2->normals_16[3*i+1]);
        data.push_back(mesh_2->normals_16[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID2_16);

    glBindVertexArray(VAO_ID2_16);
    glGenBuffers(1, &VBO_ID2_16);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID2_16);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_2->vertices_16.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_16);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_2->faces_16.size() * sizeof(int), &mesh_2->faces_16[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer5::initBuffer2_32(){
    glGenBuffers(1, &vertexBufferVBO_ID2_32);
    glGenBuffers(1, &normalBufferVBO_ID2_32);
    glGenBuffers(1, &facesBufferVBO_ID2_32);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_2->vertices_32.size()/3;i++){

        data.push_back(mesh_2->vertices_32[3*i]);
        data.push_back(mesh_2->vertices_32[3*i+1]);
        data.push_back(mesh_2->vertices_32[3*i+2]);
        data.push_back(mesh_2->normals_32[3*i]);
        data.push_back(mesh_2->normals_32[3*i+1]);
        data.push_back(mesh_2->normals_32[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID2_32);

    glBindVertexArray(VAO_ID2_32);
    glGenBuffers(1, &VBO_ID2_32);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID2_32);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_2->vertices_32.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_32);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_2->faces_32.size() * sizeof(int), &mesh_2->faces_32[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer5::initBuffer2_64(){
    glGenBuffers(1, &vertexBufferVBO_ID2_64);
    glGenBuffers(1, &normalBufferVBO_ID2_64);
    glGenBuffers(1, &facesBufferVBO_ID2_64);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_2->vertices_64.size()/3;i++){

        data.push_back(mesh_2->vertices_64[3*i]);
        data.push_back(mesh_2->vertices_64[3*i+1]);
        data.push_back(mesh_2->vertices_64[3*i+2]);
        data.push_back(mesh_2->normals_64[3*i]);
        data.push_back(mesh_2->normals_64[3*i+1]);
        data.push_back(mesh_2->normals_64[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID2_64);

    glBindVertexArray(VAO_ID2_64);
    glGenBuffers(1, &VBO_ID2_64);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID2_64);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_2->vertices_64.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_64);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_2->faces_64.size() * sizeof(int), &mesh_2->faces_64[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}
void Viewer5::initBuffer2_128(){
    glGenBuffers(1, &vertexBufferVBO_ID2_128);
    glGenBuffers(1, &normalBufferVBO_ID2_128);
    glGenBuffers(1, &facesBufferVBO_ID2_128);

    std::vector<GLfloat> data;
    data.clear();

    for(int i =0;i<(int)mesh_2->vertices_128.size()/3;i++){

        data.push_back(mesh_2->vertices_128[3*i]);
        data.push_back(mesh_2->vertices_128[3*i+1]);
        data.push_back(mesh_2->vertices_128[3*i+2]);
        data.push_back(mesh_2->normals_128[3*i]);
        data.push_back(mesh_2->normals_128[3*i+1]);
        data.push_back(mesh_2->normals_128[3*i+2]);
    }

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_ID2_128);

    glBindVertexArray(VAO_ID2_128);
    glGenBuffers(1, &VBO_ID2_128);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_ID2_128);

     glBufferData(GL_ARRAY_BUFFER,2*mesh_2->vertices_128.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
     glEnableVertexAttribArray(0);
     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float),(void*)(3*sizeof(float)));

     glEnableVertexAttribArray(1);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_ID2_128);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_2->faces_128.size() * sizeof(int), &mesh_2->faces_128[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);
}

//***********INIT VERTEX FOR THE WALL***************//

void Viewer5::initBuffer_WALL(){
    glGenBuffers(1, &vertexBufferVBO_WALL);
    glGenBuffers(1, &normalBufferVBO_WALL);
    glGenBuffers(1, &facesBufferVBO_WALL);



    glActiveTexture(GL_TEXTURE0);
    wall_tex = new QOpenGLTexture(QImage("Map/wall_tex.jpg"));
    wall_tex->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
    wall_tex->setMagnificationFilter(QOpenGLTexture::Linear);
    wall_tex->bind(0);
    glUniform1i(glGetUniformLocation(m_program[WALL].programId(), "wall_tex"), 0);
   // stbi_image_free(dataI);

    float s3=1/sqrt(3);
    std::vector<GLint> wall_faces={0,2,3,0,1,2,1,4,2,1,5,4,5,7,4,5,6,7,6,3,7,6,0,3,2,4,3,7,3,4};
    size_faces_wall=wall_faces.size();
    std::vector<GLfloat> data={1.0,0.0,1.0,s3,-s3,s3,0.0,0.0,//V0,N0,S0,TO
                               1.0,0.0,-1.0,s3,-s3,-s3,1.0,0.0,//V1,N1
                               1.0,4.0,-1.0,s3,s3,-s3,1.0,1.0,//V2,N2
                               1.0,4.0,1.0,s3,-s3,s3,0.0,1.0,//V3,N3
                               -1.0,4.0,-1.0,-s3,s3,-s3,0.0,1.0,//V4,N4
                               -1.0,0.0,-1.0,-s3,-s3,-s3,0.0,0.0,//V5,N5
                               -1.0,0.0,1.0,-s3,-s3,s3,1.0,0.0,//V6,N6
                               -1.0,4.0,1.0,-s3,s3,s3,1.0,1.0//V7,N7
                              };

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_WALL);

    glBindVertexArray(VAO_WALL);
    glGenBuffers(1, &VBO_WALL);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_WALL);

    glBufferData(GL_ARRAY_BUFFER,data.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float),(void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8*sizeof(float),(void*)(6*sizeof(float)));
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_WALL);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_faces_wall * sizeof(int), &wall_faces[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);

}
//***********INIT VERTEX FOR THE FLOOR***************//
void Viewer5::initBuffer_FLOOR(){
    glGenBuffers(1, &vertexBufferVBO_FLOOR);
    glGenBuffers(1, &normalBufferVBO_FLOOR);
    glGenBuffers(1, &facesBufferVBO_FLOOR);


    glActiveTexture(GL_TEXTURE1);
    floor_tex = new QOpenGLTexture(QImage("Map/Map.png"));
    floor_tex->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
    floor_tex->setMagnificationFilter(QOpenGLTexture::Linear);
    floor_tex->bind(1);
    glUniform1i(glGetUniformLocation(m_program[FLOOR].programId(), "floor_tex"), 1);

    std::vector<GLint> floor_faces={0,1,2,0,2,3};
    size_faces_floor=floor_faces.size();
    std::vector<GLfloat> data={61,0.0,61,0,1,0,1.0,0.0,//V0,N0
                               61,0.0,-61,0,1,0,1.0,1.0,//V1,N1
                               -61,0.0,-61,0,1,0,0.0,1.0,//V2,N2
                               -61,0.0,61,0,1,0,0.0,0.0//V3,N3
                              };

    glGenVertexArrays = (_glGenVertexArrays) QGLWidget::context()->getProcAddress("glGenVertexArrays");
    glBindVertexArray = (_glBindVertexArray) QGLWidget::context()->getProcAddress("glBindVertexArray");


    glGenVertexArrays(1, &VAO_FLOOR);

    glBindVertexArray(VAO_FLOOR);
    glGenBuffers(1, &VBO_FLOOR);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_FLOOR);

    glBufferData(GL_ARRAY_BUFFER,data.size() * sizeof(GLfloat), &data[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float),(void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8*sizeof(float),(void*)(6*sizeof(float)));
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facesBufferVBO_FLOOR);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_faces_floor * sizeof(int), &floor_faces[0], GL_STATIC_DRAW);
glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
glBindVertexArray(0);

camera_.UpdateModel(Eigen::Vector3f(-1,0,-1),Eigen::Vector3f(1,0,1));
}
