
#include "mesh_importer.h"

MI::MI(HW *hw){
    m_hw = hw;
}



void MI::AddImportExportPLY(QGroupBox*groupBox )
{
    QPushButton *buttonImport  = new QPushButton("Import model");
    QPushButton *buttonExport  = new QPushButton("Export model");
    QLabel *m_vertexText = new QLabel("Vertex:  ");
    m_vertexCount = new QLabel("0");

    QLabel *m_facesText = new QLabel("Faces:  ");
    m_facesCount = new QLabel("0");

    QFrame* line = new QFrame();
    line->setFrameShape(QFrame::HLine);
    QLabel *title = new QLabel(" PLY 3 Extension ");

    HW::connect(buttonImport, SIGNAL(clicked()), m_hw,SLOT(importModel()));
    HW::connect(buttonExport, SIGNAL(clicked()), m_hw,SLOT(exportModel()));
    auto layout = dynamic_cast<QGridLayout*>(groupBox->layout());

    int row = layout->rowCount() + 1;

    layout->addWidget(line,row,0, 1 ,layout->columnCount());
    row++;
    layout->addWidget(title,row,layout->columnCount()/2);
    row++;
    layout->addWidget(buttonImport,row,0);
    layout->addWidget(buttonExport,row,1);
    row++;
    layout->addWidget(m_vertexText,row,0);
    layout->addWidget(m_vertexCount,row,1);
    row++;
    layout->addWidget(m_facesText,row,0);
    layout->addWidget(m_facesCount,row,1);

    groupBox->setLayout(layout);
}

bool MI::Load3Models(std::vector<QString> filename,int m){

    std::unique_ptr<data_representation::TriangleMesh> mesh =
    std::unique_ptr<data_representation::TriangleMesh>(new data_representation::TriangleMesh);

    std::string file = filename[4].toUtf8().constData();
    uint pos = file.find_last_of(".");
    std::string type = file.substr(pos + 1);

    bool res = false;
    if (type.compare("ply") == 0) {
      res = data_representation::ReadFromPlyMuseum4(file, mesh.get());
    }
    file = filename[3].toUtf8().constData();
    pos = file.find_last_of(".");
    type = file.substr(pos + 1);
    if (type.compare("ply") == 0) {
      res = data_representation::ReadFromPlyMuseum3(file, mesh.get());
    }
    file = filename[2].toUtf8().constData();
    pos = file.find_last_of(".");
    type = file.substr(pos + 1);
    if (type.compare("ply") == 0) {
      res = data_representation::ReadFromPlyMuseum2(file, mesh.get());
    }
    file = filename[1].toUtf8().constData();
    pos = file.find_last_of(".");
    type = file.substr(pos + 1);
    if (type.compare("ply") == 0) {
      res = data_representation::ReadFromPlyMuseum1(file, mesh.get());
    }
    file = filename[0].toUtf8().constData();
    pos = file.find_last_of(".");
    type = file.substr(pos + 1);
    if (type.compare("ply") == 0) {
      res = data_representation::ReadFromPlyMuseum0(file, mesh.get());
    }


    if (res&&m==0) {
        m_hw->mesh_.reset(mesh.release());
        return true;
    }
    if (res&&m==1) {
        m_hw->mesh_1.reset(mesh.release());
        return true;
    }
    if (res&&m==2) {
        m_hw->mesh_2.reset(mesh.release());
        return true;
    }


    return false;

}

bool MI::LoadModel(QString filename) {
  std::string file = filename.toUtf8().constData();
  uint pos = file.find_last_of(".");
  std::string type = file.substr(pos + 1);

  std::unique_ptr<data_representation::TriangleMesh> mesh =
  std::unique_ptr<data_representation::TriangleMesh>(new data_representation::TriangleMesh);

  bool res = false;
  if (type.compare("ply") == 0) {
    res = data_representation::ReadFromPly(file, mesh.get());
  }

  if (res) {
      m_hw->mesh_.reset(mesh.release());

      m_hw->camera_.UpdateModel(m_hw->mesh_->min_, m_hw->mesh_->max_);

      m_facesCount->setText(QString(std::to_string(m_hw->mesh_->faces_.size() / 3).c_str()));
      m_vertexCount->setText(QString(std::to_string(m_hw->mesh_->vertices_.size() / 3).c_str()));

      m_hw->camera_.UpdateModel(m_hw->mesh_->min_, m_hw->mesh_->max_);

      m_hw->initVertexBuffer();
      return true;
  }

  return false;
}

bool MI::SaveModel(QString filename,int lod,data_representation::TriangleMesh &mesh) {

  std::string file = filename.toUtf8().constData();
  uint pos = file.find_last_of(".");
  std::string type = file.substr(pos + 1);

  bool res = false;
  if (type.compare("ply") == 0) {
    res = data_representation::WriteToPly(file,mesh,lod);
  }

  return res;
}
void MI::updateNumVert(int num_vert,int num_face){
    m_facesCount->setText(QString(std::to_string(num_face).c_str()));
    m_vertexCount->setText(QString(std::to_string(num_vert).c_str()));
}
