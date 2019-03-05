#ifndef MESH_IMPORTER_H
#define MESH_IMPORTER_H

#include "HW.h"

#include <QtWidgets>
#include <QGLWidget>
#include <QGLFunctions>
#include <QGLShaderProgram>
#include <QtOpenGL>
#include <Helpers/mesh_io.h>
#include <Helpers/triangle_mesh.h>

class MI
{
public:
    MI              (HW *hw);
    bool            LoadModel(QString filename);
    bool            SaveModel(QString filename,int lod,data_representation::TriangleMesh &mesh) ;
    bool            Load3Models(std::vector<QString> filename,int m);
    void            AddImportExportPLY(QGroupBox*groupBox );
    void            updateNumVert(int num_vert,int num_face);

private:
    QLabel *m_vertexCount;
    QLabel *m_facesCount;
    HW     *m_hw;

};

#endif // MESH_IMPORTER_H
