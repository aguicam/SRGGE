// Author: Marc Comino 2018

#include "./mesh_io.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <math.h>


#include "./triangle_mesh.h"

namespace data_representation {

namespace {
int countF,countInv,countOutcell=0;


template <typename T>
void Add3Items(T i1, T i2, T i3, int index, std::vector<T> *vector) {
  (*vector)[index] = i1;
  (*vector)[index + 1] = i2;
  (*vector)[index + 2] = i3;
}

bool ReadPlyHeader(std::ifstream *fin, int *vertices, int *faces) {
  char line[100];

  fin->getline(line, 100);
  if (strncmp(line, "ply", 3) != 0) return false;

  *vertices = 0;
  fin->getline(line, 100);
  while (strncmp(line, "end_header", 10) != 0) {
    if (strncmp(line, "element vertex", 14) == 0) *vertices = atoi(&line[15]);
    fin->getline(line, 100);
    if (strncmp(line, "element face", 12) == 0) *faces = atoi(&line[13]);
  }

  if (*vertices <= 0) return false;

  std::cout << "Loading triangle mesh" << std::endl;
  std::cout << "\tVertices = " << *vertices << std::endl;
  std::cout << "\tFaces = " << *faces << std::endl;

  return true;
}

void ReadPlyVertices(std::ifstream *fin, TriangleMesh *mesh) {
  const int kVertices = mesh->vertices_.size() / 3;
  for (int i = 0; i < kVertices; ++i) {
    float x, y, z;
    fin->read(reinterpret_cast<char *>(&x), sizeof(float));
    fin->read(reinterpret_cast<char *>(&y), sizeof(float));
    fin->read(reinterpret_cast<char *>(&z), sizeof(float));

    Add3Items(x, y, z, i * 3, &(mesh->vertices_));
  }
}

void ReadPlyFaces(std::ifstream *fin, TriangleMesh *mesh) {
  unsigned char vertex_per_face;

  const int kFaces = mesh->faces_.size() / 3;
  for (int i = 0; i < kFaces; ++i) {
    int v1, v2, v3;
    fin->read(reinterpret_cast<char *>(&vertex_per_face),
              sizeof(unsigned char));
    assert(vertex_per_face == 3);

    fin->read(reinterpret_cast<char *>(&v1), sizeof(int));
    fin->read(reinterpret_cast<char *>(&v2), sizeof(int));
    fin->read(reinterpret_cast<char *>(&v3), sizeof(int));
    Add3Items(v1, v2, v3, i * 3, &(mesh->faces_));
  }
}
void ReadPlyVerticesM(std::ifstream *fin, std::vector<float> *vertices_) {
  const int kVertices = vertices_->size() / 3;
  for (int i = 0; i < kVertices; ++i) {
    float x, y, z;
    fin->read(reinterpret_cast<char *>(&x), sizeof(float));
    fin->read(reinterpret_cast<char *>(&y), sizeof(float));
    fin->read(reinterpret_cast<char *>(&z), sizeof(float));

    Add3Items(x, y, z, i * 3, (vertices_));
  }
}

void ReadPlyFacesM(std::ifstream *fin, std::vector<int> *faces_) {
  unsigned char vertex_per_face;

  const int kFaces = faces_->size() / 3;
  for (int i = 0; i < kFaces; ++i) {
    int v1, v2, v3;
    fin->read(reinterpret_cast<char *>(&vertex_per_face),
              sizeof(unsigned char));
    assert(vertex_per_face == 3);

    fin->read(reinterpret_cast<char *>(&v1), sizeof(int));
    fin->read(reinterpret_cast<char *>(&v2), sizeof(int));
    fin->read(reinterpret_cast<char *>(&v3), sizeof(int));
    Add3Items(v1, v2, v3, i * 3, (faces_));
  }
}

void ComputeVertexNormals(const std::vector<float> &vertices,
                          const std::vector<int> &faces,
                          std::vector<float> *normals) {
  const int kFaces = faces.size();
  std::vector<float> face_normals(kFaces, 0);

  for (int i = 0; i < kFaces; i += 3) {

    Eigen::Vector3d v1(vertices[faces[i] * 3], vertices[faces[i] * 3 + 1],
                       vertices[faces[i] * 3 + 2]);
    Eigen::Vector3d v2(vertices[faces[i + 1] * 3],
                       vertices[faces[i + 1] * 3 + 1],
                       vertices[faces[i + 1] * 3 + 2]);
    Eigen::Vector3d v3(vertices[faces[i + 2] * 3],
                       vertices[faces[i + 2] * 3 + 1],
                       vertices[faces[i + 2] * 3 + 2]);
    Eigen::Vector3d v1v2 = v2 - v1;
    Eigen::Vector3d v1v3 = v3 - v1;
    Eigen::Vector3d normal = v1v2.cross(v1v3);

    if (normal.norm() > 0) {
      normal.normalize();
    } else {
      normal = Eigen::Vector3d(0.0, 0.0, 0.0);
    }

    for (int j = 0; j < 3; ++j) face_normals[i + j] = normal[j];
  }

  const int kVertices = vertices.size();
  normals->resize(kVertices, 0);
  for (int i = 0; i < kFaces; i += 3) {
    for (int j = 0; j < 3; ++j) {
      int idx = faces[i + j];
      Eigen::Vector3d v1(vertices[faces[i + j] * 3],
                         vertices[faces[i + j] * 3 + 1],
                         vertices[faces[i + j] * 3 + 2]);
      Eigen::Vector3d v2(vertices[faces[i + (j + 1) % 3] * 3],
                         vertices[faces[i + (j + 1) % 3] * 3 + 1],
                         vertices[faces[i + (j + 1) % 3] * 3 + 2]);
      Eigen::Vector3d v3(vertices[faces[i + (j + 2) % 3] * 3],
                         vertices[faces[i + (j + 2) % 3] * 3 + 1],
                         vertices[faces[i + (j + 2) % 3] * 3 + 2]);

      Eigen::Vector3d v1v2 = v2 - v1;
      Eigen::Vector3d v1v3 = v3 - v1;
      double angle = acos(v1v2.dot(v1v3) / (v1v2.norm() * v1v3.norm()));

      if (angle == angle) {
        for (int k = 0; k < 3; ++k) {
          (*normals)[idx * 3 + k] += face_normals[i + k] * angle;
        }
      }
    }
  }

  const int kNormals = normals->size();
  for (int i = 0; i < kNormals; i += 3) {
    Eigen::Vector3d normal((*normals)[i], (*normals)[i + 1], (*normals)[i + 2]);
    if (normal.norm() > 0) {
      normal.normalize();
    } else {
      normal = Eigen::Vector3d(0, 0, 0);
    }

    for (int j = 0; j < 3; ++j) (*normals)[i + j] = normal[j];
  }
}

void ComputeBoundingBox(const std::vector<float> vertices, TriangleMesh *mesh) {
  const int kVertices = vertices.size() / 3;
  for (int i = 0; i < kVertices; ++i) {
    mesh->min_[0] = std::min(mesh->min_[0], vertices[i * 3]);
    mesh->min_[1] = std::min(mesh->min_[1], vertices[i * 3 + 1]);
    mesh->min_[2] = std::min(mesh->min_[2], vertices[i * 3 + 2]);

    mesh->max_[0] = std::max(mesh->max_[0], vertices[i * 3]);
    mesh->max_[1] = std::max(mesh->max_[1], vertices[i * 3 + 1]);
    mesh->max_[2] = std::max(mesh->max_[2], vertices[i * 3 + 2]);
  }
}



//**************Quadratic Error***********///////
std::vector<float> quadraticErrorMetrics(std::vector<int> vertex_i,TriangleMesh *mesh, int cellID, int N){

    std::vector<float> res_vert;
    Eigen::Vector4f res;
    Eigen::Matrix4f Q;
    Q << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;

    // For each vertex in the cluster, we compute the Q matrix given by the vertex p*pt, where the 3 first components
    // are the x y z of the normal, and the forth component is the negative of the dot product of the vertex position and its normal
    // We add all the Q matrices together.
    for(int i=0;i<(int)vertex_i.size();i++){
        Eigen::Vector3f n1,v1;
        Eigen::Vector4f p1;
        Eigen::Matrix4f Qx;
        n1<<mesh->normals_[3*vertex_i[i]],mesh->normals_[3*vertex_i[i]+1],mesh->normals_[3*vertex_i[i]+2];
        v1<<mesh->vertices_[3*vertex_i[i]],mesh->vertices_[3*vertex_i[i]+1],mesh->vertices_[3*vertex_i[i]+2];
        p1<< n1[0],n1[1],n1[2],-n1.dot(v1);
        Qx=p1*p1.transpose();
        Q=Q+Qx;
    }

    //We compute the mean of all the vertices
    float vertex_mean[3]={0,0,0};
    for(int v=0;v<(int)vertex_i.size();v++){
        vertex_mean[0]=vertex_mean[0]+mesh->vertices_[3*vertex_i[v]];
        vertex_mean[1]=vertex_mean[1]+mesh->vertices_[3*vertex_i[v]+1];
        vertex_mean[2]=vertex_mean[2]+mesh->vertices_[3*vertex_i[v]+2];
    }

    vertex_mean[0]=vertex_mean[0]/vertex_i.size();
    vertex_mean[1]=vertex_mean[1]/vertex_i.size();
    vertex_mean[2]=vertex_mean[2]/vertex_i.size();

    //In order to make Q more invertible
    // We put the last row to 0 0 0 1
    Q(3,0)=0;
    Q(3,1)=0;
    Q(3,2)=0;
    Q(3,3)=1;
    // add 1 to the diagonal of the 3x3 matrix
    Q(0,0)+=1;
    Q(1,1)+=1;
    Q(2,2)+=1;
    // we substract the mean to the last column
    Q(0,3)-=vertex_mean[0];
    Q(1,3)-=vertex_mean[1];
    Q(2,3)-=vertex_mean[2];


    Eigen::Vector4f v_ (0.0,0.0,0.0,1.0);
    Eigen::Matrix4f Q_inv;
    bool inv=false;
    bool out_cell=false;
    // We try to invert Q
    Q.computeInverseWithCheck(Q_inv,inv,0.1);

    //We compute the distance between center of cell and one of the corners;
    float seg_len = sqrt(pow((mesh->max_[0]-mesh->min_[0])/(2*N),2)+pow((mesh->max_[1]-mesh->min_[1])/(2*N),2)+pow((mesh->max_[2]-mesh->min_[2])/(2*N),2));

    if(inv){ // It Q invertible
        // We compute the vertex
        res=Q_inv*v_;
        cellID=floor(cellID/6);// due to the  shape preservation algorithm
        // We check if the new vertex is within a certain distance of the cell
        int k= fmod(cellID,N);
        int j = (fmod(cellID,pow(N,2))-k)/N;
        int i =(cellID -j*N-k)/pow(N,2);

         //****CERTAIN DISTANCE FROM THE CENTER OF THE CELL DISCRIMINATION
        float minL=mesh->min_[0]+i*(mesh->max_[0]-mesh->min_[0])/N;
        float maxL=mesh->min_[0]+(i+1)*(mesh->max_[0]-mesh->min_[0])/N;
        float centerX=(minL+maxL)/2;
         minL=mesh->min_[1]+j*(mesh->max_[1]-mesh->min_[1])/N;
         maxL=mesh->min_[1]+(j+1)*(mesh->max_[1]-mesh->min_[1])/N;
        float centerY=(minL+maxL)/2;
         minL=mesh->min_[2]+k*(mesh->max_[2]-mesh->min_[2])/N;
         maxL=mesh->min_[2]+(k+1)*(mesh->max_[2]-mesh->min_[2])/N;
        float centerZ=(minL+maxL)/2;


        //distance from vertex to center of cell
        float d_r_c=sqrt(pow(centerX-res[0],2)+pow(centerY-res[1],2)+pow(centerZ-res[2],2));

        //We compare distances and decide if we keep the vertex or not, the factor is 2;
        if (d_r_c>16*seg_len/N)out_cell=true;

    }
    if(!inv) countInv++;
    if(out_cell)countOutcell++;
    if(!inv || out_cell){ // If Q not invertible or the vertex is too far from the cell we take the mean
        countF++;
        res<<vertex_mean[0],vertex_mean[1],vertex_mean[2],0;
    }

    res_vert.push_back(res(0));
    res_vert.push_back(res(1));
    res_vert.push_back(res(2));

    return res_vert;
}

int shapePreservation(int i,int j,int k, int v,TriangleMesh *mesh,int N){

    // We look at what direction is the vertex normal pointing at, and assign it to a cell subset of the original cell
    //x+ 0  // x- 1 // y+ 2 // y- 3 // z+ 4 // z- 5
    int index =6*(i*pow(N,2)+j*N+k);
    Eigen::Vector3f normal (mesh->normals_[3*v],mesh->normals_[3*v+1],mesh->normals_[3*v+2]);
    int pos=0;
    float min_angle= FLT_MAX;
    // X+
    if(std::min(min_angle,(float)acos(normal.dot(Eigen::Vector3f(1.0,0.0,0.0))/normal.norm()))!=min_angle){
        pos=0;
    }
    // X-
    if(std::min(min_angle,(float)acos(normal.dot(Eigen::Vector3f(-1.0,0.0,0.0))/normal.norm()))!=min_angle){
        pos=1;
    }
    // Y+
    if(std::min(min_angle,(float)acos(normal.dot(Eigen::Vector3f(0.0,1.0,0.0))/normal.norm()))!=min_angle){
        pos=2;
    }
    // Y-
    if(std::min(min_angle,(float)acos(normal.dot(Eigen::Vector3f(0.0,-1.0,0.0))/normal.norm()))!=min_angle){
        pos=3;
    }
    // Z+
    if(std::min(min_angle,(float)acos(normal.dot(Eigen::Vector3f(0.0,0.0,1.0))/normal.norm()))!=min_angle){
        pos=4;
    }
    // Z-
    if(std::min(min_angle,(float)acos(normal.dot(Eigen::Vector3f(0.0,0.0,-1.0))/normal.norm()))!=min_angle){
        pos=5;
    }

    return index+pos;
}


void vertexCluster(int N,TriangleMesh *mesh){

    std::vector<ce_vert> cells2;
    std::vector<int> vertex_cell;
    std::vector<int> faces_simp;
    std::vector<int> faces_temp;
    std::vector<float> vertex_simp;
    cells2.clear();
    vertex_cell.clear();
    bool in=false;
    int one_per_cell=0;
    for(int v=0;v< (int) mesh->vertices_.size()/3;v++){ //Size would be different for shape preservation
        in=false;
        //We calculate in what segment of each axis of the grid the vertex is found
        int i= fmod((N*(mesh->vertices_[3*v]-mesh->min_[0])/(mesh->max_[0]-mesh->min_[0])),N);
        int j =fmod((N*(mesh->vertices_[3*v+1]-mesh->min_[1])/(mesh->max_[1]-mesh->min_[1])),N);
        int k =fmod((N*(mesh->vertices_[3*v+2]-mesh->min_[2])/(mesh->max_[2]-mesh->min_[2])),N);
        // We calculate the index of the cell where the vertex is found
       // int index_cell = i*pow(N,2)+j*N+k; // when we don't use shape preservation
        //
        int index_cell = shapePreservation(i,j,k,v,mesh,N);// We use Shape Preservation
        int ind=0;
        //We look if the cell exists, if it does, we add the vertex index, if it doesnt, we add the new cell to the vector
        for(int s=0;s<(int)cells2.size();s++){
            if(cells2[s].cell==index_cell){cells2[s].inices_vertex.push_back(v);
            in=true;
            ind=s;}
        }
        if(!in){
            ce_vert cv;
            cv.inices_vertex.push_back(v);
            cv.cell=index_cell;
            cells2.push_back(cv);
            ind=cells2.size()-1;
        }
        vertex_cell.push_back(ind);
    }
    faces_temp.clear();
    faces_temp = mesh->faces_;
    //We modify the vertex index from the faces vector with the new vertex index (cell index)
    for(int f = 0; f<(int)faces_temp.size();f++){
        faces_temp.at(f)=vertex_cell[mesh->faces_[f]];
     }
    // We create the final faces vector by keeping all the faces that have vertices in 3 different cells
    faces_simp.clear();
    for (int f=0;f<(int)faces_temp.size()/3;f++){
        if((faces_temp[3*f]!=faces_temp[3*f+1]) && (faces_temp[3*f+2]!=faces_temp[3*f+1])&&(faces_temp[3*f+2]!=faces_temp[3*f])){

                faces_simp.push_back(faces_temp[3*f]);
                faces_simp.push_back(faces_temp[3*f+1]);
                faces_simp.push_back(faces_temp[3*f+2]);
        }

    }
    //We create the new vertices
    vertex_simp.clear();
    for(int c=0;c<(int)cells2.size();c++){

        //****************Vertex Clustering*****SESSION 2*********///

        //We mean all the vertices in each cell
 /*      float vertex[3]={0,0,0};
        for(int v=0;v<(int)cells2[c].inices_vertex.size();v++){
            vertex[0]=vertex[0]+mesh->vertices_[3*cells2[c].inices_vertex[v]];
            vertex[1]=vertex[1]+mesh->vertices_[3*cells2[c].inices_vertex[v]+1];
            vertex[2]=vertex[2]+mesh->vertices_[3*cells2[c].inices_vertex[v]+2];
        }
        float h = vertex[0]/cells2[c].inices_vertex.size();
        vertex_simp.push_back(h);
        h =vertex[1]/cells2[c].inices_vertex.size();
        vertex_simp.push_back(h);
        h=vertex[2]/cells2[c].inices_vertex.size();
        vertex_simp.push_back(h);  */

        //****************Vertex Clustering: QuadraticErrorMetrics*****SESSION 3*********///

        //We pick the vertex representative using quadric error metrics for cells with more than one vertix,
        //otherwise the vertex representative is the only vertex in the cell
        std::vector<float> vert;
        vert.clear();
        if(cells2[c].inices_vertex.size()>1){
            vert= quadraticErrorMetrics(cells2[c].inices_vertex,mesh,cells2[c].cell,N);
        }
        else  if(cells2[c].inices_vertex.size()==1){
            one_per_cell++;
            vert.push_back(mesh->vertices_[3*cells2[c].inices_vertex[0]]);
            vert.push_back(mesh->vertices_[3*cells2[c].inices_vertex[0]+1]);
            vert.push_back(mesh->vertices_[3*cells2[c].inices_vertex[0]+2]);
        }
        vertex_simp.push_back(vert[0]);
        vertex_simp.push_back(vert[1]);
        vertex_simp.push_back(vert[2]); /////

    }
    std::cout<<"N: "<<std::to_string(N)<<std::endl;
    std::cout<<"Total cells with vertices: "<<std::to_string(vertex_simp.size()/3)<<std::endl;
    std::cout<<"Cells with one vertex: "<<std::to_string(one_per_cell)<<std::endl;
    std::cout<<"Cells that go to use the quadtratic error metrics: "<<std::to_string(vertex_simp.size()/3-one_per_cell)<<std::endl;
    std::cout<<"In q_error vertex go to mean: "<<std::to_string(countF)<<std::endl;
    std::cout<<"In q_error vertex go to mean due to vertex calcluated out of cell: "<<std::to_string(countOutcell)<<std::endl;
    std::cout<<"In q_error vertex go to mean due to invertible Q: "<<std::to_string(countInv)<<std::endl;
    countF=0;
    countInv=0;
    countOutcell=0;

    std::cout<<"Total vertices: "<<std::to_string(vertex_simp.size()/3)<<std::endl;
    std::cout<<"Total faces: "<<std::to_string(faces_simp.size()/3)<<std::endl;

    // We store for each LoD the faces and the vertices.
    switch(N){
    case 16:
        mesh->faces_16=faces_simp;
        mesh->vertices_16=vertex_simp;
        break;
    case 32:
        mesh->faces_32=faces_simp;
        mesh->vertices_32=vertex_simp;
        break;
    case 64:
        mesh->faces_64=faces_simp;
        mesh->vertices_64=vertex_simp;
        break;
    case 128:
        mesh->faces_128=faces_simp;
        mesh->vertices_128=vertex_simp;
        break;
    }

}

}  // namespace

bool ReadFromPly(const std::string &filename, TriangleMesh *mesh) {
  std::ifstream fin;

  fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!fin.is_open() || !fin.good()) return false;

  int vertices = 0, faces = 0;
  if (!ReadPlyHeader(&fin, &vertices, &faces)) {
    fin.close();
    return false;
  }

  mesh->vertices_.resize(vertices * 3);
  ReadPlyVertices(&fin, mesh);

  mesh->faces_.resize(faces * 3);
  ReadPlyFaces(&fin, mesh);

  fin.close();

  ComputeVertexNormals(mesh->vertices_, mesh->faces_, &mesh->normals_);
  ComputeBoundingBox(mesh->vertices_, mesh);

  vertexCluster(128,mesh);
  ComputeVertexNormals(mesh->vertices_128, mesh->faces_128, &mesh->normals_128);
  vertexCluster(64,mesh);
  ComputeVertexNormals(mesh->vertices_64, mesh->faces_64, &mesh->normals_64);
  vertexCluster(32,mesh);
  ComputeVertexNormals(mesh->vertices_32, mesh->faces_32, &mesh->normals_32);
  vertexCluster(16,mesh);
  ComputeVertexNormals(mesh->vertices_16, mesh->faces_16, &mesh->normals_16);
  return true;
}

bool ReadFromPlyMuseum4(const std::string &filename, TriangleMesh *mesh) {
  std::ifstream fin;

  fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!fin.is_open() || !fin.good()) return false;

  int vertices = 0, faces = 0;
  if (!ReadPlyHeader(&fin, &vertices, &faces)) {
    fin.close();
    return false;
  }

  mesh->vertices_.resize(vertices * 3);
  ReadPlyVertices(&fin, mesh);

  mesh->faces_.resize(faces * 3);
  ReadPlyFaces(&fin, mesh);

  fin.close();

  ComputeVertexNormals(mesh->vertices_, mesh->faces_, &mesh->normals_);
  ComputeBoundingBox(mesh->vertices_, mesh);
  return true;
}

bool ReadFromPlyMuseum0(const std::string &filename, TriangleMesh *mesh) {
  std::ifstream fin;

  fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!fin.is_open() || !fin.good()) return false;

  int vertices = 0, faces = 0;
  if (!ReadPlyHeader(&fin, &vertices, &faces)) {
    fin.close();
    return false;
  }

  mesh->vertices_16.resize(vertices * 3);
  ReadPlyVerticesM(&fin, &mesh->vertices_16);

  mesh->faces_16.resize(faces * 3);
  ReadPlyFacesM(&fin, &mesh->faces_16);

  fin.close();

  ComputeVertexNormals(mesh->vertices_16, mesh->faces_16, &mesh->normals_16);
  return true;
}
bool ReadFromPlyMuseum1(const std::string &filename, TriangleMesh *mesh) {
  std::ifstream fin;

  fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!fin.is_open() || !fin.good()) return false;

  int vertices = 0, faces = 0;
  if (!ReadPlyHeader(&fin, &vertices, &faces)) {
    fin.close();
    return false;
  }

  mesh->vertices_32.resize(vertices * 3);
  ReadPlyVerticesM(&fin,&mesh->vertices_32);

  mesh->faces_32.resize(faces * 3);
  ReadPlyFacesM(&fin, &mesh->faces_32);

  fin.close();

  ComputeVertexNormals(mesh->vertices_32, mesh->faces_32, &mesh->normals_32);
  return true;
}
bool ReadFromPlyMuseum2(const std::string &filename, TriangleMesh *mesh) {
  std::ifstream fin;

  fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!fin.is_open() || !fin.good()) return false;

  int vertices = 0, faces = 0;
  if (!ReadPlyHeader(&fin, &vertices, &faces)) {
    fin.close();
    return false;
  }

  mesh->vertices_64.resize(vertices * 3);
  ReadPlyVerticesM(&fin, &mesh->vertices_64);

  mesh->faces_64.resize(faces * 3);
  ReadPlyFacesM(&fin, &mesh->faces_64);

  fin.close();

  ComputeVertexNormals(mesh->vertices_64, mesh->faces_64, &mesh->normals_64);
  return true;
}
bool ReadFromPlyMuseum3(const std::string &filename, TriangleMesh *mesh) {
  std::ifstream fin;

  fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!fin.is_open() || !fin.good()) return false;

  int vertices = 0, faces = 0;
  if (!ReadPlyHeader(&fin, &vertices, &faces)) {
    fin.close();
    return false;
  }

  mesh->vertices_128.resize(vertices * 3);
  ReadPlyVerticesM(&fin, &mesh->vertices_128);

  mesh->faces_128.resize(faces * 3);
  ReadPlyFacesM(&fin, &mesh->faces_128);

  fin.close();

  ComputeVertexNormals(mesh->vertices_128, mesh->faces_128, &mesh->normals_128);
  return true;
}



bool WriteToPly(const std::string &filename, TriangleMesh &mesh,const int lod) {
//  (void)filename;
 // (void)mesh;

    std::ofstream outFile( filename.c_str(), std::ios::binary );


     if ( !outFile )
     {
         std::cerr << "Error opening output file: " << filename << "!" << std::endl;
         exit( 1 );
     }

     ////
     // Header
     ////

     int pointNum =0;
     int triangleNum=0;
     switch(lod){
     case 0:
         pointNum    = ( int ) mesh.vertices_16.size()/3;
         triangleNum = ( int ) mesh.faces_16.size()/3;
         break;
     case 1:
          pointNum    = ( int ) mesh.vertices_32.size()/3;
         triangleNum = ( int ) mesh.faces_32.size()/3;
         break;
     case 2:
         pointNum    = ( int ) mesh.vertices_64.size()/3;
         triangleNum = ( int ) mesh.faces_64.size()/3;
         break;
     case 3:
         pointNum    = ( int ) mesh.vertices_128.size()/3;
         triangleNum = ( int ) mesh.faces_128.size()/3;
         break;
     case 4:
         pointNum    = ( int ) mesh.vertices_.size()/3;
         triangleNum = ( int ) mesh.faces_.size()/3;
         break;

     }


     outFile << "ply" << std::endl;
     outFile << "format ascii 1.0" << std::endl;
     outFile << "element vertex " << pointNum << std::endl;
     outFile << "property float x" << std::endl;
     outFile << "property float y" << std::endl;
     outFile << "property float z" << std::endl;
     outFile << "element face " << triangleNum << std::endl;
     outFile << "property list uchar int vertex_index" << std::endl;
     outFile << "end_header" << std::endl;

     ////
     // Points
     ////

     ///
     // Triangles
     ////

     int fo = 3;
     switch(lod){
     case 0:
         for ( int pi = 0; pi < pointNum; pi++ )
         {
             for ( int vi = 0; vi < 3; ++vi ){
                 float val = mesh.vertices_16[3*pi+vi];
                 outFile.write(reinterpret_cast<char *>(&val ), sizeof( float ) );}
         }
         for ( int ti = 0; ti < triangleNum; ++ti )
         {
             outFile.write(reinterpret_cast<char *>(&fo), sizeof(unsigned char));
             for ( int vi = 0; vi < 3; ++vi )
                 outFile.write(const_cast<char *>(reinterpret_cast<const char *>(&mesh.faces_16[ti*3 + vi])), sizeof( int ) );
         }
         break;
     case 1:
         for ( int pi = 0; pi < pointNum; pi++ )
         {
             for ( int vi = 0; vi < 3; ++vi ){
                 float val = mesh.vertices_32[3*pi+vi];
                 outFile.write(reinterpret_cast<char *>(&val ), sizeof( float ) );}
         }
         for ( int ti = 0; ti < triangleNum; ++ti )
         {
             outFile.write(reinterpret_cast<char *>(&fo), sizeof(unsigned char));
             for ( int vi = 0; vi < 3; ++vi )
                 outFile.write(const_cast<char *>(reinterpret_cast<const char *>(&mesh.faces_32[ti*3 + vi] ) ), sizeof( int ) );
         }
         break;
     case 2:
         for ( int pi = 0; pi < pointNum; pi++ )
         {for ( int vi = 0; vi < 3; ++vi ){
                 float val = mesh.vertices_64[3*pi+vi];
                 outFile.write(reinterpret_cast<char *>(&val ), sizeof( float ) );}
         }
         for ( int ti = 0; ti < triangleNum; ++ti )
         {
             outFile.write(reinterpret_cast<char *>(&fo), sizeof(unsigned char));
             for ( int vi = 0; vi < 3; ++vi )
                 outFile.write(const_cast<char *>(reinterpret_cast<const char *>(&mesh.faces_64[ti*3 + vi] ) ), sizeof( int ) );
         }
         break;
     case 3:
         for ( int pi = 0; pi < pointNum; pi++ )
         {
             for ( int vi = 0; vi < 3; ++vi ){
                 float val = mesh.vertices_128[3*pi+vi];
                 outFile.write(reinterpret_cast<char *>(&val ), sizeof( float ) );}
         }
         for ( int ti = 0; ti < triangleNum; ++ti )
         {
             outFile.write(reinterpret_cast<char *>(&fo), sizeof(unsigned char));
             for ( int vi = 0; vi < 3; ++vi )
                 outFile.write(const_cast<char *>(reinterpret_cast<const char *>(&mesh.faces_128[ti*3 + vi] ) ), sizeof( int ) );
         }
         break;
     case 4:
         for ( int pi = 0; pi < pointNum; pi++ )
         {for ( int vi = 0; vi < 3; ++vi ){
                 float val = mesh.vertices_[3*pi+vi];
                 outFile.write(reinterpret_cast<char *>(&val ), sizeof( float ) );}
         }
         for ( int ti = 0; ti < triangleNum; ++ti )
         {
             outFile.write(reinterpret_cast<char *>(&fo), sizeof(unsigned char));
             for ( int vi = 0; vi < 3; ++vi )
                 outFile.write(const_cast<char *>(reinterpret_cast<const char *>(&mesh.faces_[ti*3 + vi]  )), sizeof( int ) );
         }
         break;

     }

     return true;
    //return false;
}


}  // namespace data_representation
