// Author: Marc Comino 2018

#ifndef MESH_IO_H_
#define MESH_IO_H_

#include <Helpers/triangle_mesh.h>

#include <string>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen>

namespace data_representation {
struct ce_vert{
    int cell;
    std::vector<int> inices_vertex;
};


/**
 * @brief ReadFromPly Read the mesh stored in PLY format at the path filename
 * and stores the corresponding TriangleMesh representation
 * @param filename The path to the PLY mesh.
 * @param mesh The resulting representation with computed per-vertex normals.
 * @return Whether it was able to read the file.
 */
bool ReadFromPly(const std::string &filename, TriangleMesh *mesh);
bool ReadFromPlyMuseum0(const std::string &filename, TriangleMesh *mesh);
bool ReadFromPlyMuseum1(const std::string &filename, TriangleMesh *mesh);
bool ReadFromPlyMuseum2(const std::string &filename, TriangleMesh *mesh);
bool ReadFromPlyMuseum3(const std::string &filename, TriangleMesh *mesh);
bool ReadFromPlyMuseum4(const std::string &filename, TriangleMesh *mesh);

/**
 * @brief WriteToPly Stores the mesh representation in PLY format at the path
 * filename.
 * @param filename The path where the mesh will be stored.
 * @param mesh The mesh to be stored.
 * @return Whether it was able to store the file.
 */
bool WriteToPly(const std::string &filename, TriangleMesh &mesh, const int lod);


}  // namespace data_representation





#endif  // MESH_IO_H_
