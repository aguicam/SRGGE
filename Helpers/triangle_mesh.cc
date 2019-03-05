// Author: Marc Comino 2018

#include "./triangle_mesh.h"

#include <algorithm>
#include <limits>
#include <math.h>

namespace data_representation {

TriangleMesh::TriangleMesh() { Clear(); }

void TriangleMesh::Clear() {
  vertices_.clear();
  faces_.clear();
  normals_.clear();

  vertices_128.clear();
  vertices_64.clear();
  vertices_32.clear();
  vertices_16.clear();
  faces_128.clear();
  faces_64.clear();
  faces_32.clear();
  faces_16.clear();
  normals_128.clear();
  normals_64.clear();
  normals_32.clear();
  normals_16.clear();

  min_ = Eigen::Vector3f(std::numeric_limits<float>::max(),
                         std::numeric_limits<float>::max(),
                         std::numeric_limits<float>::max());
  max_ = Eigen::Vector3f(std::numeric_limits<float>::lowest(),
                         std::numeric_limits<float>::lowest(),
                         std::numeric_limits<float>::lowest());
}


}  // namespace data_representation


