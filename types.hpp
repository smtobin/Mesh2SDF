#ifndef __MESH2SDF_TYPES_HPP
#define __MESH2SDF_TYPES_HPP

#include <Eigen/Dense>
#include "array3.hpp"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define MESH2SDF_EPS 1e-8

namespace mesh2sdf
{

typedef Array3<float> Array3f;
typedef Array3<int> Array3i;
typedef Array3<Eigen::Vector3f> Array3Vec3f;
typedef Eigen::Matrix<float, 3, -1> VertexMat;
typedef Eigen::Matrix<int, 3, -1> TriangleMat;

typedef std::pair<VertexMat,TriangleMat> VerticesAndTriangles;

}

#endif // __TYPES_HPP