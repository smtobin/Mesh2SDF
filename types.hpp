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

#ifdef DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

typedef Array3<Real> Array3r;
typedef Array3<int> Array3i;
typedef Eigen::Vector<Real, 3> Vec3r;
typedef Array3<Vec3r> Array3Vec3r;
typedef Eigen::Matrix<Real, 3, -1> VertexMat;
typedef Eigen::Matrix<int, 3, -1> TriangleMat;

typedef std::pair<VertexMat,TriangleMat> VerticesAndTriangles;

}

#endif // __TYPES_HPP