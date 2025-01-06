#ifndef __UTIL_HPP
#define __UTIL_HPP

#include "types.hpp"

Eigen::Vector3f pointSegmentDirection(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2);
float pointSegmentDistance(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2);
Eigen::Vector3f pointTriangleDirection(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, const Eigen::Vector3f& x3);
float pointTriangleDistance(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, const Eigen::Vector3f& x3);
bool pointInTriangle2D( double x0, double y0,
                        double x1, double y1, double x2, double y2, double x3, double y3,
                        double& a, double& b, double& c);

float angleBetweenVectors(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2);
Eigen::Vector3f vectorSlerp(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, float t);
Eigen::Vector3f vectorBiSlerp(const Eigen::Vector3f& v00, const Eigen::Vector3f& v01,
                              const Eigen::Vector3f& v10, const Eigen::Vector3f& v11,
                              float t0, float t1);
Eigen::Vector3f vectorTriSlerp( const Eigen::Vector3f& v000, const Eigen::Vector3f& v100,
                                const Eigen::Vector3f& v010, const Eigen::Vector3f& v110,
                                const Eigen::Vector3f& v001, const Eigen::Vector3f& v101,
                                const Eigen::Vector3f& v011, const Eigen::Vector3f& v111,
                                float t0, float t1, float t2);

VerticesAndTriangles loadDataFromObj(const std::string& filename);

std::string formatFloat(float value, int width);
void parseVector3f(Eigen::Vector3f& vec, const char* str, int width);

#endif
