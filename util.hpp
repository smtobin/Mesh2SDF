#ifndef __MESH2SDF_UTIL_HPP
#define __MESH2SDF_UTIL_HPP

#include "types.hpp"

namespace mesh2sdf
{

// =========== GEOMETRY UTILS ==============

/** Returns the direction from a point (x0) to the closest point on a line segment x1-x2.
 * This is used to determine the gradient during the SDF computation.
 * @param x0 - the point
 * @param x1, x2 - line segment endpoints
 */
Eigen::Vector3f pointSegmentDirection(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2);

/** Returns the minimum distance from a point (x0) to a line segment x1-x2.
 * @param x0 - the point
 * @param x1, x2 - line segment endpoints
*/
float pointSegmentDistance(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2);

/** Returns the direction from a point (x0) to the closest point on a triangle x1-x2-x3.
 * This is used to determine the gradient during the SDF computation.
 * @param x0 - the point
 * @param x1, x2, x3 - triangle vertices
 */
Eigen::Vector3f pointTriangleDirection(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, const Eigen::Vector3f& x3);

/** Returns the minimum distance from a point (x0) to a triangle x1-x2-x3.
 * @param x0 - the point
 * @param x1, x2, x3 - triangle vertices
 */
float pointTriangleDistance(const Eigen::Vector3f& x0, const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, const Eigen::Vector3f& x3);

/** Determins if a 2D point (x0,y0) is inside (or on the border of) a 2D triangle (x1,y1)-(x2,y2)-(x3,y3), and computes the barycentric coordinates (a,b,c).
 * Note: some floating point tolerance is used for points near the triangle border. This becomes important in the counting of intersections for inside-outside determination.
 * @param x0,y0 - the point
 * @param x1,y1,x2,y2,x3,y3 - triangle vertices
 * @param a,b,c - barycentric coordinates of (x0,y0) w.r.t the triangle
 * @returns true if (x0,y0) is inside the triangle, false otherwise
 */
bool pointInTriangle2D( double x0, double y0,
                        double x1, double y1, double x2, double y2, double x3, double y3,
                        double& a, double& b, double& c);

/** Returns the angle between two 3D vectors, between [-pi, pi]. atan2 is used for numerical precision.
 * @param v1,v2 - the two vectors
 */
float angleBetweenVectors(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2);

/** Performs a spherical linear interpolation between two vectors.
 */
Eigen::Vector3f vectorSlerp(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, float t);

/** Performs a spherical bilinear interpolation between two vectors.
 */
Eigen::Vector3f vectorBiSlerp(const Eigen::Vector3f& v00, const Eigen::Vector3f& v01,
                              const Eigen::Vector3f& v10, const Eigen::Vector3f& v11,
                              float t0, float t1);

/** Performs a spherical trilinear interpolation between two vectors.
 * This is used to interpolate between gradient vectors to find the gradient in the middle of a grid cell in the SDF.
 */
Eigen::Vector3f vectorTriSlerp( const Eigen::Vector3f& v000, const Eigen::Vector3f& v100,
                                const Eigen::Vector3f& v010, const Eigen::Vector3f& v110,
                                const Eigen::Vector3f& v001, const Eigen::Vector3f& v101,
                                const Eigen::Vector3f& v011, const Eigen::Vector3f& v111,
                                float t0, float t1, float t2);



// =========== FILE UTILS ==============

/** Formats a float to be a constant width string. Puts the float in scientific notation using std::scientific.
 * This is useful for printing the SDF to file - when reading it back we can know exactly how many digits each float has.
 * @param value - the float value to format as a string
 * @param precision - the precision (number of decimal places) the float should have - the length of the string will be precision + 6 chars
 */
std::string formatFloat(float value, int precision);

/** Parses a float 3-Vector from a raw char array given the number of chars the float has and stores it in a Eigen Vector3f.
 * Assumes the start of the char array is the first char of the first float, and that the floats have a single space char between them.
 * @param vec - the vector to store the data in
 * @param str - the raw char buffer to read from
 * @param width - the number of chars each float has
 */
void parseVector3f(Eigen::Vector3f& vec, const char* str, int width);

#ifdef HAVE_ASSIMP

VerticesAndTriangles loadMeshDataFromFile(const std::string& filename);

#else

/** A simple parser to load vertices and triangles from a .obj file.
 * Very simple - cannot handle vertex normals or materials or anything fancy - just vertices and triangles.
 */
VerticesAndTriangles loadMeshDataFromFile(const std::string& filename);

#endif


} // namespace mesh2sdf

#endif
