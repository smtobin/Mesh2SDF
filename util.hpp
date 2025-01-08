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
Vec3r pointSegmentDirection(const Vec3r& x0, const Vec3r& x1, const Vec3r& x2);

/** Returns the minimum distance from a point (x0) to a line segment x1-x2.
 * @param x0 - the point
 * @param x1, x2 - line segment endpoints
*/
Real pointSegmentDistance(const Vec3r& x0, const Vec3r& x1, const Vec3r& x2);

/** Returns the direction from a point (x0) to the closest point on a triangle x1-x2-x3.
 * This is used to determine the gradient during the SDF computation.
 * @param x0 - the point
 * @param x1, x2, x3 - triangle vertices
 */
Vec3r pointTriangleDirection(const Vec3r& x0, const Vec3r& x1, const Vec3r& x2, const Vec3r& x3);

/** Returns the minimum distance from a point (x0) to a triangle x1-x2-x3.
 * @param x0 - the point
 * @param x1, x2, x3 - triangle vertices
 */
Real pointTriangleDistance(const Vec3r& x0, const Vec3r& x1, const Vec3r& x2, const Vec3r& x3);

/** Determins if a 2D point (x0,y0) is inside (or on the border of) a 2D triangle (x1,y1)-(x2,y2)-(x3,y3), and computes the barycentric coordinates (a,b,c).
 * Note: some Realing point tolerance is used for points near the triangle border. This becomes important in the counting of intersections for inside-outside determination.
 * @param x0,y0 - the point
 * @param x1,y1,x2,y2,x3,y3 - triangle vertices
 * @param a,b,c - barycentric coordinates of (x0,y0) w.r.t the triangle
 * @returns true if (x0,y0) is inside the triangle, false otherwise
 */
bool pointInTriangle2D( Real x0, Real y0,
                        Real x1, Real y1, Real x2, Real y2, Real x3, Real y3,
                        Real& a, Real& b, Real& c);

/** Returns the angle between two 3D vectors, between [-pi, pi]. atan2 is used for numerical precision.
 * @param v1,v2 - the two vectors
 */
Real angleBetweenVectors(const Vec3r& v1, const Vec3r& v2);

/** Performs a spherical linear interpolation between two vectors.
 */
Vec3r vectorSlerp(const Vec3r& v1, const Vec3r& v2, Real t);

/** Performs a spherical bilinear interpolation between two vectors.
 */
Vec3r vectorBiSlerp(const Vec3r& v00, const Vec3r& v01,
                              const Vec3r& v10, const Vec3r& v11,
                              Real t0, Real t1);

/** Performs a spherical trilinear interpolation between two vectors.
 * This is used to interpolate between gradient vectors to find the gradient in the middle of a grid cell in the SDF.
 */
Vec3r vectorTriSlerp( const Vec3r& v000, const Vec3r& v100,
                                const Vec3r& v010, const Vec3r& v110,
                                const Vec3r& v001, const Vec3r& v101,
                                const Vec3r& v011, const Vec3r& v111,
                                Real t0, Real t1, Real t2);



// =========== FILE UTILS ==============

/** Formats a Real to be a constant width string. Puts the Real in scientific notation using std::scientific.
 * This is useful for printing the SDF to file - when reading it back we can know exactly how many digits each Real has.
 * @param value - the Real value to format as a string
 * @param precision - the precision (number of decimal places) the Real should have - the length of the string will be precision + 6 chars
 */
std::string formatFloat(Real value, int precision);

/** Parses a Real 3-Vector from a raw char array given the number of chars the Real has and stores it in a Eigen Vector3f.
 * Assumes the start of the char array is the first char of the first Real, and that the Reals have a single space char between them.
 * @param vec - the vector to store the data in
 * @param str - the raw char buffer to read from
 * @param width - the number of chars each Real has
 */
void parseVector3f(Vec3r& vec, const char* str, int width);

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
