#include "util.hpp"

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace mesh2sdf
{

/** Returns the minimum distance from a point to a line segment
 * @param x0 - the point
 * @param x1 - one end of the line segment
 * @param x2 - the other end of the line segment
 */
Vec3r pointSegmentDirection(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2)
{
    const Vec3r &dx(x2 - x1);
    Real m2 = dx.squaredNorm();
    // find parameter value of closest point on segment
    Real s12 = (Real)dx.dot(x2 - x0) / m2;
    // cap parameter value to [0,1]
    if (s12 < 0)
        s12 = 0;
    else if (s12 > 1)
        s12 = 1;

    // and find the distance
    const Vec3r closest_point_on_line = s12 * x1 + (1 - s12) * x2;
    return (x0 - closest_point_on_line).normalized();
}

Real pointSegmentDistance(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2)
{
    const Vec3r dx = x2 - x1;
    Real m2 = dx.squaredNorm();
    // find parameter value of closest point on segment
    Real s12 = (Real)dx.dot(x2 - x0) / m2;
    // cap parameter value to [0,1]
    if (s12 < 0)
        s12 = 0;
    else if (s12 > 1)
        s12 = 1;

    // and find the distance
    const Vec3r closest_point_on_line = s12 * x1 + (1 - s12) * x2;
    return (x0 - closest_point_on_line).norm();
}

/** Returns the minimum distance from a point to a triangle.
 * @param x0 - the point
 * @param x1 - 1st triangle vertex
 * @param x2 - 2nd triangle vertex
 * @param x3 - 3rd triangle vertex
 */
Vec3r pointTriangleDirection(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2, const Vec3r &x3)
{
    // first find barycentric coordinates of closest point on infinite plane
    const Vec3r x13(x1 - x3), x23(x2 - x3), x03(x0 - x3);
    Real m13 = x13.squaredNorm(), m23 = x23.squaredNorm(), d = x13.dot(x23);
    Real invdet = 1.f / std::max(m13 * m23 - d * d, Real(1e-30));
    Real a = x13.dot(x03), b = x23.dot(x03);

    // the barycentric coordinates themselves
    Real w23 = invdet * (m23 * a - d * b);
    Real w31 = invdet * (m13 * b - d * a);
    Real w12 = 1 - w23 - w31;

    if (w23 >= 0 && w31 >= 0 && w12 >= 0) // if we're inside the triangle
    {
        const Vec3r closest_point_on_triangle = w23 * x1 + w31 * x2 + w12 * x3;
        const Vec3r diff = x0 - closest_point_on_triangle;
        const Real dist = diff.norm();
        if (dist < MESH2SDF_EPS) // check if we're on the plane - if so use the plane normal as the gradient
        {
            const Vec3r x21 = x2 - x1;
            const Vec3r x31 = x3 - x1;
            const Vec3r n = x21.cross(x31).normalized();
            return n;
        }
        return diff / dist;
    }
    else // we have to clamp to one of the edges
    {
        if (w23 > 0) // this rules out edge 2-3 for us
            if (pointSegmentDistance(x0, x1, x2) < pointSegmentDistance(x0, x1, x3))
                return pointSegmentDirection(x0, x1, x2);
            else
                return pointSegmentDirection(x0, x1, x3);
        else if (w31 > 0) // this rules out edge 1-3
            if (pointSegmentDistance(x0,x1,x2) < pointSegmentDistance(x0,x2,x3))
                return pointSegmentDirection(x0,x1,x2);
            else
                return pointSegmentDirection(x0,x2,x3);
        else // w12 must be >0, ruling out edge 1-2
            if (pointSegmentDistance(x0,x1,x3) < pointSegmentDistance(x0,x2,x3))
                return pointSegmentDirection(x0,x1,x3);
            else
                return pointSegmentDirection(x0,x2,x3);
    }
}

Real pointTriangleDistance(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2, const Vec3r &x3)
{
    // first find barycentric coordinates of closest point on infinite plane
    const Vec3r x13(x1 - x3), x23(x2 - x3), x03(x0 - x3);
    Real m13 = x13.squaredNorm(), m23 = x23.squaredNorm(), d = x13.dot(x23);
    Real invdet = 1.f / std::max(m13 * m23 - d * d, Real(1e-30));
    Real a = x13.dot(x03), b = x23.dot(x03);

    // the barycentric coordinates themselves
    Real w23 = invdet * (m23 * a - d * b);
    Real w31 = invdet * (m13 * b - d * a);
    Real w12 = 1 - w23 - w31;

    if (w23 >= 0 && w31 >= 0 && w12 >= 0) // if we're inside the triangle
    {
        const Vec3r closest_point_on_triangle = w23 * x1 + w31 * x2 + w12 * x3;
        return (x0 - closest_point_on_triangle).norm();
    }
    else // we have to clamp to one of the edges
    {
        if (w23 > 0) // this rules out edge 2-3 for us
            return std::min(pointSegmentDistance(x0, x1, x2), pointSegmentDistance(x0, x1, x3));
        else if (w31 > 0) // this rules out edge 1-3
            return std::min(pointSegmentDistance(x0, x1, x2), pointSegmentDistance(x0, x2, x3));
        else // w12 must be >0, ruling out edge 1-2
            return std::min(pointSegmentDistance(x0, x1, x3), pointSegmentDistance(x0, x2, x3));
    }
}

/** Returns true if point (x0,y0) is in the triangle defined by (x1,y1), (x2,y2), (x3,y3).
 * The barycentric coordinates (a,b,c) are computed.
 * @param x0,y0 - the point
 * @param x1,y1 - the 1st triangle vertex
 * @param x2,y2 - the 2nd triangle vertex
 * @param x3,y3 - the 3rd triangle vertex
 * @param a,b,c (OUTPUT) - the barycentric coordinates of (x0,y0) in the triangle.
 */
bool pointInTriangle2D(Real x0, Real y0,
                       Real x1, Real y1, Real x2, Real y2, Real x3, Real y3,
                       Real &a, Real &b, Real &c)
{
    a = ((y2 - y3) * (x0 - x3) + (x3 - x2) * (y0 - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
    b = ((y3 - y1) * (x0 - x3) + (x1 - x3) * (y0 - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
    c = 1 - a - b;

    return (a + MESH2SDF_EPS) >= 0 && a <= 1 && (b + MESH2SDF_EPS) >= 0 && b <= 1 && (c + MESH2SDF_EPS) >= 0 && c <= 1;
}


Real angleBetweenVectors(const Vec3r& v1, const Vec3r& v2)
{
    const Real dot = v1.dot(v2);
    const Real det = v1.cross(v2).norm();
    return std::atan2(det, dot);
}

Vec3r vectorSlerp(const Vec3r& v1, const Vec3r& v2, Real t)
{
    const Real angle = angleBetweenVectors(v1, v2);
    // in the unlikely case that v1 and v2 are colinear (i.e. have angle between them of 180 deg)
    if (M_PI - std::abs(angle) < 1e-6)
    {
        // from https://math.stackexchange.com/a/211195 - find perpendicular vector
        Vec3r perp_v(v1[2], v1[2], -v1[0]-v1[1]);
        if (perp_v.squaredNorm() < MESH2SDF_EPS) perp_v = Vec3r(-v1[1]-v1[2], v1[0], v1[0]);
        perp_v.normalize();

        // this perpendicular vector is the halfway point
        if (t < 0.5f)   return vectorSlerp(v1, perp_v, t*2.0);
        else            return vectorSlerp(perp_v, v2, (t-0.5)*2.0);
    }
    // if angle is 0 between them, just return v1
    else if (std::abs(angle) < 1e-6)
    {
        return v1;
    }

    const Real inv_sin_ang = 1 / std::sin(angle);
    return std::sin( (1-t) * angle) * inv_sin_ang * v1 + std::sin(t*angle) * inv_sin_ang * v2;
}

Vec3r vectorBiSlerp(const Vec3r& v00, const Vec3r& v10,
                              const Vec3r& v01, const Vec3r& v11,
                              Real t0, Real t1)
{
    return vectorSlerp( vectorSlerp(v00, v10, t0),
                        vectorSlerp(v01, v11, t0),
                        t1   );
}

Vec3r vectorTriSlerp( const Vec3r& v000, const Vec3r& v100,
                                const Vec3r& v010, const Vec3r& v110,
                                const Vec3r& v001, const Vec3r& v101,
                                const Vec3r& v011, const Vec3r& v111,
                                Real t0, Real t1, Real t2)
{
    return vectorSlerp( vectorBiSlerp(v000, v100, v010, v110, t0, t1),
                        vectorBiSlerp(v001, v101, v011, v111, t0, t1),
                        t2  );
}


std::string formatFloat(Real value, int precision)
{
    std::ostringstream oss;
    if (value == 0.0f) value = 0.0f;    // fix -0.0f
    if (value >= 0)
    {
        oss << std::scientific << std::setprecision(precision) << value;
    }
    else
    {
        oss << std::scientific << std::setprecision(precision-1) << value;
    }
    return oss.str();
}

void parseVector3f(Vec3r& vec, const char* str, int width)
{
    vec[0] = atof(str);
    vec[1] = atof(str+width+1);
    vec[2] = atof(str+2*width+2);
}

#ifdef HAVE_ASSIMP

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags

VerticesAndTriangles loadMeshDataFromFile(const std::string& filename)
{
    Assimp::Importer importer;

    // And have it read the given file with some example postprocessing
    // Usually - if speed is not the most important aspect for you - you'll
    // probably to request more postprocessing than we do in this example.
    const aiScene* scene = importer.ReadFile( filename,
        aiProcess_Triangulate            |
        aiProcess_JoinIdenticalVertices  |
        aiProcess_SortByPType);

    // If the import failed, report it
    if (scene == nullptr)
    {
        std::cerr << "\tAssimp::Importer could not open " << filename << std::endl;
        std::cerr << "\tEnsure that the file is in a format that assimp can handle." << std::endl;
        assert(0);
    }

    const aiMesh* ai_mesh = scene->mMeshes[0];

    // Extract vertices
    VertexMat verts(3, ai_mesh->mNumVertices);
    for (unsigned i = 0; i < ai_mesh->mNumVertices; i++)
    {
        verts(0,i) = ai_mesh->mVertices[i].x;
        verts(1,i) = ai_mesh->mVertices[i].y;
        verts(2,i) = ai_mesh->mVertices[i].z;
    }

    // Extract faces
    TriangleMat tris(3, ai_mesh->mNumFaces);
    for (unsigned i = 0; i < ai_mesh->mNumFaces; i++)
    {
        tris(0,i) = ai_mesh->mFaces[i].mIndices[0];
        tris(1,i) = ai_mesh->mFaces[i].mIndices[1];
        tris(2,i) = ai_mesh->mFaces[i].mIndices[2];
    }

    std::cout << "Loaded " << verts.cols() << " vertices and " << tris.cols() << " triangles from " << filename << "." << std::endl;

    return std::make_pair(verts, tris);
}

#else

VerticesAndTriangles loadMeshDataFromFile(const std::string &filename)
{
    // make sure .obj file was passed in
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj." << std::endl;
        assert(0);
    }

    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Failed to open file with name " << filename << ". Terminating." << std::endl;
        assert(0);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<std::array<Real, 3>> vert_list;
    std::vector<std::array<int, 3>> tri_list;
    while (!infile.eof())
    {
        std::getline(infile, line);

        //.obj files sometimes contain vertex normals indicated by "vn"
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn"))
        {
            std::stringstream data(line);
            char c;
            std::array<Real, 3> point;
            data >> c >> point[0] >> point[1] >> point[2];
            vert_list.push_back(std::move(point));
        }
        else if (line.substr(0, 1) == std::string("f"))
        {
            std::stringstream data(line);
            char c;
            std::array<int, 3> tri;
            data >> c >> tri[0] >> tri[1] >> tri[2];
            tri_list.push_back(std::move(tri));
        }
        else if (line.substr(0, 2) == std::string("vn"))
        {
            std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
            assert(0);
        }
        else
        {
            ++ignored_lines;
        }
    }
    infile.close();

    // create Eigen matrices from vert and tri lists
    VertexMat verts(3, vert_list.size());
    TriangleMat tris(3, tri_list.size());
    for (size_t i = 0; i < vert_list.size(); i++)
    {
        verts.col(i) = Vec3r(vert_list[i][0], vert_list[i][1], vert_list[i][2]);
    }

    for (size_t i = 0; i < tri_list.size(); i++)
    {
        tris.col(i) = Eigen::Vector3i(tri_list[i][0] - 1, tri_list[i][1] - 1, tri_list[i][2] - 1);
    }

    std::cout << "Loaded " << verts.cols() << " vertices and " << tris.cols() << " triangles from " << filename << "." << std::endl;

    return std::make_pair(verts, tris);
}

} // namespace mesh2sdf

#endif
