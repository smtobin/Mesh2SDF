#include "util.hpp"

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

/** Returns the minimum distance from a point to a line segment
 * @param x0 - the point
 * @param x1 - one end of the line segment
 * @param x2 - the other end of the line segment
 */
Eigen::Vector3f pointSegmentDirection(const Eigen::Vector3f &x0, const Eigen::Vector3f &x1, const Eigen::Vector3f &x2)
{
    const Eigen::Vector3f &dx(x2 - x1);
    float m2 = dx.squaredNorm();
    // find parameter value of closest point on segment
    float s12 = (float)dx.dot(x2 - x0) / m2;
    // cap parameter value to [0,1]
    if (s12 < 0)
        s12 = 0;
    else if (s12 > 1)
        s12 = 1;

    // and find the distance
    const Eigen::Vector3f closest_point_on_line = s12 * x1 + (1 - s12) * x2;
    return (x0 - closest_point_on_line).normalized();
}

float pointSegmentDistance(const Eigen::Vector3f &x0, const Eigen::Vector3f &x1, const Eigen::Vector3f &x2)
{
    const Eigen::Vector3f dx = x2 - x1;
    float m2 = dx.squaredNorm();
    // find parameter value of closest point on segment
    float s12 = (float)dx.dot(x2 - x0) / m2;
    // cap parameter value to [0,1]
    if (s12 < 0)
        s12 = 0;
    else if (s12 > 1)
        s12 = 1;

    // and find the distance
    const Eigen::Vector3f closest_point_on_line = s12 * x1 + (1 - s12) * x2;
    return (x0 - closest_point_on_line).norm();
}

/** Returns the minimum distance from a point to a triangle.
 * @param x0 - the point
 * @param x1 - 1st triangle vertex
 * @param x2 - 2nd triangle vertex
 * @param x3 - 3rd triangle vertex
 */
Eigen::Vector3f pointTriangleDirection(const Eigen::Vector3f &x0, const Eigen::Vector3f &x1, const Eigen::Vector3f &x2, const Eigen::Vector3f &x3)
{
    // first find barycentric coordinates of closest point on infinite plane
    const Eigen::Vector3f x13(x1 - x3), x23(x2 - x3), x03(x0 - x3);
    float m13 = x13.squaredNorm(), m23 = x23.squaredNorm(), d = x13.dot(x23);
    float invdet = 1.f / std::max(m13 * m23 - d * d, 1e-30f);
    float a = x13.dot(x03), b = x23.dot(x03);

    // the barycentric coordinates themselves
    float w23 = invdet * (m23 * a - d * b);
    float w31 = invdet * (m13 * b - d * a);
    float w12 = 1 - w23 - w31;

    if (w23 >= 0 && w31 >= 0 && w12 >= 0) // if we're inside the triangle
    {
        const Eigen::Vector3f closest_point_on_triangle = w23 * x1 + w31 * x2 + w12 * x3;
        const Eigen::Vector3f diff = x0 - closest_point_on_triangle;
        const float dist = diff.norm();
        if (dist < EPS) // check if we're on the plane - if so use the plane normal as the gradient
        {
            const Eigen::Vector3f x21 = x2 - x1;
            const Eigen::Vector3f x31 = x3 - x1;
            const Eigen::Vector3f n = x21.cross(x31).normalized();
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

float pointTriangleDistance(const Eigen::Vector3f &x0, const Eigen::Vector3f &x1, const Eigen::Vector3f &x2, const Eigen::Vector3f &x3)
{
    // first find barycentric coordinates of closest point on infinite plane
    const Eigen::Vector3f x13(x1 - x3), x23(x2 - x3), x03(x0 - x3);
    float m13 = x13.squaredNorm(), m23 = x23.squaredNorm(), d = x13.dot(x23);
    float invdet = 1.f / std::max(m13 * m23 - d * d, 1e-30f);
    float a = x13.dot(x03), b = x23.dot(x03);

    // the barycentric coordinates themselves
    float w23 = invdet * (m23 * a - d * b);
    float w31 = invdet * (m13 * b - d * a);
    float w12 = 1 - w23 - w31;

    if (w23 >= 0 && w31 >= 0 && w12 >= 0) // if we're inside the triangle
    {
        const Eigen::Vector3f closest_point_on_triangle = w23 * x1 + w31 * x2 + w12 * x3;
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
bool pointInTriangle2D(double x0, double y0,
                       double x1, double y1, double x2, double y2, double x3, double y3,
                       double &a, double &b, double &c)
{
    a = ((y2 - y3) * (x0 - x3) + (x3 - x2) * (y0 - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
    b = ((y3 - y1) * (x0 - x3) + (x1 - x3) * (y0 - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
    c = 1 - a - b;

    return (a + EPS) >= 0 && a <= 1 && (b + EPS) >= 0 && b <= 1 && (c + EPS) >= 0 && c <= 1;
}


float angleBetweenVectors(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2)
{
    const float dot = v1.dot(v2);
    const float det = v1.cross(v2).norm();
    return std::atan2(det, dot);
}

Eigen::Vector3f vectorSlerp(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, float t)
{
    const float angle = angleBetweenVectors(v1, v2);
    // in the unlikely case that v1 and v2 are colinear (i.e. have angle between them of 180 deg)
    if (M_PI - std::abs(angle) < 1e-6)
    {
        // from https://math.stackexchange.com/a/211195 - find perpendicular vector
        Eigen::Vector3f perp_v(v1[2], v1[2], -v1[0]-v1[1]);
        if (perp_v.squaredNorm() < EPS) perp_v = Eigen::Vector3f(-v1[1]-v1[2], v1[0], v1[0]);
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

    const float inv_sin_ang = 1 / std::sin(angle);
    return std::sin( (1-t) * angle) * inv_sin_ang * v1 + std::sin(t*angle) * inv_sin_ang * v2;
}

Eigen::Vector3f vectorBiSlerp(const Eigen::Vector3f& v00, const Eigen::Vector3f& v10,
                              const Eigen::Vector3f& v01, const Eigen::Vector3f& v11,
                              float t0, float t1)
{
    return vectorSlerp( vectorSlerp(v00, v10, t0),
                        vectorSlerp(v01, v11, t0),
                        t1   );
}

Eigen::Vector3f vectorTriSlerp( const Eigen::Vector3f& v000, const Eigen::Vector3f& v100,
                                const Eigen::Vector3f& v010, const Eigen::Vector3f& v110,
                                const Eigen::Vector3f& v001, const Eigen::Vector3f& v101,
                                const Eigen::Vector3f& v011, const Eigen::Vector3f& v111,
                                float t0, float t1, float t2)
{
    return vectorSlerp( vectorBiSlerp(v000, v100, v010, v110, t0, t1),
                        vectorBiSlerp(v001, v101, v011, v111, t0, t1),
                        t2  );
}

VerticesAndTriangles loadDataFromObj(const std::string &filename)
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
    std::vector<std::array<float, 3>> vert_list;
    std::vector<std::array<int, 3>> tri_list;
    while (!infile.eof())
    {
        std::getline(infile, line);

        //.obj files sometimes contain vertex normals indicated by "vn"
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn"))
        {
            std::stringstream data(line);
            char c;
            std::array<float, 3> point;
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
        verts.col(i) = Eigen::Vector3f(vert_list[i][0], vert_list[i][1], vert_list[i][2]);
    }

    for (size_t i = 0; i < tri_list.size(); i++)
    {
        tris.col(i) = Eigen::Vector3i(tri_list[i][0] - 1, tri_list[i][1] - 1, tri_list[i][2] - 1);
    }

    std::cout << "Loaded " << verts.cols() << " vertices and " << tris.cols() << " triangles from " << filename << "." << std::endl;

    return std::make_pair(verts, tris);
}

std::string formatFloat(float value, int width)
{
    std::ostringstream oss;
    if (value == 0.0f) value = 0.0f;    // fix -0.0f
    if (value >= 0)
    {
        oss << std::scientific << std::setprecision(width) << value;
    }
    else
    {
        oss << std::scientific << std::setprecision(width-1) << value;
    }
    return oss.str();
}

void parseVector3f(Eigen::Vector3f& vec, const char* str, int width)
{
    vec[0] = atof(str);
    vec[1] = atof(str+width+1);
    vec[2] = atof(str+2*width+2);
}