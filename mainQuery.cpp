#include "types.hpp"
#include "util.hpp"
#include "MeshSDF.hpp"

#include <iostream>

using namespace mesh2sdf;

int main(int argc, char* argv[])
{
    // this is more just meant to be a sort of playground for trying out the API
    const auto [verts, tris] = loadMeshDataFromFile("obj/cube.obj");    // cube.obj is a cube spanning from (0,0,0) to (1,1,1)
    MeshSDF sdf1(verts, tris, 128, 8, true);
    sdf1.writeToFile("output.sdf");

    // load the SDF we just saved back to make sure that it is the same as the original we just created
    MeshSDF sdf2("output.sdf");

    BoundingBox grid_bbox = sdf2.gridBoundingBox();
    BoundingBox mesh_bbox = sdf2.meshBoundingBox();
    std::cout << "Grid bbox: (" << grid_bbox.first[0] << ", " << grid_bbox.first[1] << ", " << grid_bbox.first[2] << ") to ("
        << grid_bbox.second[0] << ", " << grid_bbox.second[1] << ", " << grid_bbox.second[2] << ")" << std::endl;
    std::cout << "Mesh bbox: (" << mesh_bbox.first[0] << ", " << mesh_bbox.first[1] << ", " << mesh_bbox.first[2] << ") to ("
        << mesh_bbox.second[0] << ", " << mesh_bbox.second[1] << ", " << mesh_bbox.second[2] << ")" << std::endl;

    // test the SDFs by querying a point inside the cube
    const Vec3r p(0.6, 0.2, 0.9);
    const Real dist1 = sdf1.evaluate(p);
    const Vec3r grad1 = sdf1.gradient(p);

    const Real dist2 = sdf2.evaluate(p);
    const Vec3r grad2 = sdf2.gradient(p);

    // they should be the same
    std::cout << "Signed distance at (" << p[0] << ", " << p[1] << ", " << p[2] << "):\n  Original SDF: " << dist1 << "\n  Reloaded SDF: " << dist2 << std::endl;
    std::cout << "Gradient at (" << p[0] << ", " << p[1] << ", " << p[2] << "):\n  Original SDF: " << grad1[0] << ", " << grad1[1] << ", " << grad1[2] <<
        "\n  Reloaded SDF: " << grad2[0] << ", " << grad2[1] << ", " << grad2[2] << std::endl;
    

    

    
}