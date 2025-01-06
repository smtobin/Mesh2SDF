#include "types.hpp"
#include "util.hpp"
#include "MeshSDF.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    // this is more just meant to be a sort of playground for trying out the API
    const auto [verts, tris] = loadMeshDataFromFile("obj/cube.obj");    // cube.obj is a cube spanning from (0,0,0) to (1,1,1)
    MeshSDF sdf1(verts, tris, 128, 8, true);
    sdf1.writeToFile("output.sdf");

    // load the SDF we just saved back to make sure that it is the same as the original we just created
    MeshSDF sdf2("output.sdf");

    // test the SDFs by querying a point inside the cube
    const Eigen::Vector3f p(0.6, 0.2, 0.9);
    const float dist1 = sdf1.evaluate(p);
    const Eigen::Vector3f grad1 = sdf1.gradient(p);

    const float dist2 = sdf2.evaluate(p);
    const Eigen::Vector3f grad2 = sdf2.gradient(p);

    // they should be the same
    std::cout << "Signed distance at (" << p[0] << ", " << p[1] << ", " << p[2] << "):\n  Original SDF: " << dist1 << "\n  Reloaded SDF: " << dist2 << std::endl;
    std::cout << "Gradient at (" << p[0] << ", " << p[1] << ", " << p[2] << "):\n  Original SDF: " << grad1[0] << ", " << grad1[1] << ", " << grad1[2] <<
        "\n  Reloaded SDF: " << grad2[0] << ", " << grad2[1] << ", " << grad2[2] << std::endl;
    

    

    
}