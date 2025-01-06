#include "types.hpp"
#include "util.hpp"
#include "MeshSDF.hpp"


#include <iostream>
#include <memory>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Must specify a file name!" << std::endl;
        exit(-1);
    }

    const std::string filename(argv[1]);
    const std::string ext = filename.substr(filename.size() - 4);

    std::unique_ptr<MeshSDF> sdf;

    if (ext == std::string(".obj") || ext == std::string(".stl"))
    {
        const auto [verts, tris] = loadMeshDataFromFile(filename);  // load vertices and triangles from mesh file

        sdf = std::make_unique<MeshSDF>(verts, tris, 128, 6, true); // create SDF from mesh data
        
        sdf->writeToFile("output.sdf"); // write SDF to file
    }
    else if (ext == std::string(".sdf"))    // if user gives .sdf file, load the SDF from it
    {
        sdf = std::make_unique<MeshSDF>(filename);
    }

    // test the SDF by querying a point
    const Eigen::Vector3f p(0.6, 0.2, 0.9);
    const float dist = sdf->evaluate(p);
    const Eigen::Vector3f grad = sdf->gradient(p);

    std::cout << "Signed distance at (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << dist << std::endl;
    std::cout << "Gradient at (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << grad[0] << ", " << grad[1] << ", " << grad[2] << std::endl;
    

    

    
}