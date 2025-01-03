#include "types.hpp"
#include "util.hpp"
#include "MeshSDF.hpp"


#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Must specify a file name!" << std::endl;
        exit(-1);
    }

    const std::string filename(argv[1]);
    const auto [verts, tris] = loadDataFromObj(filename);

    MeshSDF sdf(verts, tris, 256, 10, true);
    const Eigen::Vector3f p(1.02,-0.02,1.01);
    const float dist = sdf.evaluate(p);
    const Eigen::Vector3f grad = sdf.gradient(p);

    std::cout << "Signed distance at (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << dist << std::endl;
    std::cout << "Gradient at (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << grad[0] << ", " << grad[1] << ", " << grad[2] << std::endl;
}