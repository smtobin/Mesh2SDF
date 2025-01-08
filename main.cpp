#include "types.hpp"
#include "util.hpp"
#include "MeshSDF.hpp"

#include <iostream>

using namespace mesh2sdf;

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cout << "Mesh2SDF - a command line tool for converting closed oriented triangle meshes into grid-based SDFs. Based on SDFGen by Christopher Batty.\n";
        std::cout << "Usage: Mesh2SDF <filename> <num_cells> <padding> <with_gradients>\n\n";
        std::cout << "\t<filename> specifies a .obj or .stl file representing an oriented triangle mesh.\n";
        std::cout << "\t<num_cells> specifies the number of grid cells per side in the grid.\n";
        std::cout << "\t<padding> specifies the number of cells worth of padding (applied on both sides) between the object bounding box and the boundary of the SDF grid. E.g. if num_cells=128 and padding=5, 118 grid cells will span the length, width, and depth of the object.\n";
        std::cout << "\t<with_gradients> (0 or 1) specifies whether or not to compute the gradients of the SDF along with the distances.\n";
        exit(-1);
    }

    const std::string filename(argv[1]);
    const std::string ext = filename.substr(filename.size() - 4);

    const int num_cells = atoi(argv[2]);
    const int padding = atoi(argv[3]);
    const bool with_gradients = atoi(argv[4]);

    assert(num_cells > padding*2);

    if (ext == std::string(".obj") || ext == std::string(".stl"))
    {
        const auto [verts, tris] = loadMeshDataFromFile(filename);  // load vertices and triangles from mesh file

        MeshSDF sdf(verts, tris, num_cells, padding, with_gradients); // create SDF from mesh data
        
        sdf.writeToFile("output.sdf"); // write SDF to file
    }
    else
    {
        std::cerr << "Input file must be .obj or .stl!" << std::endl;
        exit(-1);
    }

    // test the SDF by querying a point
    // const Eigen::Vector3f p(0.6, 0.2, 0.9);
    // const float dist = sdf->evaluate(p);
    // const Eigen::Vector3f grad = sdf->gradient(p);

    // std::cout << "Signed distance at (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << dist << std::endl;
    // std::cout << "Gradient at (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << grad[0] << ", " << grad[1] << ", " << grad[2] << std::endl;
    

    

    
}