#ifndef __MESH2SDF_MESH_SDF_HPP
#define __MESH2SDF_MESH_SDF_HPP

#include "types.hpp"
#include "array3.hpp"

namespace mesh2sdf
{

class MeshSDF
{
    public:
    /** Construct a SDF from input vertices and triangles.
     * @param verts - the vertices of the mesh
     * @param tris - the triangles of the mesh, specified as columns of 3 integers
     * @param grid_size - the number of cells per side in the grid
     * @param padding - the number of cells of padding around the mesh - by default 1
     * @param with_gradient - whether or not to store and evaluate the gradient when computing the SDF - by default true
     */
    MeshSDF(const VertexMat& verts, const TriangleMat& tris, int grid_size, int padding=1, bool with_gradient = true);

    /** Loads a SDF that was written to file by this class. Expects a .sdf extension. */
    MeshSDF(const std::string& filename);

    /** Evaluates the SDF at the given point. */
    float evaluate(const Eigen::Vector3f& p) const;
    
    /** Computes the gradient of the SDF at the given point. */
    Eigen::Vector3f gradient(const Eigen::Vector3f& p) const;

    /** Writes the SDF to specified file. */
    void writeToFile(const std::string& filename) const;

    private:

    /** Helper function that loads an SDF from file. Expects a .sdf file that was produced by this class. */
    void _loadSDFFromFile(const std::string& filename);

    /** Computes the SDF given vertices and faces.
     * Adapted from Christopher Batty's code.
     */
    void _makeSDF(const VertexMat& verts, const TriangleMat& tris, int padding, bool with_gradient);

    /** Helper function for SDF computation.
     * Adapted from Christopher Batty's code.
     */
    void _checkNeighbor(const VertexMat& verts, const TriangleMat& tris, 
                        Array3i& closest_tri,
                        const Eigen::Vector3f& grid_point,
                        int i0, int j0, int k0, int i1, int j1, int k1);

    /** Helper function for SDF computation.
     * Adapted from Christopher Batty's code.
     */
    void _sweep(const VertexMat& verts, const TriangleMat& tris,
                Array3i& closest_tri,
                int di, int dj, int dk);

    /** Returns the global coordinates of a point in the SDF grid given (i,j,k) grid indices. */
    Eigen::Vector3f _gridPointFromIJK(int i, int j, int k) const;

    /** Returns DECIMAL (i,j,k) grid indices given a point in global coordinates. */
    Eigen::Vector3f _gridIJKFromPoint(const Eigen::Vector3f& p) const;

    private:
    static constexpr int FLOAT_PRECISION = 10;  // the number of places after the decimal to use when printing the SDF to file

    int _N;  // number of voxels per side in the grid
    Eigen::Vector3f _cell_size;  // size of each voxel in the grid

    Eigen::Vector3f _bbox_min;  // bounding box minimum for the SDF grid
    Eigen::Vector3f _bbox_max;  // bounding box maximum for the SDF grid

    bool _with_gradient; // whether or not the gradient of the SDF was computed and stored along with the distance

    Array3<float> _distance_grid;  // stores the distances in a grid
    Array3<Eigen::Vector3f> _gradient_grid;  // stores the gradients in a grid - this will be empty if with_gradient=false
};

} // namespace mesh2sdf

#endif // __MESH_SDF_HPP