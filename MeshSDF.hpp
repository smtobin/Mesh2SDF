#ifndef __MESH2SDF_MESH_SDF_HPP
#define __MESH2SDF_MESH_SDF_HPP

#include "types.hpp"
#include "array3.hpp"

namespace mesh2sdf
{

class MeshSDF
{
    public:
    /** Empty constructor */
    MeshSDF();

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
    Real evaluate(const Vec3r& p) const;
    
    /** Computes the gradient of the SDF at the given point. */
    Vec3r gradient(const Vec3r& p) const;

    /** Writes the SDF to specified file. */
    void writeToFile(const std::string& filename) const;

    /** Returns the bounding box around the SDF grid. */
    BoundingBox gridBoundingBox() const { return BoundingBox(_grid_bbox_min, _grid_bbox_max); }

    /** Returns the bounding box around the mesh itself. */
    BoundingBox meshBoundingBox() const { return BoundingBox(_mesh_bbox_min, _mesh_bbox_max); }

    /** Returns the mass center for the mesh. */
    Vec3r meshMassCenter() const { return _mesh_mass_center; }

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
                        const Vec3r& grid_point,
                        int i0, int j0, int k0, int i1, int j1, int k1);

    /** Helper function for SDF computation.
     * Adapted from Christopher Batty's code.
     */
    void _sweep(const VertexMat& verts, const TriangleMat& tris,
                Array3i& closest_tri,
                int di, int dj, int dk);

    /** Returns the global coordinates of a point in the SDF grid given integer (i,j,k) grid indices. */
    Vec3r _gridPointFromIJK(int i, int j, int k) const;

    /** Returns the global coordinates of a point in the SDF grid given decimal (i,j,k) grid indices. */
    Vec3r _gridPointFromIJK(Real i, Real j, Real k) const;

    /** Returns DECIMAL (i,j,k) grid indices given a point in global coordinates. */
    Vec3r _gridIJKFromPoint(const Vec3r& p) const;

    /** Performs trilinear interpolation on the distance grid for a cell with corners (i0,j0,k0) to (i1,j1,k1).
     * (id,jd,kd) are each in [0,1] and are the interpolation parameter.
     * Returns the interpolated distance in the middle of the cell.
     */
    Real _interpolateDistanceGrid(int i0, int j0, int k0, int i1, int j1, int k1, int id, int jd, int kd) const;

    /** Performs trilinear spherical interpolation on the gradient grid for a cell with corners (i0,j0,k0) to (i1,j1,k1).
     * (id,jd,kd) are each in [0,1] and are the interpolation parameters.
     * Reteurns the interpolated gradient in the middle of the cell.
     */
    Vec3r _interpolateGradientGrid(int i0, int j0, int k0, int i1, int j1, int k1, int id, int jd, int kd) const;

    private:
    static constexpr int FLOAT_PRECISION = 10;  // the number of places after the decimal to use when printing the SDF to file

    int _N;  // number of voxels per side in the grid
    Vec3r _cell_size;  // size of each voxel in the grid

    Vec3r _grid_bbox_min;  // bounding box minimum for the SDF grid
    Vec3r _grid_bbox_max;  // bounding box maximum for the SDF grid

    Vec3r _mesh_bbox_min;   // bounding box minimum for the mesh itself - this will likely be larger than the grid bbox min due to padding cells
    Vec3r _mesh_bbox_max;   // bounding box maximum for the mesh itself - this will likely be smaller than the grid bbox max due to padding cells

    Vec3r _mesh_mass_center; // the center of mass for the mesh

    bool _with_gradient; // whether or not the gradient of the SDF was computed and stored along with the distance

    Array3<Real> _distance_grid;  // stores the distances in a grid
    Array3<Vec3r> _gradient_grid;  // stores the gradients in a grid - this will be empty if with_gradient=false
};

} // namespace mesh2sdf

#endif // __MESH_SDF_HPP