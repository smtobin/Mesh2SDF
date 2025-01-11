#include "MeshSDF.hpp"
#include "util.hpp"

#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>

// for fast file parsing
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

namespace mesh2sdf
{

MeshSDF::MeshSDF()
    : _N(0), _cell_size(),
     _grid_bbox_min(Vec3r::Zero()), _grid_bbox_max(Vec3r::Zero()),
     _mesh_bbox_min(Vec3r::Zero()), _mesh_bbox_max(Vec3r::Zero()),
     _with_gradient(false)
{}

MeshSDF::MeshSDF(const VertexMat& verts, const TriangleMat& tris, int grid_size, int padding, bool with_gradient)
    : _N(grid_size), _with_gradient(with_gradient)
{
    _makeSDF(verts, tris, padding, with_gradient);
}

MeshSDF::MeshSDF(const std::string& filename)
    : _N(0), _cell_size(),
     _grid_bbox_min(Vec3r::Zero()), _grid_bbox_max(Vec3r::Zero()),
     _mesh_bbox_min(Vec3r::Zero()), _mesh_bbox_max(Vec3r::Zero()),
     _with_gradient(false)
{
    _loadSDFFromFile(filename);
}

Real MeshSDF::evaluate(const Vec3r& p) const
{
    const Vec3r& ijk = _gridIJKFromPoint(p);
    
    const int i0 = std::floor(ijk[0]);    const int i1 = i0+1;
    const int j0 = std::floor(ijk[1]);    const int j1 = j0+1;
    const int k0 = std::floor(ijk[2]);    const int k1 = k0+1;

    // the decimal part of the (i,j,k) coordinates - becomes the interpolation parameters
    const Real id = ijk[0] - i0;
    const Real jd = ijk[1] - j0;
    const Real kd = ijk[2] - k0;
 
    // check to make sure we are inside the grid boundaries
    if (i0 >= 0 && i0 < _N-1 && j0 >= 0 && j0 < _N-1 && k0 >= 0 && k0 < _N-1)
    {
        return _interpolateDistanceGrid(i0, j0, k0, i1, j1, k1, id, jd, kd);
    }

    // we are outside the bounds of the grid - so clamp to the grid border
    const int clamped_i0 = std::clamp(i0, 0, _N-1);     const int clamped_i1 = std::clamp(i1, 0, _N-1);
    const int clamped_j0 = std::clamp(j0, 0, _N-1);     const int clamped_j1 = std::clamp(j1, 0, _N-1);
    const int clamped_k0 = std::clamp(k0, 0, _N-1);     const int clamped_k1 = std::clamp(k1, 0, _N-1);

    // get the distance and gradient from the border point
    const Real dist_from_border = _interpolateDistanceGrid(clamped_i0, clamped_j0, clamped_k0, clamped_i1, clamped_j1, clamped_k1, id, jd, kd);
    const Vec3r border_point = _gridPointFromIJK(   (clamped_i0 == i0) ? ijk[0] : (Real)clamped_i0,
                                                    (clamped_j0 == j0) ? ijk[1] : (Real)clamped_j0,
                                                    (clamped_k0 == k0) ? ijk[2] : (Real)clamped_k0  );

    // if we don't have the SDF gradient, we can't get the closest point on the mesh
    // so the best we can do is an estimate ==> (Distance From p to SDF Grid border) + (Distance from SDF Grid border to mesh)
    if (!_with_gradient)
    {
        return dist_from_border + (p - border_point).norm();    // THIS IS NOT EXACT! But without the gradient (I think) this is the best we can do
    }

    // if we have the gradient, we can be more exact by getting the closest point on the mesh
    const Vec3r grad_from_border = _interpolateGradientGrid(clamped_i0, clamped_j0, clamped_k0, clamped_i1, clamped_j1, clamped_k1, id, jd, kd);
    
    // find the closest point on the mesh from the grid border
    const Vec3r closest_point = border_point - dist_from_border * grad_from_border;

    // now that we have the closest point on the mesh, we can find the distance
    return (p - closest_point).norm();
}

Vec3r MeshSDF::gradient(const Vec3r& p) const
{
    assert(_with_gradient);
    const Vec3r& ijk = _gridIJKFromPoint(p);
    const int i0 = std::floor(ijk[0]);    const int i1 = i0+1;
    const int j0 = std::floor(ijk[1]);    const int j1 = j0+1;
    const int k0 = std::floor(ijk[2]);    const int k1 = k0+1;

    // the decimal part of the (i,j,k) coordinates - becomes the interpolation parameters
    const Real id = ijk[0] - i0;
    const Real jd = ijk[1] - j0;
    const Real kd = ijk[2] - k0;

    // check to make sure we are inside the grid boundaries
    if (i0 >= 0 && i0 < _N-1 && j0 >= 0 && j0 < _N-1 && k0 >= 0 && k0 < _N-1)
    {
        return _interpolateGradientGrid(i0, j0, k0, i1, j1, k1, id, jd, kd);
    }

    // we are outside the bounds of the grid - so clamp to the grid border
    const int clamped_i0 = std::clamp(i0, 0, _N-1);     const int clamped_i1 = std::clamp(i1, 0, _N-1);
    const int clamped_j0 = std::clamp(j0, 0, _N-1);     const int clamped_j1 = std::clamp(j1, 0, _N-1);
    const int clamped_k0 = std::clamp(k0, 0, _N-1);     const int clamped_k1 = std::clamp(k1, 0, _N-1);

    // get the distance and gradient from the border point
    const Real dist_from_border = _interpolateDistanceGrid(clamped_i0, clamped_j0, clamped_k0, clamped_i1, clamped_j1, clamped_k1, id, jd, kd);
    const Vec3r grad_from_border = _interpolateGradientGrid(clamped_i0, clamped_j0, clamped_k0, clamped_i1, clamped_j1, clamped_k1, id, jd, kd);
    const Vec3r border_point = _gridPointFromIJK(   (clamped_i0 == i0) ? ijk[0] : (Real)clamped_i0,
                                                    (clamped_j0 == j0) ? ijk[1] : (Real)clamped_j0,
                                                    (clamped_k0 == k0) ? ijk[2] : (Real)clamped_k0  );
    // find the closest point on the mesh from the grid border
    const Vec3r closest_point = border_point - dist_from_border * grad_from_border;

    return (p - closest_point).normalized();
}

void MeshSDF::writeToFile(const std::string& filename) const
{
    std::cout << "Writing SDF to " << filename << "..." << std::endl;

    std::stringstream ss;
    ss << _N << " " << _N << " " << _N << "\n"; // print grid size
    ss << formatFloat(_cell_size[0],FLOAT_PRECISION) << " " << formatFloat(_cell_size[1],FLOAT_PRECISION) << " " << formatFloat(_cell_size[2],FLOAT_PRECISION) << "\n";    // print cell size
    ss << formatFloat(_grid_bbox_min[0],FLOAT_PRECISION) << " " << formatFloat(_grid_bbox_min[1],FLOAT_PRECISION) << " " << formatFloat(_grid_bbox_min[2],FLOAT_PRECISION) << "\n";
    ss << formatFloat(_grid_bbox_max[0],FLOAT_PRECISION) << " " << formatFloat(_grid_bbox_max[1],FLOAT_PRECISION) << " " << formatFloat(_grid_bbox_max[2],FLOAT_PRECISION) << "\n";
    ss << formatFloat(_mesh_bbox_min[0],FLOAT_PRECISION) << " " << formatFloat(_mesh_bbox_min[1],FLOAT_PRECISION) << " " << formatFloat(_mesh_bbox_min[2],FLOAT_PRECISION) << "\n";  
    ss << formatFloat(_mesh_bbox_max[0],FLOAT_PRECISION) << " " << formatFloat(_mesh_bbox_max[1],FLOAT_PRECISION) << " " << formatFloat(_mesh_bbox_max[2],FLOAT_PRECISION) << "\n"; 
    ss << formatFloat(_mesh_mass_center[0],FLOAT_PRECISION) << " " << formatFloat(_mesh_mass_center[1],FLOAT_PRECISION) << " " << formatFloat(_mesh_mass_center[2],FLOAT_PRECISION) << "\n";

    // print distances
    for (int i = 0; i < _N; i++)    for (int j = 0; j < _N; j++)    for (int k = 0; k < _N; k++)
    {
        ss << formatFloat(_distance_grid.at(i,j,k),FLOAT_PRECISION) << "\n";
    }

    // print gradients (if applicable)
    if (_with_gradient)
    {
        for (int i = 0; i < _N; i++)    for (int j = 0; j < _N; j++)    for (int k = 0; k < _N; k++)
        {
            Vec3r grad = _gradient_grid.at(i,j,k);
            ss << formatFloat(grad[0],FLOAT_PRECISION) << " " << formatFloat(grad[1],FLOAT_PRECISION) << " " << formatFloat(grad[2],FLOAT_PRECISION) << "\n";
        }
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error opening file " << filename << "!" << std::endl;
        assert(0);
    }

    outfile << ss.str();
    outfile.close();

    std::cout << "Done!" << std::endl;

}

void MeshSDF::_loadSDFFromFile(const std::string& filename)
{
    const int REAL_CHAR_WIDTH = FLOAT_PRECISION + 6;   // the number of characters a Real takes up in the file

    // fast read from file using mmap - https://stackoverflow.com/a/17925197
    struct stat sb;
    long cntr = 0;
    int fd, line_length;
    int line_num = 0;
    char* data;
    char* line;

    int grid_index_1D = 0;
    Real f;
    bool reading_distances = true;
    bool reading_gradients = false;
    // map the file
    fd = open(filename.c_str(), O_RDONLY);
    fstat(fd, &sb);

    data = (char*)mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    line = data;

    // read the first few lines

    _N = atoi(line);    // read the number of grid cells (grid has the same number of cells in each direction)
    _distance_grid.resize(_N, _N, _N);  // resize the distance grid accordingly
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line
    
    parseVector3f(_cell_size, data, REAL_CHAR_WIDTH);  // read cell size
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line

    
    parseVector3f(_grid_bbox_min, data, REAL_CHAR_WIDTH);   // read minimum grid bounding box point
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line
    parseVector3f(_grid_bbox_max, data, REAL_CHAR_WIDTH);   // read maximum grid bounding box point
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line

    parseVector3f(_mesh_bbox_min, data, REAL_CHAR_WIDTH);   // read minimum mesh bounding box point
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line
    parseVector3f(_mesh_bbox_max, data, REAL_CHAR_WIDTH);   // read maximum mesh bounding box point
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line

    parseVector3f(_mesh_mass_center, data, REAL_CHAR_WIDTH);   // read mesh mass center
    while (*data++ != '\n' && cntr++ < sb.st_size);     // move to next line

    // read distances and gradients (when applicable)
    while (cntr < sb.st_size)
    {
        line = data;

        // move to the next line
        while (cntr++ < sb.st_size && *data++ != '\n');
        

        /** Process the line */

        // find 3D grid index (i.e. i,j,k)
        const int i = grid_index_1D / (_N*_N);
        const int j = (grid_index_1D / _N) % _N;
        const int k = grid_index_1D % _N;

        // if we are reading distances, we only expect a single Real per line
        if (reading_distances)
        {
            // read Real from buffer    
            _distance_grid(i,j,k) = atof(line);
        }

        // if we are reading gradients, we expect 3 Reals per line separated by spaces
        if (reading_gradients)
        {
            // read Real 3-vector from char buffer
            parseVector3f(_gradient_grid(i,j,k), line, REAL_CHAR_WIDTH);
        }

        grid_index_1D++;

        // when we reach the end of distances, switch over to reading gradients
        if (grid_index_1D >= _N*_N*_N)
        {
            if (reading_distances)
            {
                reading_distances = false;

                // check if there are gradients after the distances
                if (sb.st_size - cntr > _N*_N*_N*(REAL_CHAR_WIDTH + 3))
                {
                    // if so, resize the gradient grid appropriately
                    _gradient_grid.resize(_N, _N, _N);
                    _with_gradient = true;

                    // and indicate that we are expecting to read gradients now
                    reading_gradients = true;
                }
                
            }
            else if (reading_gradients)
            {
                reading_gradients = false;
            }

            grid_index_1D = 0;
            continue;
        }
    }
    
    std::cout << "Finished reading SDF from file..." << std::endl;
    std::cout << "\tGrid size: " << _N << " x " << _N << " x " << _N << std::endl;
    std::cout << "\tCell size: " << _cell_size[0] << ", " << _cell_size[1] << ", " << _cell_size[2] << std::endl;
    std::cout << "\tGrid Bounding box: (" << _grid_bbox_min[0] << ", " << _grid_bbox_min[1] << ", " << _grid_bbox_min[2] << ") to (" << _grid_bbox_max[0] << ", " << _grid_bbox_max[1] << ", " << _grid_bbox_max[2] << ")" << std::endl;
    std::cout << "\tMesh Bounding box: (" << _mesh_bbox_min[0] << ", " << _mesh_bbox_min[1] << ", " << _mesh_bbox_min[2] << ") to (" << _mesh_bbox_max[0] << ", " << _mesh_bbox_max[1] << ", " << _mesh_bbox_max[2] << ")" << std::endl;
    std::cout << "\tWith gradients: " << (_with_gradient ? "True" : "False") << std::endl;
}

void MeshSDF::_makeSDF(const VertexMat& verts, const TriangleMat& tris, int padding, bool with_gradient)
{
    // set up the grids based on the input parameters
    _distance_grid.resize(_N, _N, _N);

    // compute the bounding box around the vertices
    const Vec3r vertex_mins = verts.rowwise().minCoeff();
    const Vec3r vertex_maxs = verts.rowwise().maxCoeff();

    // compute the cell size for the grid in absolute units, accomodating for padding cells on all sides
    _cell_size = (vertex_maxs - vertex_mins) / (_N - padding*2);

    // compute the SDF bounding box limits
    _grid_bbox_min = vertex_mins - _cell_size * padding;
    _grid_bbox_max = vertex_maxs + _cell_size * padding;

    // store the mesh bounding box limits
    _mesh_bbox_min = vertex_mins;
    _mesh_bbox_max = vertex_maxs;

    _mesh_mass_center = massCenter(verts, tris);

    // initialize distances with really large value
    _distance_grid.assign(std::numeric_limits<Real>::max());


    Array3i closest_tri(_N, _N, _N, -1);    // keeps track of index of closest triangle to each grid point
    Array3i intersection_count(_N, _N, _N, 0); // intersection_count(i,j,k) is # of tri intersections in (i-1,i]x{j}x{k}


    // tracks the number of edge/vertex intersections occur while doing intersection testing
    // by tracking when we intersect with an edge or vertex, we can be sure that we only count an intersection once,
    // leading to accurate inside/outside information
    //
    // the key is the grid index of the grid point (mapped from 3D to 1D)
    // the value is a vector that stores the indices of vertices (or the first vertex in the edge for an edge intersection) directly intersected
    std::map<int, std::vector<unsigned>> direct_hits;

    // no idea what this does, but it was in the original implementation
    const int exact_band = 1;

    for(unsigned ti = 0; ti < tris.cols(); ti++){
        // extract triangle vertices
        const int p = tris(0,ti);
        const int q = tris(1,ti);
        const int r = tris(2,ti);

        const Vec3r& vp = verts.col(p);
        const Vec3r& vq = verts.col(q);
        const Vec3r& vr = verts.col(r);

        // coordinates in grid to high precision
        Real fip=((Real)vp[0]-_grid_bbox_min[0])/_cell_size[0], fjp=((Real)vp[1]-_grid_bbox_min[1])/_cell_size[1], fkp=((Real)vp[2]-_grid_bbox_min[2])/_cell_size[2];
        Real fiq=((Real)vq[0]-_grid_bbox_min[0])/_cell_size[0], fjq=((Real)vq[1]-_grid_bbox_min[1])/_cell_size[1], fkq=((Real)vq[2]-_grid_bbox_min[2])/_cell_size[2];
        Real fir=((Real)vr[0]-_grid_bbox_min[0])/_cell_size[0], fjr=((Real)vr[1]-_grid_bbox_min[1])/_cell_size[1], fkr=((Real)vr[2]-_grid_bbox_min[2])/_cell_size[2];

        // do distances nearby
        int i0=std::clamp(int(std::min({fip,fiq,fir}))-exact_band, 0, _N-1), i1=std::clamp(int(std::max({fip,fiq,fir}))+exact_band+1, 0, _N-1);
        int j0=std::clamp(int(std::min({fjp,fjq,fjr}))-exact_band, 0, _N-1), j1=std::clamp(int(std::max({fjp,fjq,fjr}))+exact_band+1, 0, _N-1);
        int k0=std::clamp(int(std::min({fkp,fkq,fkr}))-exact_band, 0, _N-1), k1=std::clamp(int(std::max({fkp,fkq,fkr}))+exact_band+1, 0, _N-1);

        for (int k = k0; k <= k1; k++)  for (int j = j0; j <= j1; j++)  for (int i = i0; i <= i1; i++)
        {
            const Vec3r grid_point(i*_cell_size[0]+_grid_bbox_min[0], j*_cell_size[1]+_grid_bbox_min[1], k*_cell_size[2]+_grid_bbox_min[2]);
            Real dist = pointTriangleDistance(grid_point, vp, vq, vr);
            if(dist < _distance_grid(i,j,k)){
                _distance_grid(i,j,k) = dist;
                closest_tri(i,j,k) = ti;
            }
        }
        // and do intersection counts
        j0=std::clamp((int)std::ceil(std::min({fjp,fjq,fjr})), 0, _N-1);
        j1=std::clamp((int)std::floor(std::max({fjp,fjq,fjr})), 0, _N-1);
        k0=std::clamp((int)std::ceil(std::min({fkp,fkq,fkr})), 0, _N-1);
        k1=std::clamp((int)std::floor(std::max({fkp,fkq,fkr})), 0, _N-1);
        for (int k = k0; k <= k1; k++)  for (int j = j0; j <= j1; j++)
        {
            // check if this grid point is inside the triangle projected in the YZ-plane
            Real a, b, c;
            bool in_triangle = pointInTriangle2D(j, k, fjp, fkp, fjq, fkq, fjr, fkr, a, b, c);
            
            if (in_triangle)
            {
                Real fi=a*fip+b*fiq+c*fir; // intersection i coordinate
                int i_interval=std::max(0, int(std::ceil(fi))); // intersection is in (i_interval-1,i_interval]

                // SPECIAL CASE: the ray directly intersects a vertex or edge of the triangle
                //  - since vertices and edges are shared by multiple triangles, we must keep track of the vertices/edges we hit directly so that
                //    we only count 1 intersection
                if ( (a-MESH2SDF_EPS) <= 0 || (b-MESH2SDF_EPS) <= 0 || (c-MESH2SDF_EPS) <= 0)
                {
                    unsigned direct_hit = -1;
                    if ( (a-MESH2SDF_EPS) <= 0 && (b-MESH2SDF_EPS) <= 0)      direct_hit = r; // intersects vertex R
                    else if ( (a-MESH2SDF_EPS) <= 0 && (c-MESH2SDF_EPS) <= 0) direct_hit = q; // intersects vertex Q
                    else if ( (b-MESH2SDF_EPS) <= 0 && (c-MESH2SDF_EPS) <= 0) direct_hit = p; // intersects vertex P
                    else if ( (a-MESH2SDF_EPS) <= 0) direct_hit = std::min(q,r); // intersects edge QR
                    else if ( (b-MESH2SDF_EPS) <= 0) direct_hit = std::min(p,r); // intersects edge PR
                    else if ( (c-MESH2SDF_EPS) <= 0) direct_hit = std::min(p,q); // intersects edge PQ
                    else    assert(0);  // shouldn't ever get to here

                    // convert (i,j,k) to a unique integer key for the map
                    int key = i_interval*_N*_N + j*_N + k;
                    
                    // if key does not exist in map, create an empty vector there
                    if (direct_hits.count(key) == 0)
                    {
                        direct_hits[key] = std::vector<unsigned>();
                    }

                    std::vector<unsigned>& vec = direct_hits[key];
                    // if we've already counted an intersection for this edge or vertex, skip the intersection counting
                    if (std::find(vec.begin(), vec.end(), direct_hit) != vec.end())
                        continue;
                    // otherwise register that we've hit this edge or vertex, and go on to count the intersection
                    else
                        vec.push_back(direct_hit);
                }

                // count the intersection
                if(i_interval < 0) ++intersection_count(0, j, k); // we enlarge the first interval to include everything to the -x direction
                else if(i_interval < _N) ++intersection_count(i_interval,j,k);
                // we ignore intersections that are beyond the +x side of the grid
            }
        }
    }

    // and now we fill in the rest of the distances with fast sweeping
    for(unsigned int pass=0; pass<2; ++pass){
        _sweep(verts, tris, closest_tri, +1, +1, +1);
        _sweep(verts, tris, closest_tri, -1, -1, -1);
        _sweep(verts, tris, closest_tri, +1, +1, -1);
        _sweep(verts, tris, closest_tri, -1, -1, +1);
        _sweep(verts, tris, closest_tri, +1, -1, +1);
        _sweep(verts, tris, closest_tri, -1, +1, -1);
        _sweep(verts, tris, closest_tri, +1, -1, -1);
        _sweep(verts, tris, closest_tri, -1, +1, +1);
    }

    // then figure out signs (inside/outside) from intersection counts
    for (int k=0; k<_N; ++k)    for (int j=0; j<_N; ++j)
    {
        int total_count=0;
        for(int i=0; i<_N; ++i){
            total_count+=intersection_count(i,j,k);
            if(total_count%2==1){ // if parity of intersections so far is odd,
                _distance_grid(i,j,k) = -_distance_grid(i,j,k); // we are inside the mesh
            }
        }
    }

    // evaluate gradients if required
    if (with_gradient)
    {
        _gradient_grid.resize(_N, _N, _N);
        // do gradients from closest tris
        for (int k=0; k<_N; ++k) for (int j=0; j<_N; ++j) for (int i=0; i<_N; ++i)
        {
            const int ti = closest_tri(i,j,k);
            // extract triangle vertices
            const int p = tris(0,ti);
            const int q = tris(1,ti);
            const int r = tris(2,ti);

            const Vec3r& vp = verts.col(p);
            const Vec3r& vq = verts.col(q);
            const Vec3r& vr = verts.col(r);

            const Vec3r grid_point = _gridPointFromIJK(i,j,k);
            const Vec3r grad = pointTriangleDirection(grid_point, vp, vq, vr);

            if (_distance_grid.at(i,j,k) < 0)   _gradient_grid.at(i,j,k) = -grad;
            else                                _gradient_grid.at(i,j,k) = grad;
            
        }
    }

    std::cout << "Finished creating SDF from mesh..." << std::endl;
    std::cout << "\tGrid size: " << _N << " x " << _N << " x " << _N << std::endl;
    std::cout << "\tCell size: " << _cell_size[0] << ", " << _cell_size[1] << ", " << _cell_size[2] << std::endl;
    std::cout << "\tBounding box: (" << _grid_bbox_min[0] << ", " << _grid_bbox_min[1] << ", " << _grid_bbox_min[2] << ") to (" << _grid_bbox_max[0] << ", " << _grid_bbox_max[1] << ", " << _grid_bbox_max[2] << ")" << std::endl;
    std::cout << "\tWith gradients: " << (_with_gradient ? "True" : "False") << std::endl;
    
}

void MeshSDF::_checkNeighbor(const VertexMat& verts, const TriangleMat& tris, 
                        Array3i& closest_tri,
                        const Vec3r& grid_point,
                        int i0, int j0, int k0, int i1, int j1, int k1)
{
    const int ti = closest_tri(i1,j1,k1);
    if (ti >= 0)
    {
        // extract triangle vertices
        const int p = tris(0,ti);
        const int q = tris(1,ti);
        const int r = tris(2,ti);

        const Vec3r& vp = verts.col(p);
        const Vec3r& vq = verts.col(q);
        const Vec3r& vr = verts.col(r);
        
        Real dist = pointTriangleDistance(grid_point, vp, vq, vr);
        if (dist < _distance_grid(i0,j0,k0))
        {
            _distance_grid(i0,j0,k0) = dist;
            closest_tri(i0,j0,k0) = ti;
        }
    }
}

void MeshSDF::_sweep(const VertexMat& verts, const TriangleMat& tris,
                Array3i& closest_tri,
                int di, int dj, int dk)
{
    int i0, i1;
    if(di>0){ i0=1; i1=_N; }
    else{ i0=_N-2; i1=-1; }
    int j0, j1;
    if(dj>0){ j0=1; j1=_N; }
    else{ j0=_N-2; j1=-1; }
    int k0, k1;
    if(dk>0){ k0=1; k1=_N; }
    else{ k0=_N-2; k1=-1; }

    for(int k=k0; k!=k1; k+=dk) for(int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di){
      const Vec3r grid_point = _gridPointFromIJK(i,j,k);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j,    k   );
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i,    j-dj, k   );
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j-dj, k   );
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i,    j,    k-dk);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j,    k-dk);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i,    j-dj, k-dk);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j-dj, k-dk);
   }

}

Vec3r MeshSDF::_gridPointFromIJK(int i, int j, int k) const
{
    return Vec3r(i*_cell_size[0]+_grid_bbox_min[0], j*_cell_size[1]+_grid_bbox_min[1], k*_cell_size[2]+_grid_bbox_min[2]);
}

Vec3r MeshSDF::_gridPointFromIJK(Real i, Real j, Real k) const
{
    return Vec3r(i*_cell_size[0]+_grid_bbox_min[0], j*_cell_size[1]+_grid_bbox_min[1], k*_cell_size[2]+_grid_bbox_min[2]);
}

Vec3r MeshSDF::_gridIJKFromPoint(const Vec3r& p) const
{
    return (p - _grid_bbox_min).array() / _cell_size.array();
}

Real MeshSDF::_interpolateDistanceGrid(int i0, int j0, int k0, int i1, int j1, int k1, int id, int jd, int kd) const
{
    // trilinear interpolation
    const Real c000 = _distance_grid.at(i0,j0,k0);    const Real c100 = _distance_grid.at(i1,j0,k0);
    const Real c001 = _distance_grid.at(i0,j0,k1);    const Real c101 = _distance_grid.at(i1,j0,k1);
    const Real c010 = _distance_grid.at(i0,j1,k0);    const Real c110 = _distance_grid.at(i1,j1,k0);
    const Real c011 = _distance_grid.at(i0,j1,k1);    const Real c111 = _distance_grid.at(i1,j1,k1);
    const Real c00 = c000*(1-id) + c100*id;
    const Real c01 = c001*(1-id) + c101*id;
    const Real c10 = c010*(1-id) + c110*id;
    const Real c11 = c011*(1-id) + c111*id;
    const Real c0 = c00*(1-jd) + c10*jd;
    const Real c1 = c01*(1-jd) + c11*jd;
    const Real c = c0*(1-kd) + c1*kd;

    return c;
}

Vec3r MeshSDF::_interpolateGradientGrid(int i0, int j0, int k0, int i1, int j1, int k1, int id, int jd, int kd) const
{
    const Vec3r& c000 = _gradient_grid.at(i0,j0,k0);    const Vec3r& c100 = _gradient_grid.at(i1,j0,k0);
    const Vec3r& c001 = _gradient_grid.at(i0,j0,k1);    const Vec3r& c101 = _gradient_grid.at(i1,j0,k1);
    const Vec3r& c010 = _gradient_grid.at(i0,j1,k0);    const Vec3r& c110 = _gradient_grid.at(i1,j1,k0);
    const Vec3r& c011 = _gradient_grid.at(i0,j1,k1);    const Vec3r& c111 = _gradient_grid.at(i1,j1,k1);
    return vectorTriSlerp(c000, c100, c010, c110, c001, c101, c011, c111, id, jd, kd);
}

} // mesh2sdf