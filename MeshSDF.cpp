#include "MeshSDF.hpp"
#include "util.hpp"

#include <map>
#include <iostream>

MeshSDF::MeshSDF(const VertexMat& verts, const TriangleMat& tris, int grid_size, int padding, bool with_gradient)
    : _N(grid_size), _with_gradient(with_gradient)
{
    _makeSDF(verts, tris, padding, with_gradient);
}

float MeshSDF::evaluate(const Eigen::Vector3f& p) const
{
    const Eigen::Vector3f& ijk = _gridIJKFromPoint(p);
    // trilinear interpolation
    const int i0 = std::floor(ijk[0]);    const int i1 = i0+1;
    const int j0 = std::floor(ijk[1]);    const int j1 = j0+1;
    const int k0 = std::floor(ijk[2]);    const int k1 = k0+1;
    const double id = ijk[0] - i0;
    const double jd = ijk[1] - j0;
    const double kd = ijk[2] - k0;

    const float c000 = _distance_grid.at(i0,j0,k0);    const float c100 = _distance_grid.at(i1,j0,k0);
    const float c001 = _distance_grid.at(i0,j0,k1);    const float c101 = _distance_grid.at(i1,j0,k1);
    const float c010 = _distance_grid.at(i0,j1,k0);    const float c110 = _distance_grid.at(i1,j1,k0);
    const float c011 = _distance_grid.at(i0,j1,k1);    const float c111 = _distance_grid.at(i1,j1,k1);
    const float c00 = c000*(1-id) + c100*id;
    const float c01 = c001*(1-id) + c101*id;
    const float c10 = c010*(1-id) + c110*id;
    const float c11 = c011*(1-id) + c111*id;
    const float c0 = c00*(1-jd) + c10*jd;
    const float c1 = c01*(1-jd) + c11*jd;
    const float c = c0*(1-kd) + c1*kd;
    return c;
}

Eigen::Vector3f MeshSDF::gradient(const Eigen::Vector3f& p) const
{
    assert(_with_gradient);
    const Eigen::Vector3f& ijk = _gridIJKFromPoint(p);
    const int i0 = std::floor(ijk[0]);    const int i1 = i0+1;
    const int j0 = std::floor(ijk[1]);    const int j1 = j0+1;
    const int k0 = std::floor(ijk[2]);    const int k1 = k0+1;
    const double id = ijk[0] - i0;
    const double jd = ijk[1] - j0;
    const double kd = ijk[2] - k0;

    const Eigen::Vector3f& c000 = _gradient_grid.at(i0,j0,k0);    const Eigen::Vector3f& c100 = _gradient_grid.at(i1,j0,k0);
    const Eigen::Vector3f& c001 = _gradient_grid.at(i0,j0,k1);    const Eigen::Vector3f& c101 = _gradient_grid.at(i1,j0,k1);
    const Eigen::Vector3f& c010 = _gradient_grid.at(i0,j1,k0);    const Eigen::Vector3f& c110 = _gradient_grid.at(i1,j1,k0);
    const Eigen::Vector3f& c011 = _gradient_grid.at(i0,j1,k1);    const Eigen::Vector3f& c111 = _gradient_grid.at(i1,j1,k1);
    return vectorTriSlerp(c000, c100, c010, c110, c001, c101, c011, c111, id, jd, kd);
}

void MeshSDF::_makeSDF(const VertexMat& verts, const TriangleMat& tris, int padding, bool with_gradient)
{
    // set up the grids based on the input parameters
    _distance_grid.resize(_N, _N, _N);

    // compute the bounding box around the vertices
    const Eigen::Vector3f vertex_mins = verts.rowwise().minCoeff();
    const Eigen::Vector3f vertex_maxs = verts.rowwise().maxCoeff();

    // compute the cell size for the grid in absolute units, accomodating for padding cells on all sides
    _cell_size = (vertex_maxs - vertex_mins) / (_N - padding*2);

    // compute the SDF bounding box limits
    _bbox_min = vertex_mins - _cell_size * padding;
    _bbox_max = vertex_maxs + _cell_size * padding;

    // initialize distances with really large value
    _distance_grid.assign(std::numeric_limits<float>::max());


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

        const Eigen::Vector3f& vp = verts.col(p);
        const Eigen::Vector3f& vq = verts.col(q);
        const Eigen::Vector3f& vr = verts.col(r);

        // coordinates in grid to high precision
        double fip=((double)vp[0]-_bbox_min[0])/_cell_size[0], fjp=((double)vp[1]-_bbox_min[1])/_cell_size[1], fkp=((double)vp[2]-_bbox_min[2])/_cell_size[2];
        double fiq=((double)vq[0]-_bbox_min[0])/_cell_size[0], fjq=((double)vq[1]-_bbox_min[1])/_cell_size[1], fkq=((double)vq[2]-_bbox_min[2])/_cell_size[2];
        double fir=((double)vr[0]-_bbox_min[0])/_cell_size[0], fjr=((double)vr[1]-_bbox_min[1])/_cell_size[1], fkr=((double)vr[2]-_bbox_min[2])/_cell_size[2];

        // do distances nearby
        int i0=std::clamp(int(std::min({fip,fiq,fir}))-exact_band, 0, _N-1), i1=std::clamp(int(std::max({fip,fiq,fir}))+exact_band+1, 0, _N-1);
        int j0=std::clamp(int(std::min({fjp,fjq,fjr}))-exact_band, 0, _N-1), j1=std::clamp(int(std::max({fjp,fjq,fjr}))+exact_band+1, 0, _N-1);
        int k0=std::clamp(int(std::min({fkp,fkq,fkr}))-exact_band, 0, _N-1), k1=std::clamp(int(std::max({fkp,fkq,fkr}))+exact_band+1, 0, _N-1);

        for (int k = k0; k <= k1; k++)  for (int j = j0; j <= j1; j++)  for (int i = i0; i <= i1; i++)
        {
            const Eigen::Vector3f grid_point(i*_cell_size[0]+_bbox_min[0], j*_cell_size[1]+_bbox_min[1], k*_cell_size[2]+_bbox_min[2]);
            float dist = pointTriangleDistance(grid_point, vp, vq, vr);
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
            double a, b, c;
            bool in_triangle = pointInTriangle2D(j, k, fjp, fkp, fjq, fkq, fjr, fkr, a, b, c);
            
            if (in_triangle)
            {
                double fi=a*fip+b*fiq+c*fir; // intersection i coordinate
                int i_interval=std::max(0, int(std::ceil(fi))); // intersection is in (i_interval-1,i_interval]

                // SPECIAL CASE: the ray directly intersects a vertex or edge of the triangle
                //  - since vertices and edges are shared by multiple triangles, we must keep track of the vertices/edges we hit directly so that
                //    we only count 1 intersection
                if ( (a-EPS) <= 0 || (b-EPS) <= 0 || (c-EPS) <= 0)
                {
                    unsigned direct_hit = -1;
                    if ( (a-EPS) <= 0 && (b-EPS) <= 0)      direct_hit = r; // intersects vertex R
                    else if ( (a-EPS) <= 0 && (c-EPS) <= 0) direct_hit = q; // intersects vertex Q
                    else if ( (b-EPS) <= 0 && (c-EPS) <= 0) direct_hit = p; // intersects vertex P
                    else if ( (a-EPS) <= 0) direct_hit = std::min(q,r); // intersects edge QR
                    else if ( (b-EPS) <= 0) direct_hit = std::min(p,r); // intersects edge PR
                    else if ( (c-EPS) <= 0) direct_hit = std::min(p,q); // intersects edge PQ
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

            const Eigen::Vector3f& vp = verts.col(p);
            const Eigen::Vector3f& vq = verts.col(q);
            const Eigen::Vector3f& vr = verts.col(r);

            const Eigen::Vector3f grid_point = _gridPointFromIJK(i,j,k);
            const Eigen::Vector3f grad = pointTriangleDirection(grid_point, vp, vq, vr);

            if (_distance_grid.at(i,j,k) < 0)   _gradient_grid.at(i,j,k) = -grad;
            else                                _gradient_grid.at(i,j,k) = grad;
            
        }
    }
    
}

void MeshSDF::_checkNeighbor(const VertexMat& verts, const TriangleMat& tris, 
                        Array3i& closest_tri,
                        const Eigen::Vector3f& grid_point,
                        int i0, int j0, int k0, int i1, int j1, int k1)
{
    const int ti = closest_tri(i1,j1,k1);
    if (ti >= 0)
    {
        // extract triangle vertices
        const int p = tris(0,ti);
        const int q = tris(1,ti);
        const int r = tris(2,ti);

        const Eigen::Vector3f& vp = verts.col(p);
        const Eigen::Vector3f& vq = verts.col(q);
        const Eigen::Vector3f& vr = verts.col(r);
        
        float dist = pointTriangleDistance(grid_point, vp, vq, vr);
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
      const Eigen::Vector3f grid_point = _gridPointFromIJK(i,j,k);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j,    k   );
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i,    j-dj, k   );
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j-dj, k   );
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i,    j,    k-dk);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j,    k-dk);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i,    j-dj, k-dk);
      _checkNeighbor(verts, tris, closest_tri, grid_point, i, j, k, i-di, j-dj, k-dk);
   }

}

Eigen::Vector3f MeshSDF::_gridPointFromIJK(int i, int j, int k) const
{
    return Eigen::Vector3f(i*_cell_size[0]+_bbox_min[0], j*_cell_size[1]+_bbox_min[1], k*_cell_size[2]+_bbox_min[2]);
}

Eigen::Vector3f MeshSDF::_gridIJKFromPoint(const Eigen::Vector3f& p) const
{
    return (p - _bbox_min).array() / _cell_size.array();
}