# Mesh2SDF
Takes an input .obj or .stl mesh and computes its signed distance function (SDF).

The core algorithm is the same as [Christopher Batty's SDFGen](https://github.com/christopherbatty/SDFGen), with some augmentations/fixes:
- The gradient of the SDF is now calculated, which is stored in addition to the signed distance (though can be disabled by using `with_gradients=0`).
- The SDF is now queryable with global coordinates (instead of SDF grid indices). Trilinear interpolation is used when the queried point does not fall exactly on a grid point.
- The SDF inside/outside testing was made slightly more robust by handling the edge case where the tested intersection ray directly intersects a mesh edge or vertex.
- The algorithms were adapted to use [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), for cross-compatibility with other projects of mine that use this library.
- All extraneous code from SDFGen was removed


## Building
[`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) is used under the hood for representing vectors and matrices. This can be installed for CMake easily on Linux with:
```
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build
cmake ..
make
sudo make install
```

[`assimp`](https://github.com/assimp/assimp) is an optional requirement that is used to load mesh files. On Ubuntu, it can be installed simply with
```
sudo apt-get install libassimp-dev
```
If the `assimp` library is not installed on your system, only simple `.obj` files (with only vertices and triangular faces, no vertex normals) can be parsed and loaded.

To build, perform the following steps from the home directory of the repo:
```
mkdir build
cd build
cmake ..
make
```
By default, the library is built to perform computations using floating point arithmetic. To use double precision, pass the command line argument `-DUSE_DOUBLE_PRECISION=True` with the CMake command, i.e.
```
cmake .. -DUSE_DOUBLE_PRECISION=True
```

## Executables
Building will produce two binaries that let you generate and query grid-based SDFs of triangle meshes.
The first binary, `Mesh2SDF`, generates a SDF from a mesh file of your choice. For example:
```
./Mesh2SDF   <filename>   <num_cells>   <padding>   <with_gradients>
./Mesh2SDF   obj/cube.obj 128           8           1
```
will create a size 128<sup>3</sup> SDF that has 8 cells of padding between the object bounding box and the boundary of the SDF grid. It will also compute the gradient of the SDF at each grid point. Finally, it will write this SDF (and its gradients) to the file `output.sdf`.

The second binary, `QuerySDF`, models useful API features such as creating a SDF from a mesh, saving it to file, loading it back from file, and querying distances and gradients at various points.

## API
The core functionality is contained in the class `MeshSDF`, demonstrated in the code snippet below.
```C++
// === Writing to File ===
const auto [verts, tris] = loadMeshDataFromFile('cube.obj');  // load vertices and triangles from mesh file
MeshSDF sdf(verts, tris, num_cells, padding, with_gradients); // create SDF from mesh data
sdf.writeToFile("output.sdf");

// === Loading from File ===
MeshSDF sdf("output.sdf");

// === Querying the SDF ===
Eigen::Vector3f p(0.2, 0.1, 0.4);
float dist = sdf.evaluate(p); // query the SDF at some arbitrary point - will do trilinear interpolation
Eigen::Vector3f grad = sdf.gradient(p) // get the gradient of SDF at some arbitrary - will do trilinear interpolation
```
