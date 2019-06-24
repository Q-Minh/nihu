Mesh building {#tut_mesh_building}
=============

\page tut_mesh_building

[Eigen]:http://eigen.tuxfamily.org/index.php?title=Main_Page
[mexFunction]:http://www.mathworks.com/help/matlab/apiref/mexfunction.html
[MEX]:http://www.mathworks.com/help/matlab/create-mex-files.html
[Matlab frontend]:\ref install_matlab_interface
[OFF]:http://people.sc.fsu.edu/~jburkardt/data/off/off.html
[gmsh]:http://geuz.org/gmsh/

[TOC]

Introduction {#tut_mesh_intro}
============

The C++ NiHu core does not contain versatile mesh building methods for creating parametric
meshes or to import different mesh formats. NiHu provides a single mesh building
interface, the NiHu::create_mesh function that reads two mesh description matrices:
one for the vertex locations and one for the element nodal connectivities.

\note For versatile mesh generation and mesh importing functionalities, refer to the [Matlab frontend].

The create_mesh function {#tut_mesh_create_mesh}
========================

The first example {#tut_mesh_first}
-----------------

The NiHu::create_mesh function takes two matrices as its first arguments:
- `nodes` is the nodal location matrix
- `elements` is the nodal connectivity matrix.

As NiHu uses the linear algebra library [Eigen] for matrix manipulations, we use two `typedef`-s for [Eigen] matrices.

\snippet mesh_building.mex.cpp Typedefs

- `dMatrix` denotes a dynamically resizeable matrix with real (`double`) coefficients,
- `uMatrix` denotes the same with `unsigned` coefficients.

We instantiate the two matrices by allocating their coefficients, and initialising them using [Eigen]'s comma initialiser syntax:

\snippet mesh_building.mex.cpp Matrices

- The rows of matrix `nodes` contain 3D locations.
- The rows or matrix `elements` contain an element id followed by as many element nodal indices as needed by the element type. Unused coefficients are left zero.

Apparently, we are building a heterogeneous mesh consisting of linear quadrangles (NiHu::quad_1_elem) and triangles (NiHu::tria_1_elem).
The mesh is finally built by calling the NiHu::create_mesh function

\snippet mesh_building.mex.cpp Creation

The additional arguments NiHu::quad_1_tag () and NiHu::tria_1_tag () are additional type information telling the compiler which element types are used in the mesh.
From this information, the compiler builds a specific mesh type that stores the specified element types separately, in an optimised way.

\note The function's return type is hidden from the library user by the `auto` keyword. In fact, it is
~~~~~~~~~~~~
mesh<tmp::vector<quad_1_elem, tria_1_elem> >
~~~~~~~~~~~~
but this information is not needed to use the library.

Building a 2D mesh {#tut_mesh_line}
------------------

The next example builds a circle from line elements in 2D:

\snippet mesh_building.mex.cpp 2D Circle

Working with other matrix formats {#tut_mesh_othermatrix}
---------------------------------

The NiHu::create_mesh function is a template that can be called with any indexable matrix types.
The only requirement against the matrix types is that they should be indexable with integers, using a `(row,col)`-type syntax.

This feature, together with the NiHu::mex::real_matrix class, makes easy Matlab [MEX] integration possible.

\snippet mesh_building.mex.cpp Matlab example

\note For more information on the Matlab MEX interface and on the arguments of the [mexFunction], refer to the Matlab [MEX] site.

Reading a mesh from an `OFF` file {#tut_mesh_offfile}
=================================

`NiHu` can import meshes stored in [OFF].
This can easily be done using the NiHu::read_off_mesh function.
The usage of this function is demonstrated below.

~~~~~~~~~~~~
auto msh = read_off_mesh("a_mesh_file.off", tria_1_tag(), quad_1_tag())
~~~~~~~~~~~~

The element types to be read from the file are listed by passing the NiHu::read_off_mesh function instances of the corresponding tags.

\note Only NiHu::tria_1_tag and quad_1_tag are currently supported by the NiHu::read_off_mesh function.

