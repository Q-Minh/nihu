Mesh building {#tut_mesh_building}
=============

\page tut_mesh_building

[Eigen]:http://eigen.tuxfamily.org/index.php?title=Main_Page
[mexFunction]:http://www.mathworks.com/help/matlab/apiref/mexfunction.html
[MEX]:http://www.mathworks.com/help/matlab/create-mex-files.html

[TOC]

Introduction {#tut_mesh_intro}
============

The C++ NiHu core does not contain versatile mesh building methods for creating parametric
meshes or to import different mesh formats. NiHu provides a single mesh building
interface, the ::create_mesh function that reads two mesh description matrices,
`nodes` and `elements`, where `nodes` contains the vertex locations and
`elements` contains the connectivity information.

For versatile mesh generation and mesh importing functionalities, refer to the Matlab frontend.

The create_mesh function {#tut_mesh_create_mesh}
========================

The first example {#tut_mesh_first}
-----------------

The ::create_mesh function takes two matrices as its first arguments:
- `nodes` is the nodal location matrix
- `elements` is the nodal connectivity matrix.

As NiHu uses the linear algebra library [Eigen] for matrix manipulations, we use two `typedef`-s for [Eigen] matrices.

\snippet mesh_building.mex.cpp Typedefs

- `dMatrix` denotes a dynamically resizable matrix with real (`double`) coefficients,
- `uMatrix` denotes the same with `unsigned` coefficients.

We instantiate the two matrices by allocating their coefficients, and initialising them using [Eigen]'s comma initialiser syntax:

\snippet mesh_building.mex.cpp Matrices

- The rows of matrix `nodes` contains 3D locations.
- The rows or matrix `elements` contain an element id followed by element nodal indices.

Apparently, we are building a heterogeneous mesh consisting of linear quadrangles (::quad_1_elem) and triangles (::tria_1_elem).

The mesh is finally built by calling the ::create_mesh function

\snippet mesh_building.mex.cpp Creation

The additional arguments ::_quad_1_tag () and ::_tria_1_tag () are used to tell the C++ compiler the element types stored in the mesh.
From this information, the compiler creates its specific mesh type that handles the two element types separately, in optimised way.
The function's return type is hidden from the library user by the `auto` keyword. In fact, it is

	mesh<tmp::vector<quad_1_elem, tria_1_elem> >

but this information is not needed to use the library.


Building a 2D mesh {#tut_mesh_line}
------------------

The next example builds a circle from line elements in 2D:

\snippet mesh_building.mex.cpp 2D Circle


Working with other matrix formats {#tut_mesh_othermatrix}
---------------------------------

The ::create_mesh function is a template that can be called with any indexable matrix types.
The only requirement against the matrix types is that they should be indexable with `(row,col)`-type indexing, where `row` and `col` are integers.

This feature, together with the ::mex::real_matrix class, makes easy Matlab [MEX] integration possible.

\snippet mesh_building.mex.cpp Matlab example

For more information on the Matlab [MEX] interface and on the arguments of the [mexFunction], refer to the Matlab [MEX] site.

