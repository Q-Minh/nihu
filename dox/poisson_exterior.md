An exterior Poisson problem {#poisson_extrior_tutorial}
===========================


[Matlab MEX]:http://www.mathworks.com/help/matlab/matlab_external/c-c-source-mex-files.html

[TOC]

Introduction {#intro}
============

This tutorial explains how to use NiHu to generate the BEM system matrices used to solve an exterior poisson problem governed by the equation

\f$\nabla^2 p = 0\f$

The problem mesh will be generated in Matlab, the system matrices are going to be built in the C++ core, and the system is going to be solved in Matlab,

Mesh generation {#mesh_generation}
===============

Matrix generation {#matrix_generation}
=================

Now we are going to build the C++ code that imports the Matlab-borne mesh,
generates the system matrices and passes them back to Matlab.

Libraries and typedefs {#includes}
----------------------

We need to include three modules:
- util/mex_matrix.hpp is included to import and export Matlab MEX matrices
- bem/weighted_residual.hpp is the main module of the NiHu C++ core
- library/poisson_kernel.hpp defines the kernel specific to our Poisson problem.

As we are going to work with double scalars, we define a convenient typedef for the corresponding Matlab MEX matrix `dMatrix`

\snippet potential_bem.mex.cpp includes

The C++ MEX header {#header}
------------------

The C++ starting point is the conventional `mexFunction` with parameters defining the number of left and right hand side Matlab arguments, as well as their pointer arrays. For the sake of brevity, we simply check the appropriate number of input and output arguments now. More sophisticated argument check methods can be found on the [Matlab MEX] site.

\snippet potential_bem.mex.cpp header

Integral operators {#integral_operators}
------------------

For the solution of our problem, we need three integral operators:
- \f$\mathcal{L}\f$ is the single layer potential operator based on the Green's function of our Poisson problem.
- \f$\mathcal{M}\f$ is the double layer potential operator based on the Green's function's normal derivative.
- \f$\mathcal{B}\f$ is the identity operator scaled by a factor of \f$-1/2\f$.

These integral operators are simply imported from the library module library/poisson_kernel.hpp. Operators \f$\mathcal{L}\f$ and \f$\mathcal{M}\f$ are instantiated using instances of their kernels. Operator \f$\mathcal{B}\f$ is default constructed.

\snippet potential_bem.mex.cpp integral operators

The boundary mesh {#surface_mesh}
-----------------

We now build the problem's boundary mesh \f$S\f$ from our mesh description matrices passed as the first and second right hand side Matlab arguments. The Matlab-borne MEX matrices are imported into C++ by constructing `dMatrix` objects from their Matlab matrix pointers.

The mesh is simply constructed from the Matlab-borne matrices. The third (and possibly subsequent) argument(s) of the factory function `create_mesh` define the element types the mesh consists of. This a-priori information is important for the compiler to optimize the C++ code for each different element type of our problem.

After the mesh is built, we generate a piecewise constant function space above the problem's geometry, the so-called `constant_view` of our boundary.

\snippet potential_bem.mex.cpp surface mesh

The boundary matrices {#surface_system}
---------------------

We generate our boundary integral equations based on the Galerkin method: 
The boundary system matrices of our problem are defined as

\f${\bf L}_s = \left<v_i, \mathcal{L}v_j\right>\f$,

\f${\bf M}_s = \left<v_i, \mathcal{M}v_j\right>\f$,

\f${\bf B}_s = \left<v_i, \mathcal{B}v_j\right> = -1/2 \cdot \left<v_i, v_i\right> \delta_{ij}\f$

where \f$v_i\f$ denote the basis functions of our piecewise constant function space.

The output Matlab MEX matrices are allocated by constructing `dMatrix` objects of given dimensions. The matrix diemenisons are determined by the number of degrees of freedoms in our function space, returned by function `get_num_dofs()`. The matlab matrices are immediately allocated in the Matlab memory, by referring to their pointers as the third arguments of the `dMatrix` constructors. The constructors automatically set the matrix elements to zero.

For simplicity, we immediately sum up matrices \f${\bf M}_s\f$ and \f${\bf B}_s\f$ into a single matrix.


\snippet potential_bem.mex.cpp surface system

The field point mesh {#field_mesh}
--------------------

Radiation from the boundary surface to field points is described by the field point matrices defined as

\f${\bf L}_f = \left<\delta_i, \mathcal{L}v_j\right>\f$,

\f${\bf M}_f = \left<\delta_i, \mathcal{M}v_j\right>\f$.

where \f$\delta_i\f$ denote Dirac delta functions at the field point locations.

Although the matrices are simply determined by a set of field points, at least for plotting purposes it is convenient to define a field point surface mesh consisting of elements and nodes. Therefore, the field point mesh, consisting of linear quadrilateral elements, is generated similar to the boundary surface mesh. The constant view again means that the field points are located at the element centers, and the additional dirac wrapper indicates that the field point system matrices will be generated by evaluating single integrals, where the kernels are bound to the nodal locations of the field point function space.

\snippet potential_bem.mex.cpp field point mesh

The boundary matrices {#surface_system}
---------------------

The field point matrices are no longer square. Their dimensions are defined by the number of degrees of freedoms in the field point and surface function spaces.

\snippet potential_bem.mex.cpp field point system

We are ready. The boundary system and field point matrices are exported into Matlab.


System solution {#system_solution}
===============


