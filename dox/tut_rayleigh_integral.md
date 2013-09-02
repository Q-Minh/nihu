Rayleigh integral {#tut_rayleigh_integral}
=================

\page tut_rayleigh_integral

[TOC]

Introduction {#tut_rayleigh_intro}
============

This tutorial explains how to use NiHu to evaluate the Rayleigh integral.
The Rayleigh integral describes acoustic radiation from an infinite vibrating plane.

In this tutorial, we learn how to connect NiHu with Matlab.
The problem geometry and parameters will be imported from, and the results will be passed back to Matlab using NiHu's Matlab interface.

Theory {#tut_rayleigh_theory}
======

The Rayleigh integral {#tut_rayleigh_definition}
---------------------

The Rayleigh integral computes the acoustic pressure \f$ p({\bf x}) \f$ radiated from a vibrating planar surface \f$ S \f$ embedded into an infinite rigid plane.
The vibrating surface is characterised by its normal velocity function \f$ v({\bf y}) \f$.
The radiated pressure is evaluated as

\f$
\displaystyle
p({\bf x})
= \int_S G({\bf x}, {\bf y}) v({\bf y}) dS_y
= \left(\mathcal{G}v\right)_S({\bf x})
\f$

where \f$ G \f$ denotes the fundamental solution of the Helmholtz equation in 3D:

\f$ G({\bf x}, {\bf y}) = \exp(-ikr)/4\pi r, \quad r = |{\bf y} - {\bf x}|\f$

Discretisation {#tut_rayleigh_discretisation}
--------------

The velocity field on the surface \f$ S \f$ is discretised using piecewise constant or isoparametric interpolation functions:

\f$ v({\bf y}) = \sum_j w_j({\bf y}) v_j \f$

The radiated pressure (the Rayleigh integral) is evaluated in a set of field points \f$ {\bf x}_i \f$.
For later convenience, we consider the field points as locations in a so called field point domain \f$ F \f$.
Evaluating the pressure in the field points is equivalent to pre-multiplying the Rayleigh integral with Dirac delta functions located at the field points \f$ {\bf x}_i \f$ and integrating with respect to the variable \f$ {\bf x} \f$ over the field point domain:

\f$
\displaystyle
p({\bf x}_i)
= \sum_j \left(\mathcal{G}w_j \right)_S({\bf x}_i) \cdot v_j
= \sum_j \left< \delta_{{\bf x}_i}, \left(\mathcal{G}w_j \right)_S \right>_F \cdot v_j
= \sum_j Z_{ij} v_j
\f$
 
where \f$ Z_{ij} \f$ denotes the transfer impedance matrix.


Program structure {#tut_rayleigh_structure}
-----------------

Our program is organised as follows:
- The two meshes (the radiating surface and the field point mesh) will be generated in Matlab.
- The mesh description matrices (nodes and elements for both meshes) and the wave number \f$ k \f$ will be passed to the C++ program from Matlab.
- The C++ program generates the function spaces and evaluates the transfer impedance matrix. The program also computes the radiating surface's radiation impedance matrix that relates the normal velocity to the sound pressure on the surface.
- The transfer and radiation impedance matrices are passed back to Matlab


The C++ code {#tut_rayleigh_Cpp_code}
============

Included modules and typedefs {#tut_rayleigh_include}
-----------------------------

We need three modules from the NiHu library
- core/weighted_residual.hpp is the main module of NiHu
- library/helmholtz_kernel.hpp defines the fundamental solution of the Helmholtz equation
- util/mex_matrix.hpp includes the Matlab interface functions

We further define two `typedef`s for Matlab matrices ::mex::real_matrix and ::mex::complex_matrix.
These matrices are allocated in Matlab's memory, in Matlab format, but are used by NiHu just as if they were C++ Eigen matrices.

\snippet rayleigh_integral_3d.mex.cpp Header

The MEX function
----------------

Our executable code will be called from Matlab using the syntax

	[Z_trans, Z_rad] = rayleigh_integral_3d(surf_nodes, surf_elements, field_nodes, field_elements, wave_number);
	
The input and output arguments are catched by the C++ function `mexFunction`

\snippet rayleigh_integral_3d.mex.cpp Mex function

whose parameters are
- `nlhs` the number of left hand side (output) parameters (two in our case)
- `lhs` array of pointers pointing to the Matlab memory, where the output arguments (`Z_trans` and `Z_rad`) are allocated
- `nrhs` the number of right hand side input parameters (five in our case)
- `rhs` array of pointers pointing to the Matlab memory, where the five input arguments are allocated

Mesh generation {#tut_rayleigh_mesh_generation}
---------------

The radiator surface mesh and the field point mesh is instantiated in C++ as follows:

\snippet rayleigh_integral_3d.mex.cpp Meshes

- The Matlab mesh description matrices are simply imported into C++ by the library class ::mex::real_matrix (`dMatrix`).
The class' constructor refers to the Matlab input pointer.
The resulting `dMatrix` objects are light-weight interfaces (_views_) providing convenient indexing capabilities to the Matlab data.

Function space generation {#tut_rayleigh_function_space}
-------------------------

The function spaces are generated from the two meshes as follows:

\snippet rayleigh_integral_3d.mex.cpp Function spaces

- We use an isoparametric function space view of the radiator surface mesh to discretise the velocity field. This means that the velocity nodes are located in the elements' corners, and the velocity is interpolated using linear functions within an element.
- The transfer impedance is evaluated in the element centers of the field point mesh. For this reason, we create a Dirac-view of a piecewise constant function space generated from the field point mesh.
- The radiation impedance is evaluated in the element centers of the radiating surface mesh. For this reason, we create a Dirac-view of a piecewise constant function space generated from the radiating surface mesh.

Kernel and weighted residual {#tut_rayleigh_kernel}
----------------------------

We instantiate the integral operator \f$ \mathcal{G} \f$ using the Helmholtz kernel as follows:

\snippet rayleigh_integral_3d.mex.cpp Kernel

- The real wave number is imported from Matlab with the standard Matlab interface function `mxGetPr` that returns a pointer to the real data.
- The Helmholtz kernel is instantiated using the class ::helmholtz_3d_SLP_kernel templated to the wave number's type (`double`).
The abbreviation SLP refers to the single layer potential.

We allocate memory for the output complex matrices:

\snippet rayleigh_integral_3d.mex.cpp Matrices

- The three-argument constructor (rows, columns, output pointer) of class ::mex::complex_matrix (`cMatrix`) allocates a Matlab output matrix of given dimensions in Matlab memory.

The weighted double integrals are evaluated using the operator notations:

\snippet rayleigh_integral_3d.mex.cpp Weighted residual

The two last lines of code are syntactically identical, but there is a great difference.
- The transfer impedance matrix is defined between two separate meshes, and is therefore computed by evaluating regular integrals only
- The radiation impedance matrix is defined on one mesh, so its evaluation involves singular integration.

The complete source of the C++ code is found here: tutorial/rayleigh_integral_3d.mex.cpp

The Matlab code {#tut_rayleigh_Matlab_code}
===============

The compiled executable can be called from Matlab, as demonstrated in the example below:

\include rayleigh_integral_3d_test.m

The generated plot looks like

\image html tut_rayleigh_integral.png

