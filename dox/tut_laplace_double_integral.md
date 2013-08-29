Double integral of the Laplace kernel {#tut_laplace_double_integral}
=====================================

\page tut_laplace_double_integral

[TOC]

Introduction {#tut_laplace_intro}
============

This tutorial explains how to use NiHu to compute a double integral of the form

\f$ I = \int_{S} \int_{S} K({\bf x},{\bf y}) dS_y dS_x \f$

where \f$ S \f$ is a square surface

\f$ S = \left\{{\bf x} = (\xi,\eta,0) : -1 \le \xi \le 1, -1 \le \eta \le 1 \right\} \f$

and the integrand kernel \f$ K \f$ is the fundamental solution of the Laplace equation in 3D

\f$ K({\bf x},{\bf y}) = 1 / 4\pi r, \quad r = |{\bf x}-{\bf y}| \f$

The analytical value of the integral is

\f$ I = \frac{32}{4\pi} \left[ \log \left( 1+\sqrt{2} \right) - \frac{\sqrt{2}-1}{3} \right] \f$


The C++ code {#tut_laplace_double_Cpp_code}
============

Libraries and typedefs {#tut_laplace_double_includes}
----------------------

We need to include two modules:
- bem/weighted_residual.hpp is the main module of the NiHu C++ core
- library/laplace_kernel.hpp defines our kernel.

\snippet laplace_double_integral.cpp Includes

We are going to work with dynamically resizable Eigen matrices. We define two convenient typedefs for double and unsigned matrices.

\snippet laplace_double_integral.cpp Typedefs


Mesh generation {#tut_laplace_double_mesh_generation}
---------------

Our problem domain \f$ S \f$ is represented by a mesh consisting of linear triangular and quadrilateral elements.

	+1 +---+---+
	   |   |4 /|
	   | 2 | / |
	   |   |/ 3|
	 0 +---+---+
	   |   |   |
	   | 0 | 1 |
	   |   |   |
	-1 +---+---+
	  -1   0  +1

Note that the integral is singular on each element pair of our mesh, and all singularity types
- ::FACE_MATCH
- ::EDGE_MATCH
- ::CORNER_MATCH

are included

We define a function `build_mesh` to build the mesh from scratch.

\snippet laplace_double_integral.cpp Mesh generation

The mesh is built using two mesh description matrices.
- `nodes` contains the nodal coordinates in its rows,
- `elements` contains the element description data in its rows.

Both matrices are initialised using Eigen's convenient comma operator syntax.
The element description data consists of an element type identifier (::quad_1_elem::id or ::tria_1_elem::id) followed by the element's nodal indices.
Unused entries in the matrix are left zero.

The library Function ::create_mesh creates the mesh from these inputs.
In order to generate a mesh optimised to our input data, the compiler needs to know the element types used in the mesh description matrices.
For this reason, the function takes additional type parameters (::_tria_1_tag (), ::_quad_1_tag ()).

The compiler can deduce the function's return type from the return expression, using the `auto - decltype` syntax (C++11 feature).
In the `decltype` expression, the matrices are replaced by their defaults.

The main function {#tut_laplace_double_main_function}
-----------------

As NiHu works in a more general context, we reformulate our original problem by extending our integration domain \f$ S \f$ by a weighting function

\f$ 1(x) = \sum_{i} W_i(x) \f$,

and performing the double integral with the weighting functions.

\f$ I = \int_{S} \sum_i W_i(x) \int_{S} K(x,y) \sum_j W_j(y) dS_y dS_x = \sum_{i} \sum_{j} I_{ij}\f$

where

\f$ I_{ij} = \int_{S} W_i(x) \int_{S} K(x,y) W_j(y) dS_y dS_x \f$

or with operator notation

\f$ I_{ij} = \left< w_i, \mathcal{K} w_j \right> \f$

where the set \f$ w_i \f$ denote a function space, and \f$ \mathcal{K} \f$ denotes an integral operator based on the Laplace kernel.

The function space is generated from the mesh in the main function:

\snippet laplace_double_integral.cpp Main

The function uses two formalisms:
- ::constant_view means that the function space consists of piecewise constant functions
- ::isoparametric_view (in our special case) means that the function space consists of piecewise linear weighting functions (for more general information, see \todo LINK HERE)

The terminology `_view` indicates that these function spaces contain no additional data, only indicate that the mesh is considered to contain constant or isoparametric elements.


The testing function {#tut_laplace_double_test_function}
--------------------

As the different function spaces are of different types, the tester function is written as a function template that can receive any function space type as parameter.

\snippet laplace_double_integral.cpp Test

The function first allocates the result matrix. The matrix size is determined by the total number of DOF in the funcion space `w`, computed by the member function ::function_space_base::get_num_dofs.

As a next step, the integral operator \f$ \mathcal{K} \f$ is instantiated from the Laplace kernel, and the double integral is evaluated.

The tester function prints the matrix, its element sum, and compares the sum to the analytical integral.

Summary {#tut_laplace_double_summary}
=======

This introductory tutorial explained how to use NiHu to evaluate singular double integrals.
It was shown
- how heterogoneous meshes are created from scratch
- how light weight function space views are created from meshes
- how integral operators and weighted residuals are defined and evaluated using an abstract syntax

The complete source of the tutorial is found here: tutorial/laplace_double_integral.cpp
