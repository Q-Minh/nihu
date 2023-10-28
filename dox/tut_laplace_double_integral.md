Double integral of the Laplace kernel {#tut_laplace_double_integral}
=====================================

\page tut_laplace_double_integral Double integral of the Laplace kernel

[Eigen]:http://eigen.tuxfamily.org/index.php?title=Main_Page

[TOC]

Introduction {#tut_laplace_double_intro}
============

This tutorial explains how to use `NiHu` to compute a singular double integral of the form

\f$
\displaystyle
I = \int_{S} \int_{S} K({\bf x},{\bf y}) dS_y dS_x,
\f$

where \f$ S \f$ is a square surface

\f$ S = \left\{{\bf x} = (\xi,\eta,0) : -1 \le \xi \le 1, -1 \le \eta \le 1 \right\} \f$

and the integrand kernel \f$ K \f$ is the fundamental solution of the Laplace equation in 3D

\f$ K({\bf x},{\bf y}) = 1 / 4\pi r, \quad r = |{\bf x}-{\bf y}| \f$.

The integrand contains an \f$ O(1/r) \f$ singularity when \f$ {\bf x} = {\bf y} \f$.
The analytical result is

\f$
\displaystyle
I = \frac{32}{4\pi} \left[ \log \left( 1+\sqrt{2} \right) - \frac{\sqrt{2}-1}{3} \right]
\f$.


The C++ code {#tut_laplace_double_Cpp_code}
============

Libraries and typedefs {#tut_laplace_double_includes}
----------------------

We need to include two modules:
- core/weighted_residual.hpp is the main module of the NiHu C++ core
- library/laplace_kernel.hpp defines our kernel.

\snippet laplace_double_integral.cpp Includes

We are going to work with dynamically resizeable Eigen matrices. We define two convenient typedefs for double and unsigned matrices.

\snippet laplace_double_integral.cpp Typedefs


The main function {#tut_laplace_double_main_function}
-----------------

Our problem domain \f$ S \f$ is represented by a mesh consisting of a single square element. For more information on mesh building, refer to the tutorial \ref tut_mesh_building.

\snippet laplace_double_integral.cpp Mesh

We reformulate our original problem in the form

\f$
\displaystyle
I = \sum_{i,j} I_{ij}, \quad 
\f$
where
\f$ 
\quad
I_{ij} = \left< w_i, \left(\mathcal{K} w_j\right)_S \right>_S,
\f$

where the set \f$ w_i \f$ denote a piecewise constant or isoparametric function space, and \f$ \mathcal{K} \f$ denotes the integral operator based on the Laplace kernel.

Two function spaces are generated from the mesh, in order to test both field generation options:

\snippet laplace_double_integral.cpp Operators


The testing function {#tut_laplace_double_test_function}
--------------------

As the different function space views are of different types, the tester function is written as a template that can receive any function space type as parameter.

\snippet laplace_double_integral.cpp Test

The function first allocates the result matrix. The matrix size is determined by the total number of DOF in the funcion space `w`, computed by the member function NiHu::function_space_base::get_num_dofs.

After instantiating the integral operator and evaluating the weighted double integral into the matrix, the tester function prints the matrix, its element sum, and compares the sum to the analytical integral.

The complete source of the tutorial is found here: tutorial/laplace_double_integral.cpp

