Single integral of the Laplace kernel {#tut_laplace_single_integral}
=====================================

\page tut_laplace_single_integral

[TOC]

Introduction {#tut_laplace_intro}
============

This tutorial explains how to use NiHu to compute a single integral of the form

\f$ I = \int_{S} K({\bf x}_0,{\bf y}) dS_y \f$

where \f$ S \f$ is a square surface

\f$ S = \left\{{\bf x} = (\xi,\eta,0) : 0 \le \xi \le 1, 0 \le \eta \le 1 \right\} \f$,

\f$ {\bf x}_0 \f$ is a point located at the surafce's center \f$ (.5, .5, 0) \f$,
and the integrand kernel \f$ K \f$ is the fundamental solution of the Laplace equation in 3D

\f$ K({\bf x},{\bf y}) = 1 / 4\pi r, \quad r = |{\bf x}-{\bf y}| \f$

The analytical result of the singular integral is

\f$ I = \log \left(1+\sqrt{2}\right) / \pi \f$


The C++ code {#tut_laplace_single_Cpp_code}
============

Mesh generation {#tut_laplace_single_mesh_generation}
---------------

The mesh of our surface domain will consist of only one 4-noded quadrilateral (::quad_1_elem) element:

\snippet laplace_single_integral.cpp Mesh

The library function ::create_mesh creates a mesh that contains only ::quad_1_elem elements.


Function space and weighted integral {#tut_laplace_single_function_space}
------------------------------------

The single integral is generalised as a double weighted integral where one of the weighting functions is a Dirac delta function and the other is constant:

\f$ I = \int_{S} K({\bf x}_0,{\bf y}) dS_y =
\int_{S} \delta({\bf x}_0) \int_{S} K({\bf x},{\bf y}) 1({\bf x}) dS_y dS_y \f$

\snippet laplace_single_integral.cpp Function spaces

As in the previous tutorial, we convert our mesh into a function space by creating a constant function space view using the library function ::constant_view.
As a second step, we create a Dirac-version of our function space using the library function ::dirac.
This function space consists of Dirac delta functions located at the nodal coordinates of the argument function space.
In our special case, this is the element center.

Note that the function spaces are stored in references.
This shows that the created objects are light-weight _views_ of the arguments (the mesh object).

As our mesh contains one element, and the piecewise constant function space contains one DOF per element, the resulting
weighted integral matrix will contain one entry.

\snippet laplace_single_integral.cpp Weighted residual

Finally, we instantiated the integral operator and performed the weighted integral using the abstract operator notation.


Results {#tut_laplace_single_results}
-------

The integration result is printed and compared to the analytical solution.

\snippet laplace_single_integral.cpp Results


Summary {#tut_laplace_single_summary}
=======

This introductory tutorial explained how to use NiHu to evaluate singular single integrals.
It was shown
- how a simple mesh is created from scratch
- how light weight function space views are created from meshes
- that single integrals are weighted double integrals with Dirac delta weighting functions
- how a function space is converted into a light weight Dirac-view
- how integral operators and weighted residuals are defined and evaluated using an abstract syntax

The complete source of the tutorial is found here: tutorial/laplace_single_integral.cpp

