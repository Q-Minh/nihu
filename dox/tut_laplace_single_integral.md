Single integral of the Laplace kernel {#tut_laplace_single_integral}
=====================================

\page tut_laplace_single_integral Single integral of the laplace kernel

[TOC]

Introduction {#tut_laplace_single_intro}
============

This tutorial explains how to use NiHu to compute a single singular integral of the form

\f$
\displaystyle
I = \int_{S} K({\bf x}_0,{\bf y}) dS_y,
\f$

where \f$ S \f$ is a square surface

\f$
\displaystyle
S = \left\{{\bf x} = (\xi,\eta,0) : 0 \le \xi \le 1, 0 \le \eta \le 1 \right\}
\f$,

\f$ {\bf x}_0 \f$ is a point located at the surface's center \f$ (.5, .5, 0) \f$,
and the integrand kernel \f$ K \f$ is the fundamental solution of the Laplace equation in 3D

\f$
\displaystyle
K({\bf x},{\bf y}) = 1 / 4\pi r, \quad r = |{\bf x}-{\bf y}| .
\f$

The analytical result of the singular integral is

\f$
\displaystyle
I = \log \left(1+\sqrt{2}\right) / \pi .
\f$


The C++ code {#tut_laplace_single_Cpp_code}
============

Mesh generation {#tut_laplace_single_mesh_generation}
---------------

The mesh of our surface domain will consist of only one 4-noded quadrilateral (NiHu::quad_1_elem) element:

\snippet laplace_single_integral.cpp Mesh


Function space and weighted integral {#tut_laplace_single_function_space}
------------------------------------

The single integral is generalised as a double weighted integral where one of the weighting functions is a Dirac delta function and the other is constant:

\f$
\displaystyle
I = \int_{S} K({\bf x}_0,{\bf y}) dS_y =
\int_{S} \delta({\bf x}_0) \int_{S} K({\bf x},{\bf y}) 1({\bf x}) dS_y dS_x
= \left< \delta_{{\bf x}_0}, \left(\mathcal{K}1\right)_S \right>_S
\f$

The corresponding function spaces are instantiated as

\snippet laplace_single_integral.cpp Function spaces

As our mesh contains one element, and the piecewise constant function space contains one DOF per element, the resulting
weighted integral matrix will contain one entry.

\snippet laplace_single_integral.cpp Weighted residual


Results {#tut_laplace_single_results}
-------

The integration result is printed and compared to the analytical solution.

\snippet laplace_single_integral.cpp Results


The complete source of the tutorial is found here: tutorial/laplace_single_integral.cpp

