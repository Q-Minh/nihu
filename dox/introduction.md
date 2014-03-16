Introduction {#introduction}
============

\page introduction

[TOC]

What is NiHu? {#introduction_whatis}
=============

NiHu is an open source C++-Matlab toolbox used for solving boundary value problems of partial differential equations by means of the boundary element method (BEM).
Specifically, the toolbox provides a general framework to numerically evaluate weighted residual integrals of the form

\f$
\displaystyle
W_{ij} = \left<t_i, \left(\mathcal{K}d_j\right)_S\right>_F
=
\int_F t_i({\bf x})
\int_S K({\bf x}, {\bf y}) d_j({\bf y}) \mathrm{d} S_{{\bf y}}
\mathrm{d} F_{{\bf x}}
\f$

where
\f$ S \f$ denotes a boundary surface,
\f$ t_i \f$ and \f$ d_j \f$ denote base functions of discretised function spaces, and
\f$ \mathcal{K} \f$ denotes an integral operator defined by its kernel function \f$ K \f$.
Such integrals frequently arise when the numerical solution of boundary value problems is formalised using the finite or boundary element methods (see for example \ref theo_bem_example).

Unified programming interface {#introduction_unified}
=============================

The main motivation of developing NiHu is to provide a unified open source software framework for BEM problems of different kinds.
Our implementation is capable to generate BEM executables
- with Direct and Indirect BEM formalism
- with Galerkin, Collocational and General solution methods
- in 2D and 3D
- with element types defined in a general manner

The main C++ core of NiHu can be considered as a skeleton that defines different BEM problems in a general way, by reflecting mathematics behind boundary elements in the C++ code.
The C++ code snippet below demonstrates the capabilities of NiHu's application layer:
~~~~~~~~~
auto msh = create_mesh(nodes, elements, quad_1_tag(), tria_1_tag());
auto func_sp = constant_view(msh);
auto L = create_integral_operator(laplace_3d_SLP_kernel());
auto M = create_integral_operator(laplace_3d_DLP_kernel());
auto I = identity_integral_operator();
L_matrix << dirac(func_sp) * L[func_sp];
M_matrix << dirac(func_sp) * M[func_sp] + dirac(func_sp) * (-.5*I)[func_sp];
~~~~~~~~~
The 7-line code snippet defines a collocational 3D BEM for the Laplace equation, implemented with constant triangular and quadrilateral elements, and computes the system matrices of the discretised boundary integrals.

Customisable library {#introduction_library}
====================

Specific aspects of BEM implementations are developed in C++ libraries designed to be easily customisable by researchers and developers.
Developers may customise
- Element types
- Kernels
- Quadratures
- Singular integration techniques specific to kernels, formalisms and element types
- and more...

in order to adapt NiHu to their specific problems.

Matlab layer {#introduction_matlab}
============

The Matlab shell layer provides extensive pre- and post processing features like
- Parametric mesh generation for academic benchmarks
- Standard mesh format import and export
- Direct and iterative solutions of the resulting matrix equations
- Easy plot generation

and makes integration with other software tools straightforward.

Summary {#introduction_summary}
=======

The introduced concept allows a beneficial application of the toolbox in
- **education**, where general mathematical aspects need to be emphasized,
- **research**, where the developer can focus on small scale details, while having his implementation working as a part of a general system;
- and **engineering practice**, as the implementation is capable of solving problems of industrial size.
