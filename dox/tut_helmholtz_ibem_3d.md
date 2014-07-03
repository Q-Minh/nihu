An indirect collocational Helmholtz BEM example  {#tut_helmholtz_ibem_3d}
===============================================

\page tut_helmholtz_ibem_3d

[TOC]

Introduction {#tut_helmholtz_ibem_3d_intro}
============

This tutorial demonstrates the usage of NiHu for solving an exterior Dirichlet problem in 3D by means of an indirect collocational type boundary element method.

Theory {#tut_helmholtz_ibem_3d_theory}
======

The indirect BEM, applied to an exterior Dirichlet problem, expresses the radiated pressure from a vibrating surface by a double layer potential:

\f$ \displaystyle
p({\bf x}) = \left(\mathcal{M}\sigma\right)_S({\bf x}) \quad {\bf x} \in F
\f$

where \f$ \sigma({\bf y}) \f$ denotes an appropriate surface source density function without direct physical meaning.

If the point \f$ \bf x \f$ approaches the boundary surface \f$ S \f$, a boundary integral equation is obtained:

\f$
\displaystyle
p({\bf x}) = \left(\left[\mathcal{M} + \frac{1}{2}\mathcal{I}\right]\sigma\right)_S({\bf x}) \quad {\bf x} \in S
\f$

The boundary integral equation can be used to obtain the source density function on the surface, once a prescribed surface pressure is defined.
In a second step, the radiated pressure is computed by evaluating the double layer potential integrals for external field points.


Program structure {#tut_helmholtz_ibem_3d_structure}
=================

We are going to implement a Matlab-C++ NiHu application, where
- the C++ executable is responsible for assembling the system matrices from the boundary surface mesh and the field point mesh
- the Matlab part defines the meshes, calls the C++ executable, solves the systems of equations, computes the radiated pressure and quantifies the error of the solution.

The C++ code {#tut_helmholtz_ibem_3d_cpp}
============

The C++ code is going to be called from Matlab as

	[Ms, Mf] = helmholtz_ibem_3d(
	    surf_nodes, surf_elements, field_nodes, field_elements, wave_number);
	
\note The computed system matrix `Ms` will contain the discretised identity operator too.

Without going into details regarding the C++ code, we examplify the lines where the operators are discretised using collocation:

\snippet helmholtz_ibem_3d.mex.cpp System


The Matlab code {#tut_helmholtz_ibem_3d_matlab}
=============== 

The example creates a sphere surface of radius \f$ R = 1\,\mathrm{m} \f$ centred at the origin, consisting of triangular elements only.
The field point mesh is a surface obtained by revolving a line around the radiator.

The Dirichlet boundary conditions are defined by a point source located at the point \f$ \mathbf{x}_0 \f$.

The compiled mex file can be called from Matlab, as demonstrated in the example below:

\include helmholtz_ibem_3d_test.m

The generated plot looks like

\image html tut_helmholtz_ibem_3d.png

And the resulting errors read as

	log10 mean error on field = -2.42

The full codes of this tutorial are available here:
- C++ code: tutorial/helmholtz_ibem_3d.mex.cpp
- Matlab code: tutorial/helmholtz_ibem_3d_test.m

