Parallel operator evaluations {#tut_couple_radiation}
=============================

\page tut_couple_radiation

[TOC]

Introduction {#tut_couple_radiation_intro}
============

In most BEM applications two or four integral operators need to be discretised to form the system of equations.
For the case of an acoustic BEM, the single and double layer potential kernels of the Helmholtz equation need to be discretised simultaneously to obtain the sytem matrices.
If the Burton-Miller method (see e.g. \ref tut_helmholtz_bem_3d_fict) is applied, then four integral operators need to be discretised using the same test and trial function spaces.

In most cases, the parallel evaluation of the different kernels is much cheaper than sequential evaluations, because
- the kernel inputs (locations, normal vectors) need to be computed only once
- the kernels can be computed from each other: for example, if \f$ G \f$ denotes the 3D acoustic Green's function, then \f$ G' = -G/r \cdot (1+\mathrm{i}kr) \f$
- the weighting shape functions and the Jacobian need to be computed only once.

Couple kernels, operators and matrices  {#tut_couple_radiation_couples}
======================================

The expression _couple_, which is extensively used in the following (and in the nomenclature of the source code) refers to a compile-time collection of classes derived from the same base classes, using static polymorphism (see \ref tech_poly).
Without the need of going into details, it is enough for our understanding that certain types can be collected into couples and the couple container also provides some functionality applyable on the whole container.
For example indexing a couple of matrices of the same type will result in a couple of the indexed elements, whereas the couple of a number of kernel functions can be evaluated resulting in a couple object consisting of the results of each kernel evaluated with the same arguments.
In the following couple objects will be traited as transparent layers over the contained data objects.

NiHu provides a transparent mechanism to automatically optimise parallel kernel evaluations.

- Kernels can be encapsulated into couple kernels.
- Couple kernels can be used to define couple integral operators
- Couple integral operators can be discretised into couple matrices, just like ordinary operators.

Without going into details with respect to the implementation, we introduce an acoustic application.

Application  {#tut_couple_radiation_application}
===========

We compute the pressure and velocity fields radiated from a closed surface \f$ S \f$.
We assume that both the surface pressure and velocity fields are known on the boundary suraface \f$ S \f$, and we are interested in the radiated pressure and velocioty fields on a field point surface \f$ F \f$.

The pressure field is described as

\f$
\displaystyle
p({\bf x}_i) = \left(\mathcal{L}q\right)_S({\bf x}_i) - \left(\mathcal{M}p\right)_S({\bf x}_i), \quad {\bf x}_i \in F.
\f$

The velocity field is computed by means of differentiation:

\f$
\displaystyle
q({\bf x}_i) = \left(\mathcal{M}^{\mathrm{T}}q\right)_S({\bf x}_i) - \left(\mathcal{N}p\right)_S({\bf x}_i) \quad {\bf x_i} \in F.
\f$

Program structure  {#tut_couple_radiation_structure}
=================

The C++ application is going to compute the discretised operators two times.
1. First with the sequential evaluation
2. And then with couple kernel evaluations


The C++ code  {#tut_couple_radiation_cpp}
============

Function spaces  {#tut_couple_radiation_funcspace}
---------------

The C++ code receives the surface and field point meshes from Matlab.
The trial function space is an isoparametric view of the mesh, the test function space is a dirac-view.
The implementation is straightforward:

\snippet couple_radiation.mex.cpp Straightforward

Sequential evaluations  {#tut_couple_radiation_sequential}
----------------------

The sequential evaluation of the matrices is done as usual:

\snippet couple_radiation.mex.cpp Without coupling

The C++ application measures the time needed to compute the matrices.
The time result is passed to Matlab as a 1x1 real matrix.

Parallel evaluations  {#tut_couple_radiation_parallel}
--------------------

The parallel evaluation of the matrices is done by creating a couple integral operator by passing multiple kernel functions as arguments to the library function ::create_integral_operator.
The result matrices are encapsulated into a couple using the library function ::create_couple

\snippet couple_radiation.mex.cpp With coupling

The Matlab tester  {#tut_couple_radiation_matlab}
================

The results can be compared using a Matlab tester like shown below:

\include couple_radiation_test.m

Possible printed results are

	L log10error:	-4.003
	M log10error:	-4.569
	Mt log10error:	-4.292
	N log10error:	-Inf
	time gain:	0.4496

Reults and discussion  {#tut_couple_radiation_results}
=====================

The results are interpreted as follows:
- The couple kernel evaluation accelerates the radiation integrals by a factor of approximately two.
- The N matrices are identical, the other matrices are slightly different.
This is because parallel kernel evaluations means that the four kernels need to be evaluated in the same integration points.
The number of integration points is defined by the (worst case) hypersingular kernel, so the other integrals are computed with unnecessarily large accuracy.

The full codes of this tutorial are available here:
- C++ code: tutorial/couple_radiation.mex.cpp
- Matlab code: tutorial/couple_radiation_test.m
