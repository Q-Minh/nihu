Scattering of a plane wave from an infinite cylinder {#tut_helmholtz_bem_2d_cyl}
====================================================

\page tut_helmholtz_bem_2d_cyl

[TOC]

Introduction {#tut_helmholtz_bem_2d_cyl_intro}
============

This tutorial demonstrates the usage of NiHu for solving an exterior acoustic scattering problem in 2D by means of a collocational type boundary element method (BEM).
The modeled problem is scattering of an incident plane wave from a rigid cylinder of radius \f$R_0\f$, centered at the origin of the 2D coordinate system.

Theory {#tut_helmholtz_bem_2d_cyl_theory}
======

The incident wave field is given as

\f$
\displaystyle
p_{\mathrm{inc}}({\bf x}) = \exp(-\mathrm{i}k {\bf x} \cdot {\bf d})
\f$

where \f$ {\bf d} \f$ denotes the unit direction vector of propagation.
The scattering is considered as superposing a reflected field to the incident pressure field so that the sum of the two fields results a zero velocity on the scatterer surface:

\f$
\displaystyle
p_{\mathrm{tot}}({\bf x}) = p_{\mathrm{inc}}({\bf x}) + p_{\mathrm{ref}}({\bf x}) \\
p'_{\mathrm{tot}}({\bf x}) = 0 \quad {\bf x} \in S
\f$

This consideration results in the Neumann boundary condition to the reflected wave field:

\f$
\displaystyle
p'_{\mathrm{ref}}({\bf x}) = -p'_{\mathrm{inc}}({\bf x}) \quad {\bf x} \in S
\f$

The problem is easily formulated with conventional BEM operator notations:

\f$
\displaystyle
q_{\mathrm{ref}} = -p'_{\mathrm{inc}}, \quad {\bf x} \in S, \\
\frac{1}{2} p_{\mathrm{ref}}({\bf x}) = \left(\mathcal{M}p_{\mathrm{ref}}\right)_S({\bf x}) - \left(\mathcal{L}q_{\mathrm{ref}}\right)_S({\bf x}), \quad {\bf x} \in S, \\
p_{\mathrm{ref}}({\bf x}) = \left(\mathcal{M}p_{\mathrm{ref}}\right)_S({\bf x}) - \left(\mathcal{L}q_{\mathrm{ref}}\right)_S({\bf x}), \quad {\bf x} \in F
\f$

and can be discretised using either collocation or Galerkin formalisms.

Program structure {#tut_helmholtz_bem_2d_cyl_structure}
=================

We are going to implement a Matlab-C++ NiHu application, where
- the C++ executable is responsible for assembling the system matrices from the boundary surface mesh and the field point mesh,
- the Matlab part defines the meshes, calls the C++ executable, solves the systems of equations, and compares the solutions by quantifying their errors.

The C++ code {#tut_helmholtz_bem_2d_cyl_cpp}
============

The C++ code is going to be called from Matlab as

	[Ls, Ms, Lf, Mf] = helmholtz_bem_2d_cyl(
	    surf_nodes, surf_elements, field_nodes, field_elements, wave_number);
	
\note The computed surface system matrix `Ms` will contain the discretised identity operator too, as required by the conventional BEM formalisms.

Header and mesh creation {#tut_helmholtz_bem_2d_cyl_header}
------------------------

The header of our C++ source file includes the necessary header files and defines some basic types for convenience.
We further include library/helmholtz_singular_integrals.hpp that defines the specialised singular integrals of the Helmholtz kernel for constant lines.

\snippet helmholtz_bem_2d_cyl.mex.cpp Header

Two meshes will be created, both meshes contain line elements only:
- One for the radiator surface (`surf_mesh`)
- One for the field points (`field_mesh`)

\snippet helmholtz_bem_2d_cyl.mex.cpp Mesh

Definition of function spaces {#tut_helmholtz_bem_2d_cyl_spaces}
-----------------------------

Next, we create the discretised function spaces of our mesh.
The definitions are straightforward, as only constant line boundary elements are dealt with.
For the definition of the field point function space, we immediately apply the dirac view.

\snippet helmholtz_bem_2d_cyl.mex.cpp Function spaces

The number of degrees of freedom on the radiating cylinder \f$ S \f$ and in the field mesh \f$ F \f$ are denoted by \f$ n \f$ and \f$ m \f$, respectively.
After the function spaces are built, the number of DOFs are known, and the memory for the system matrices can be preallocated.
The surface system matrices \f$ \mathbf{L}_s \f$, \f$ \mathbf{M}_s \f$ are of size \f$ n \times n \f$, whereas the matrices describing the field point radiation \f$ \mathbf{L}_f \f$ and \f$ \mathbf{M}_f \f$ are of size \f$ m \times n \f$.
The four system matrices are preallocated by means of the following lines of code.

\snippet helmholtz_bem_2d_cyl.mex.cpp Matrices

Integral operators and their evaluation {#tut_helmholtz_bem_2d_cyl_intop}
---------------------------------------

In the next steps, the three integral operators \f$ \mathcal{L} \f$, \f$ \mathcal{M} \f$, and \f$ \mathcal{I} \f$ are defined.
Since the kernel functions depend on the wave number \f$ k \f$, which is passed as the fifth right hand side parameter of the mex function, the integral operators are created as follows.

\snippet helmholtz_bem_2d_cyl.mex.cpp Integral operators

Finally, the system matrices are obtained by the evaluation of the integral operators on the function spaces defined above.
Note that as a collocational formalism is applied, the Dirac-view of the test function spaces is taken for the boundary integrals.

\snippet helmholtz_bem_2d_cyl.mex.cpp System matrices

The Matlab code {#tut_helmholtz_bem_2d_cyl_matlab}
=============== 

The example creates a sphere surface of radius \f$ R = 1\,\mathrm{m} \f$ centred at the origin, consisting of triangular elements only.
The CHIEF points are located on the surface of a cube, shifted from the origin.

The Neumann boundary conditions are defined by a point source located at the point \f$ \mathbf{x}_0 \f$.
By this definition the geometry is considered transparent, and the resulting acoustic fields on \f$ S \f$ are obtained as the acoustic field of the point source.

The compiled mex file can be called from Matlab, as demonstrated in the example below:

\include helmholtz_bem_2d_cyl_test.m

The generated plot looks like

\image html tut_helmholtz_bem_2d_cyl.png

And the resulting errors read as

	Log10 mean error: -2.484395

The full codes of this tutorial are available here:
- C++ code: tutorial/helmholtz_bem_2d_cyl.mex.cpp
- Matlab code: tutorial/helmholtz_bem_2d_cyl_test.m

