Customising singular integrals {#tut_custom_singular_integrals}
==============================

\page tut_custom_singular_integrals

[BEM example]:\ref theo_bem_example
[operator notation]:\ref bem_example_op
[weighted residuals]:\ref bem_example_res
[Rayleigh integral tutorial]:\ref tut_rayleigh_integral
[function space views]:\ref tut_funcspace_const_iso_view
[dirac view]:\ref tut_funcspace_dirac_view
[TOC]

Introduction {#tut_custom_singular_integrals_intro}
============

This tutorial explains how to customise the evaluation of a singular integral in NiHu. Our examples of demonstration are the collocational singular integrals of the 3D Helmholtz kernels on planar triangle elements.

Theory {#tut_custom_singular_integrals_theory}
======

The singular integrals and the method of static subtraction {#tut_custom_singular_integrals_subtraction}
-----------------------------------------------------------

The collocational singular integral of the single layert potential kernel on a constant triangular element reads as

\f$
I = \int_S \frac{1}{4\pi |{\bf y} - {\bf x_0}|}\mathrm{d}S_{\bf y}
\f$

where \f$ {\bf x}_0 \f$ is the singular collocation point in the center of the element.

The Helmholtz boundary integral equation

\f$ \displaystyle
\frac{1}{2} p({\bf x}) = \left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}), \quad \bf{x} \in S
\f$

when applied to exterior radiation or scattering problems, does not have a unique solution at the eigenfrequencies of the dual interior problem.
It can be shown that the interior mode shapes \f$ p^*_k({\bf y}) \f$ and its normal derivative \f$ q^*_k({\bf y}) = \partial p^*_k({\bf y}) / \partial n_{\bf y} \f$ satisfy the boundary integral equation, although, obviously, they have no phyiscal relation to the exterior problem.

The two discussed solution methods mitigate the problem by extending the integral equation by other equations that are not satisfied by the modal solution.

The CHIEF method  {#tut_custom_singular_integrals_chief}
----------------

The CHIEF method (the name comes from Combined Helmholtz Integral Equation Formalism) utilises that the mode shapes \f$ p^*_k({\bf y}) \f$ do not satisfy the Helmholtz integral with source point \f$ {\bf x} \f$ inside the interior volume, if the source point coincides with a nodal location of the mode shape.
Therefore, the boundary integral equation is extended by a set of additional equations as follows:

\f$
\displaystyle
\frac{1}{2} p({\bf x}) = \left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}), \quad {\bf x} \in S, \\
\displaystyle
0 = \left(\mathcal{M}p\right)_S(\tilde{\bf x}_l) - \left(\mathcal{L}q\right)_S(\tilde{\bf x}_l), \quad \tilde{\bf x}_l \in V^{\text{interior}}, \quad l = 1 \dots K
\f$

Theoretically, one additional equation would be enough to exclude the fictitious solution.
However, as it is difficult to estimate where the nodal locations of the interior problem are, the method is usually applied by selecting a number of randomly chosen interior points \f$ {\bf x}_l \f$.

After discretisation, the following overdetermined linear system of equations is to be solved:

\f$ \displaystyle
\left[ \begin{matrix} {\bf M} - {\bf I}/2 \\ {\bf M}' \end{matrix} \right]
{\bf p}
=
\left[ \begin{matrix} {\bf L} \\ {\bf L}' \end{matrix} \right]
{\bf q}
\f$

where

\f$
\displaystyle  L_{ij} = \left< \delta_{{\bf x}_i}, \left(\mathcal{L} w_j\right)_S \right>_S, \quad
\displaystyle  M_{ij} = \left< \delta_{{\bf x}_i}, \left(\mathcal{M} w_j\right)_S \right>_S \\
\f$

are the conventional system matrices (\f$ w_j \f$ denotes the weighting shape function on the surface) of a collocational BEM, and

\f$
\displaystyle  L'_{lj} = \left< \delta_{\tilde{\bf x}_l}, \left(\mathcal{L} w_j\right)_S \right>_S, \quad
\displaystyle  M'_{lj} = \left< \delta_{\tilde{\bf x}_l}, \left(\mathcal{M} w_j\right)_S \right>_S
\f$

are the additional matrices originating from the discretised CHIEF equations.


The Burton and miller Formalism  {#tut_custom_singular_integrals_bm}
-------------------------------

The Burton and Miller method searches for the solution of the original boundary integral equation and its normal derivative with respect to the normal at \f$ {\bf x} \f$:

\f$
\displaystyle
\frac{1}{2} p({\bf x}) = \left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}), \quad {\bf x} \in S, \\
\displaystyle
\frac{1}{2} q({\bf x}) = \left(\mathcal{N}p\right)_S({\bf x}) - \left(\mathcal{M}^{\mathrm{T}}q\right)_S({\bf x}), \quad {\bf x} \in S
\f$

where

\f$
\displaystyle
\left(\mathcal{M}^{\mathrm{T}}q\right)_S({\bf x}) = \frac{\partial}{\partial n_{\bf x}}\int_S G({\bf x}, {\bf y}) q({\bf y})\mathrm{d}S_{\bf y}
\f$

is the transpose of \f$ \left(\mathcal{M}q\right)({\bf x})\f$, as differentiation is performed with respect to the other variable, and

\f$ \displaystyle
\left(\mathcal{N}p\right)_S({\bf x}) = \frac{\partial}{\partial n_{\bf x}}\int_S \frac{\partial}{\partial n_{\bf y}} G({\bf x}, {\bf y}) p({\bf y})\mathrm{d}S_{\bf y}
\f$

is the hypersingular operator.

The differentiated equation suffers from the problem of fictitious eigenfrequencies too.
However, as the eigenfrequencies and mode shapes of the original and differentiated equations never coincide, the simultaneous solution is unique.
This solution is found by solving the superposition of the two equations using an appropriate coupling constant \f$ \alpha \f$

\f$
\displaystyle
\left(\left[\mathcal{M} - \frac{1}{2}\mathcal{I} + \alpha \mathcal{N} \right]p\right)_S({\bf x})
=
\left(\left[\mathcal{L} + \alpha \mathcal{M}^{\mathrm{T}} + \alpha \frac{1}{2} \mathcal{I}\right] q\right)_S({\bf x}),
\quad {\bf x} \in S
\f$

The collocational discretisation of the above equation is straightforward.

However, we should take into consideration that the hypersingular operator's kernel has an \f$ O(1/r^3) \f$ type singularity, and needs special treatment.
This topic is discussed in details in the tutorial \ref tut_custom_singular_integrals


Program structure {#tut_custom_singular_integrals_structure}
=================

We are going to implement a Matlab-C++ NiHu application, where
- the C++ executable is responsible for assembling the system matrices from the boundary surface mesh and the CHIEF points,
- the Matlab part defines the meshes, calls the C++ executable, solves the systems of equations, and compares the solutions by quantifying their error

The C++ code {#tut_custom_singular_integrals_cpp}
============

The C++ code is going to be called from Matlab as

	[L, M, Lchief, Mchief, Mtrans, N] = helmholtz_bem_3d_fict(surf_nodes, surf_elements, chief_nodes, chief_elements, wave_number);
	
\note The CHIEF points are defined as element centres of a NiHu mesh for convenience.
With this choice, the applied discretisation formalism is identical to that of a radiation integral in tutorial \ref tut_rayleigh_integral.

\note The computed system matrices `M` and `Mtrans` will contain the discretised identity operator too, as required by the conventional and Burton-Miller formalisms.

Header and mesh creation {#tut_custom_singular_integrals_header}
------------------------

The header of our C++ source file includes the necessary header files and defines some basic types for convenience.
We further include library/helmholtz_singular_integrals.hpp that defines the specialised singular integrals of the Helmholtz kernel for constant triangles.

\snippet helmholtz_bem_3d_fict.mex.cpp Header

Two meshes will be created:
- One for the radiator surface (`surf_mesh`) that contains only triangular elements
- One for the CHIEF points (`chief_mesh`) that contains only quadrangle elements
  \note Only the element centres of the CHIEF point mesh are important

\snippet helmholtz_bem_3d_fict.mex.cpp Mesh

Definition of function spaces {#tut_custom_singular_integrals_spaces}
-----------------------------

Next, we create the discretised function spaces of our mesh.
The definitions are straightforward, as only constant triangular boundary elements are dealt with.
For the definition of the CHIEF mesh, we immediately apply the dirac view.

\snippet helmholtz_bem_3d_fict.mex.cpp Function spaces

The number of degrees of freedom on the radiating surface \f$ S \f$ and in the field mesh \f$ F \f$ are denoted by \f$ n \f$ and \f$ m \f$, respectively.
After the function spaces are built, the number of DOFs are known, and the memory for the system matrices can be preallocated.
The surface system matrices \f$ \mathbf{L} \f$, \f$ \mathbf{M} \f$, , \f$ \mathbf{M}_\mathrm{transpose} \f$ and \f$ \mathbf{N} \f$ are of size \f$ n \times n \f$, whereas the matrices describing the CHIEF equations \f$ \mathbf{M}' \f$ and \f$ \mathbf{L}' \f$ are of size \f$ m \times n \f$.
The six system matrices are preallocated by means of the following lines of code.

\snippet helmholtz_bem_3d_fict.mex.cpp Matrices

Integral operators and their evaluation {#tut_custom_singular_integrals_intop}
---------------------------------------

In the next steps, the five integral operators \f$ \mathcal{L} \f$, \f$ \mathcal{M} \f$, \f$ \mathcal{M}^\mathrm{T} \f$, \f$ \mathcal{N} \f$ and \f$ \mathcal{I} \f$ are defined.
Since the kernel functions depend on the wave number \f$ k \f$, which is passed as the fifth right hand side parameter of the mex function, the integral operators are created as follows.

\snippet helmholtz_bem_3d_fict.mex.cpp Integral operators

Finally, the system matrices are obtained by the evaluation of the integral operators on the function spaces defined above.
Note that as a collocational formalism is applied, the Dirac-view of the test function spaces is taken for the boundary integrals.

\snippet helmholtz_bem_3d_fict.mex.cpp System matrices

\note It should be realised that the evaluation of matrices `M_surf` and `Mt_surf` contain the evaluation of two integral operators, \f$ \mathcal{M} \f$ and \f$ \mathcal{I} \f$, or \f$ \mathcal{M}^{\mathrm{T}} \f$ and \f$ \mathcal{I} \f$ and stores the summed results.

The Matlab code {#tut_custom_singular_integrals_matlab}
=============== 

The example creates a sphere surface of radius \f$ R = 1\,\mathrm{m} \f$ centred at the origin, consisting of triangular elements only.
The CHIEF points are located on the surface of a cube, shifted from the origin.

The Neumann boundary conditions are defined by a point source located at the point \f$ \mathbf{x}_0 \f$.
By this definition the geometry is considered transparent, and the resulting acoustic fields on \f$ S \f$ are obtained as the acoustic field of the point source.

The compiled mex file can be called from Matlab, as demonstrated in the example below:

\include helmholtz_bem_3d_fict_test.m

The generated plot looks like

\image html tut_custom_singular_integrals.png

And the resulting errors read as

	Conventional  log10 error:  0.025 
	CHIEF         log10 error: -2.344 
	Burton-Miller log10 error: -1.680 

The full codes of this tutorial are available here:
- C++ code: tutorial/helmholtz_bem_3d_fict.mex.cpp
- Matlab code: tutorial/helmholtz_bem_3d_fict_test.m
