Fictitious eigenfrequencies with a collocational Helmholtz BEM {#tut_helmholtz_colloc_bem_fict}
==============================================================

\page tut_helmholtz_colloc_bem_fict

[BEM example]:\ref theo_bem_example
[operator notation]:\ref bem_example_op
[weighted residuals]:\ref bem_example_res
[Rayleigh integral tutorial]:\ref tut_rayleigh_integral
[function space views]:\ref tut_funcspace_const_iso_view
[dirac view]:\ref tut_funcspace_dirac_view
[TOC]

Introduction {#tut_helmholtz_colloc_bem_fict_intro}
============

This tutorial demonstrates the usage of NiHu for solving an exterior acoustic radiation problem in 3D by means of a collocational type boundary element method (BEM).
It is well known that the standard BEM formaulation, when applied to an exterior radiation problem, does not have a unique solution at the eigenfrequencies of the dual interior problem. This tutorial presents how two mitigation methods
- The CHIEF method
- and the Burton Miller formalism
are applied to handle the problem.

If you are familiar with the CHIEF and Burton Miller formalisms, you should be able to implement the C++ code based on the previous tutorial \ref tut_helmholtz_galerkin_bem. However, it is worth reading this tutorial, because it will demonstrate how singular integrals can be specialised for specific formalis, element types and kernels.

Theory {#tut_helmholtz_colloc_bem_fict_theory}
======

The eigenfrequency problem {#tut_helmholtz_colloc_bem_fict_integrals}
--------------------------

The boundary integral equation

\f$ \displaystyle
\frac{1}{2} p({\bf x}) = \left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}), \quad \bf{x} \in S
\f$

when applied to exterior radiation or scattering problems, does not have a unique solution at the eigenfrequencies of the dual interior problem.
It can be shown that the interior mode shapes \f$ p^*_k({\bf y}) \f$ and its normal derivative \f$ q^*_k({\bf y}) = \partial p^*_k({\bf y}) / \partial n_{\bf y} \f$ satisfy the boundary integral equation, although, obviously, they have no phyiscal relation to the exterior problem.

The two discussed solution methods mitigate the problem by extending the integral equation by other equations that are not satisfied by the modal ''solution''.

The CHIEF method  {#tut_helmholtz_colloc_bem_fict_chief}
----------------

The CHIEF method (the name comes from Combined Helmholtz Integral Equation Formalism) utilises that the mode shapes do not satisfy the Helmholtz integral if the source point \f$ {\bf x} \f$ is inside the interior volume, and does not coincide with a nodal location of the mode shape. Therefore, the boundary integral equation is extended by a set of additional equations as follows:

\f$
\displaystyle
\frac{1}{2} p({\bf x}) = \left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}), \quad {\bf x} \in S, \\
\displaystyle
0 = \left(\mathcal{M}p\right)_S({\bf x}_l) - \left(\mathcal{L}q\right)_S({\bf x}), \quad {\bf x}_l \in V^{\text{interior}}, l = 1 \dots K
\f$

Theoretically, one additional equation would be enough to exclude the fictitious solution. However, as it is difficult to estimate where the nodal locations of the interior problem are, the method is usually applied by selecting a number of randomly chosen interior points \f$ {\bf x}_l \f$.

After discretisation, the following overdetermined linear system of equations is to be solved:

\f$ \displaystyle
\left[ \begin{matrix} {\bf M} - {\bf I}/2 \\ {\bf M}' \end{matrix} \right]
{\bf p}
=
\left[ \begin{matrix} {\bf L} \\ {\bf L}' \end{matrix} \right]
{\bf q}
\f$


The Burton and miller Formalism  {#tut_helmholtz_colloc_bem_fict_bm}
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

The differentiated equation suffers from the problem of fictitious eigenfrequencies too. However, as the eigenfrequencies and mode shapes of the original and differentiated equations never coincide, the simultaneous solution is unique. This solution is found by solving the superposition of the two equations using an appropriate coupling constant \f$ \alpha \f$

\f$
\displaystyle
\left(\left[\mathcal{M} - \frac{1}{2}\mathcal{I} + \alpha \mathcal{N} \right]p\right)_S({\bf x})
=
\left(\left[\mathcal{L} + \alpha \mathcal{M}^{\mathrm{T}} + \alpha \frac{1}{2} \mathcal{I}\right] q\right)_S({\bf x})
\f$


In the sequel, the above Matrix equations are determined and solved using NiHu in the following steps
1. Assembling the system matrices 
2. Solution of the surface equation
3. Evaluation of the radiation equation for the field points
4. Displaying the results and quantifying the error

Step (1) is performed in C++, by means of a mex function, whereas steps (2-4) are carried out utilising the Matlab interface of NiHu.

The C++ code {#tut_helmholtz_colloc_bem_fict_cpp}
============

The objective of our C++ code is to assemble the systems matrices in a function callable from Matlab by means of `mex`.

Header and mesh creation {#tut_helmholtz_colloc_bem_fict_header}
------------------------

As it was done in the [Rayleigh integral tutorial] the header of our C++ source file includes the necessary header files and defines some basic types for convenience.

\snippet helmholtz_bem_3d.mex.cpp Header

Since we have to calculate four different matrices, that are \f$ \mathbf{M}^{(s)} \f$, \f$ \mathbf{L}^{(s)} \f$, \f$ \mathbf{M}^{(f)} \f$ and \f$ \mathbf{L}^{(f)} \f$, our mex function will return four left hand side parameters.

Similar to the [Rayleigh integral tutorial], two meshes will be created.
Thus, our mex function header and the mesh generation commands will read as

\snippet helmholtz_bem_3d.mex.cpp Mex and mesh

Definition of function spaces {#tut_helmholtz_colloc_bem_fict_spaces}
-----------------------------

Next, we create the discretised function spaces of our mesh.
In this example we will use homogeneous function spaces, which are easily created by means of [function space views] of the meshes.
On the radiating surface, a piecewise constant approximation of trial and test functions are applied, whereas on the field point mesh, a piecewise constant approximation of the trial function with Dirac test functions is applied.
The Dirac delta test functions are created by using the simple [dirac view] approach.

\snippet helmholtz_bem_3d.mex.cpp Function spaces

The number of degrees of freedom on the radiating surface \f$ S \f$ and in the field mesh \f$ F \f$ are denoted by \f$ n \f$ and \f$ m \f$, respectively.
After the function spaces are built, the number of DOFs are known, and the memory for the system matrices can be preallocated.
The surface system matrices \f$ \mathbf{M}^{(s)} \f$ and \f$ \mathbf{L}^{(s)} \f$ are of size \f$ n \times n \f$, whereas the field system matrices \f$ \mathbf{M}^{(f)} \f$ and \f$ \mathbf{L}^{(f)} \f$ are of size \f$ m \times n \f$.
The four system matrices are preallocated by means of the following lines of code.

\snippet helmholtz_bem_3d.mex.cpp Matrices

Integral operators and their evaluation {#tut_helmholtz_colloc_bem_fict_intop}
---------------------------------------

In the next steps, the thre integral operators \f$ \mathcal{M} \f$, \f$ \mathcal{L} \f$ and \f$ \mathcal{I} \f$ are defined.
The operators \f$ \mathcal{M} \f$ and \f$ \mathcal{L} \f$ are created as the double and single layer potential kernels of the three-dimensional Helmholtz equations, while operator \f$ \mathcal{I} \f$ is simply the identity operator.
Since the kernel functions depend on the wave number \f$ k \f$, which is passed as the fifth right hand side parameter of the mex function, the integral operators are created as follows.

\snippet helmholtz_bem_3d.mex.cpp Integral operators

Finally, the system matrices are obtained by the evaluation of the integral operators on the function spaces defined above.

\snippet helmholtz_bem_3d.mex.cpp System matrices

It should be realised that the second line of the above code snippet contains the evaluation two integral operators, \f$ \mathcal{M} \f$ and \f$ \mathcal{I} \f$ and stores the summed result in the matrix \f$ \mathbf{M}^{(s)} \f$.

The resulting mex function is capable of computing the four system matrices for arbitatry radiating surface meshes consisting of quadrilateral and triangular elements and field meshes consisting of quadrilateral elements only at a give wave number \f$ k \f$.

The Matlab code {#tut_helmholtz_colloc_bem_fict_matlab}
=============== 

The example creates a sphere surface of radius \f$ R = 1\,\mathrm{m} \f$ centered at the origin, consisting of linear quadrilateral and triangular elements.
The field mesh is defined as a disc segment, parametrised as \f$ 1.125\,\mathrm{m} \le r \le 3.625\,\mathrm{m} \f$, \f$ 0 \le \varphi \le \pi / 2 \f$ and \f$ \vartheta = 0 \f$ .

The Neumann boundary conditions are defined by a point source located at the point \f$ \mathbf{x}_0 \f$. 
By this definition the geometry is considered transparent, and the resulting acoustic fields both on \f$ S \f$ and \f$ F \f$ are obtained as the field of the point source.
Therefore it is easy to compare the calculated numerical results to the analytical ones.

The compiled mex file can be called from Matlab, as demonstrated in the example below:

\include helmholtz_bem_3d_test.m

The generated plot looks like

\image html tut_helmholtz_colloc_bem_fict.png

And the resulting errors read as

	Surface log10 error: -1.780843 
	Field   log10 error: -1.964462 

The full codes of this tutorial are available here:
- C++ code: tutorial/helmholtz_bem_3d.mex.cpp
- Matlab code: tutorial/helmholtz_bem_3d_test.m
