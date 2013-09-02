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

This tutorial demonstrates the usage of NiHu for solving the Helmholtz equation in 3D by means of a Galerkin type boundary element method (BEM).
In this tutorial the Helmholtz equation will be applied for the solution of an acoustic wave propagation problem.
In this tutorial we will use NiHu's Matlab interface for the mesh definition and the solution of the resulting matrix equations.
This tutorial is based on the theory and notations introduced previously in the theoretical introduction [BEM example].

Theory {#tut_helmholtz_colloc_bem_fict_theory}
======

In this example the Helmholtz integrals are utilised in order to solve an external radiation problem.
Our radiating surface is a closed surface in the three-dimensional space \f$ S \subset \mathbb{R}^3 \f$.
Neumann boundary conditions are specified on the entire surface, given as \f$ q(\mathbf{x}) = \bar{q}(\mathbf{x}) \f$ if \f$ \mathbf{x} \in S \f$.
We need to determine the radiated acoustic field \f$ p(\mathbf{x}) \f$ on the surface \f$ S \f$ and in an other set of points \f$ F \subset \mathbb{R}^3 \setminus S \f$ referred to as field points.

The Helmholtz integrals {#tut_helmholtz_colloc_bem_fict_integrals}
-----------------------

The sound pressure radiated by the closed vibrating surface \f$ S \f$ is obtained by means of the Helmholtz integrals

\f$ \displaystyle \frac{1}{2} p({\bf x}) = \int_S H({\bf x}, {\bf y}) p({\bf y}) dS_y - \int_S G({\bf x}, {\bf y}) q({\bf y}) dS_y \qquad \f$ if \f$ \quad x \in S \f$

\f$ \phantom{\displaystyle \frac{1}{2}} p({\bf x}) = \displaystyle \int_S H({\bf x}, {\bf y}) p({\bf y}) dS_y - \int_S G({\bf x}, {\bf y}) q({\bf y}) dS_y \qquad \f$ if \f$ \quad x \in F , \f$

where \f$ G({\bf x}, {\bf y}) \f$ and \f$ H({\bf x}, {\bf y}) \f$ denote the fundamental solution of the inhomogeneous Helmholtz equation and its normal derivative with respect to the source point.
The kernel functions are given as

\f$ G({\bf x}, {\bf y}) = \displaystyle \frac{\mathrm{e}^{-\mathrm{i}kr}}{4\pi r} \quad \f$
and 
\f$ \quad H({\bf x}) = \displaystyle \frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial n_y} = -G({\bf x}, {\bf y}) \left(\displaystyle \frac{1 + \mathrm{i} k r}{r} \frac{{\bf r} \cdot {\bf n}_y }{r} \right) , \f$

with the distance vector \f$ \bf r \f$ and the distance \f$ r \f$ defined as
\f$ \mathbf{r} = \mathbf{y} - \mathbf{x} \f$ 
and
\f$ r = \left| \mathbf{r} \right| . \f$
The normal vector \f$ \mathbf{n}_y \f$ is pointing outward from the surface \f$ S \f$ at the point \f$ \mathbf{y} \f$.

Galerkin BEM formalism {#tut_helmholtz_colloc_bem_fict_formalism}
----------------------

In order to solve the Helmholtz integral equations numerically, the radiating surface \f$ S \f$ is divided into a number of non-overlapping elements.
The sound pressure function \f$ p(\mathbf{x}) \f$ and its normal derivative \f$ q(\mathbf{x}) \f$ are approximated by means of weighting functions \f$ w(\mathbf{x}) \f$, also called trial functions, such that

\f$ p(\mathbf{x}) \approx \displaystyle \sum_{j=1}^n w_j (\mathbf{x}) p_j \qquad \f$
and
\f$ \qquad q(\mathbf{x}) \approx \displaystyle \sum_{j=1}^n w_j (\mathbf{x}) q_j , \f$

with \f$ p_j \f$ and \f$ q_j \f$ denoting the corresponding weights and \f$ n \f$ representing the total number of weighting functions, also known as the number of degrees of freedom.

As a next step, the weak form of the Helmholtz integral equation is obtained by multiplying the equations with a set of test functions \f$ t_j(\mathbf{x}) \f$ and integrating the result over the complete surface \f$ S \f$ with respect to the variable \f$ \mathbf{x} \f$.
In the Galerkin formulation, the test function are equivalent to the weighting functions, i.e. \f$ t_j(\mathbf{x}) \equiv w_j(\mathbf{x}) \f$.
Hence, the Galerkin formalism is a special case of weighted residuals with the weighting functions being equivalent to the trial functions.

Applying the [operator notation] introduced for [weighted residuals] , the double integral for the surface is expressed in the form

\f$ \left< w_i , \left(\left[ \mathcal{M} - \frac{1}{2} \mathcal{I} \right] w_j \right)_S \right>_S p^{(s)}_j = \left< w_i, \left(\mathcal{L} w_j \right)_S \right>_S q^{(s)}_j . \f$

The upper inidices \f$ \cdot^{(s)} \f$ and \f$ \cdot^{(f)} \f$ symbolise the variables on the surface \f$ S \f$ and the field \f$ F \f$ meshes, respectively.
For the calculation of the acoustic field of the radiated surface a Dirac-delta test function is applied herein, which gives

\f$ p^{(f)}_i = \left< \delta_i , \left( \mathcal{M} w_j \right)_S \right>_F p^{(s)}_j - \left< \delta_i, \left(\mathcal{L} w_j \right)_S \right>_F q^{(s)}_j . \f$

The equations can be rewritten by utilizing matrix notations 

\f$ M^{(s)}_{ij} = \left< w_i , \left(\left[ \mathcal{M} - \frac{1}{2} \mathcal{I} \right] w_j \right)_S \right>_S \f$,
\f$ \qquad L^{(s)}_{ij} = \left< w_i, \left(\mathcal{L} w_j \right)_S \right>_S \f$,

\f$ M^{(f)}_{ij} = \left< \delta_i , \left( \mathcal{M}  w_j \right)_S \right>_F \qquad \f$
and
\f$ \qquad L^{(f)}_{ij} = \left< \delta_i , \left( \mathcal{L}  w_j \right)_S \right>_F \f$

in the forms

\f$ \mathbf{M}^{(s)} \mathbf{p}^{(s)} = \mathbf{L}^{(s)} \mathbf{q}^{(s)} \qquad \f$ and

\f$ \phantom{\mathbf{M}^{(s)}} \mathbf{p}^{(f)} = \mathbf{M}^{(f)} \mathbf{p}^{(s)} - \mathbf{L}^{(f)} \mathbf{q}^{(s)} .\f$

The vectors \f$ \mathbf{p} \f$ and \f$ \mathbf{q} \f$ denote the column vectors of the weights for the degrees of freedom of the pressure and its normal derivative, respectively.

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
