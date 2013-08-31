Integral operators and weighted residuals {#tut_integral_operators}
=========================================

\page tut_integral_operators

[TOC]

Introduction {#tut_intop_intro}
============

This tutorial explains how to implement weighted double integral matrices of the form

\f$
\displaystyle
W_{ij} =
\left<t,\left(\mathcal{K}d\right)_S\right>_F =
\int_{F} t_i({\bf x}) \int_{S} K({\bf x}, {\bf y}) d_j({\bf y}) \mathrm{d}S_y \mathrm{d}S_x
\f$

in the NiHu library.

Integral operators  {#tut_intop_intop}
==================

Integral operators are defined by their kernel function. NiHu creates integral operators using the function ::create_integral_operator, called with a kernel function instance. The above lines, for example
~~~~~~~~~~~
#include "library/laplace_kernel"
auto K = create_integral_operator(laplace_3d_DLP_kernel());
~~~~~~~~~~~
instantiates the double layer potential kernel of the Laplace equation in 3D from the library, and transform the kernel to an integral operator.

A special integral operator is the identity operator with the kernel function \f$ I({\bf x}, {\bf y}) = \delta({\bf y}-{\bf x})\f$. This is instantiated without argument, as
~~~~~~~~~~~
auto I = identity_integral_operator();
~~~~~~~~~~~

Integral operations {#tut_intop_intoperation}
===================

Projection {#tut_intop_projection}
----------

The most important operation with an integral operator is letting it act on a function \f$ d({\bf y}) \f$ as:

\f$
\displaystyle
p({\bf x}) = (\mathcal{K}d)_S({\bf x})
\f$

NiHu implements this operation by indexing the integral operator with a function space. The integration domain \f$ S \f$ (the mesh) is included in the function space's definition:
~~~~~~~~~~~~~~
auto const &d = constant_view(my_mesh);
auto K = create_integral_operator(my_kernel_instance);
auto p = K[d];
~~~~~~~~~~~~~~

Inner product {#tut_intop_innerproduct}
-------------

The inner product

\f$
\displaystyle
W = \left< t, p \right> _F= \int_{F} t({\bf x}) p({\bf x}) \mathrm{d} F_x
\f$

is implemented as a multiplication between a function space and the projection:
~~~~~~~~~~~
aut const &t = constant_view(other_mesh);
auto W = t * K[d]; 		// or	W = t * p;
~~~~~~~~~~~
where the integration domain \f$ F \f$ is contained by the function space `t`;

Finally, the result of the weighted double integral can be evaluated into a matrix using the `<<` operator:
~~~~~~~~~~~
myMatrix << ( t * K[d] );	// or	myMatrix << W;
~~~~~~~~~~~
or simply
~~~~~~~~~~~
myMatrix << t * K[d];	// * is stronger than <<   :)
~~~~~~~~~~~

The inner product multiplication only works between a function space and a projection.
If we need to compose the inner product of two function spaces (defined on the same mesh, of course),
then the identity operator needs to be applied:
~~~~~~~~~~~
auto const &v = constant_view(my_mesh);
auto const &w = isoparametric_view(my_mesh);
auto I = identity_integral_operator();
matrix << v * I[w];
~~~~~~~~~~~