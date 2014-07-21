Working with quadratures {#tut_quadrature}
==========================================

\page tut_quadrature

[TOC]

Introduction {#tut_quadrature_intro}
============

The purpose of this tutorial is to explain how you can work with quadratures in NiHu.

Defining and traversing a Gaussian quadrature {#tut_quadrature_traversing}
=============================================

1D example {#tut_quadrature_traversing_1d}
----------

Quadratures are used to perform numerical integration over intrinsic domains.
The most common regular quadratures in finite and boundary element methods are Gaussian quadratures.
A Gaussian quadrature can be defined by specialising the template class ::gaussian_quadrature with a specific domain class, and providing a quadrature order to the constructor:

\snippet quadrature.cpp define

In the above example, a line quadrature is created that integrates over the domain \f$ -1 \le \xi \le +1\f$ with a 7-th order quadrature.
A 7-th order quadrature can accurately integrate polynomials of order up to 7.
For the case of Gaussian quadratures, an \f$ n \f$-th order quadrature consists of \f$ (n+1)/2 \f$ quadrature points.

A quadrature is a standard container of ::quadrature_elem instances, where each ::quadrature_elem contains a base point and a weight.
As a consequence, a quadrature can be traversed using standard C++ traversing methods. For example, the quadrature locations and weights can be printed as

\snippet quadrature.cpp traverse

or, equivalently:

\snippet quadrature.cpp traverse 2

and the result reads as:

	-0.861136	0.347855
	-0.339981	0.652145
	0.339981	0.652145
	0.861136	0.347855

2D example {#tut_quadrature_traversing_2d}
----------

2 and 3-dimensional quadratures, as well as quadratures over triangles are defined similarly:

\snippet quadrature.cpp 3D quad

with the output

	-0.861136 -0.861136	0.121003
	-0.861136 -0.339981	0.226852
	-0.861136  0.339981	0.226852
	-0.861136  0.861136	0.121003
	-0.339981 -0.861136	0.226852
	-0.339981 -0.339981	0.425293
	-0.339981  0.339981	0.425293
	-0.339981  0.861136	0.226852
	 0.339981 -0.861136	0.226852
	 0.339981 -0.339981	0.425293
	0.339981 0.339981	0.425293
	0.339981 0.861136	0.226852
	 0.861136 -0.861136	0.121003
	 0.861136 -0.339981	0.226852
	0.861136 0.339981	0.226852
	0.861136 0.861136	0.121003

