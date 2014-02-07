An introduction to NiHu's mesh representation? {#tech_geometry}
=============================================

\page tech_geometry

[TOC]

Introduction {#tech_geometry_intro}
============

This document explains how NiHu stores its geometries.


Coordinate spaces {#tech_geometry_space}
=================

The lowest level of geometrical representations is the definition of coordinate spaces like the three-dimensional real space \f$ \mathbb{R}^{3} \f$.
Coordinate spaces are defined in class ::space, implemented in the source file space.hpp.
The class has two obvious template parameters, \c scalar, representing the type of a coordinate and \c dimension representing the dimensionality of the coordinate space.

The class - like most of NiHu's template classes - is a self-returning metafunction.
This means that class ::space contains a public \c typedef called \c type that returns class ::space itself.
This technique is useful for simplifying NiHu's syntax, as will be demonstrated later.
As an other typical NiHu concept, class ::space contains its template arguments in the form of a public \c typedef and an \c enum, making the arguments available from the class.

Class ::space further defines the location vector type \c location_t of the coordinate space in the form of an Eigen vector type.


A few commonly used coordinate spaces are predefined in NiHu's component library.
These typedefs are rather straightforward:

~~~~~
typedef space<double, 1> space_1d;
typedef space<double, 2> space_2d;
typedef space<double, 3> space_3d;
~~~~~

For a simple usage example consider the following code snippet:
~~~~~
space_3d::location_t x, y;
x << 0.0, 1.0, 1.0;
y << 1.0, 1.0, 2.0;
space_3d::scalar_t d = (y - x).norm();
std::cout << "The distance between " << x << " and " << y << " is " << d << std::endl;
~~~~~

Domains {#tech_geometry_domain}
=======


Interpolation functions {#tech_geometry_shapefun}
=======================


Elements {#tech_geometry_elements}
========

