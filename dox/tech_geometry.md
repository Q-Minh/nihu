An introduction to NiHu's mesh representation {#tech_geometry}
=============================================

\page tech_geometry

[TOC]

Introduction {#tech_geometry_intro}
============

This document explains how NiHu stores its geometries.


Coordinate spaces {#tech_geometry_space}
=================

Core definition {#tech_geometry_space_core}
---------------

The lowest level of geometrical representations is the definition of coordinate spaces like the three-dimensional real space \f$ \mathbb{R}^{3} \f$.
Coordinate spaces are defined in class ::space, implemented in the source file space.hpp.
The class has two obvious template parameters, \c scalar, representing the type of a coordinate and \c dimension representing the dimensionality of the coordinate space.
~~~~~
template <class scalar, unsigned dimension>
class space;
~~~~~

Class ::space - like most of NiHu's template classes - implements a self-returning metafunction.
This means that the class contains a public \c typedef called \c type that returns class ::space itself.
This technique is useful for simplifying NiHu's TMP syntax, as will be demonstrated later.
As an other typical NiHu concept, class ::space contains its template arguments in the form of a public \c typedef and an \c enum, making the arguments available from the class.

Class ::space further defines the location vector type \c location_t of the coordinate space in the form of an Eigen vector type.

Library reference {#tech_geometry_space_lib}
-----------------

A few commonly used coordinate spaces are predefined in NiHu's component library.
These typedefs are rather straightforward:

~~~~~
typedef space<double, 1> space_1d;
typedef space<double, 2> space_2d;
typedef space<double, 3> space_3d;
~~~~~

Example {#tech_geometry_space_example}
-------

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

Core definition {#tech_geometry_domain_core}
---------------

Class ::domain, implemented in domain.hpp represents polygon-shaped subdomains of a coordinate space.
These domains are used in NiHu to define the parameter domains \f$ \mathcal{D} \f$ of coordinate transforms that describe element geometries.

The class is a template with two arguments:
~~~~~
template <class Space, unsigned NumCorners>
class domain;
~~~~~
\c Space describes the coordinate space of the domain, and \c NumCorners describes the number of domain corners.
The class is a self-returning metafunction, and similar to class ::space, it defines its template arguments as a nested type and an enum.

The class further stores its corners (locations in the parameter space), as well as the domain's center in static members.
The class provides static methods to return
- the begin iterator to the domain's corners array ::domain::get_corners
- the domain's specific corner ::domain::get_corner
- the domain's center ::domain::get_center
- and the domain's volume ::domain::get_volume

Each domain is assigned a numeric identifier by the metafunction ::domain_traits::id. The default value of this id is 10 \c dimension + NumCorners, this the id of \c line_domain is 12, and the id of tria_domain is 23.
Each domain is assigned a textual identifier by metafunction ::domain_traits::name.
This textual id is useful for debugging and performance diagnostics.

Library reference {#tech_geometry_domain_lib}
-----------------

The component library contains the definition of a few commonly used standard domains:
~~~~~
typedef domain<space_1d, 2> line_domain;
typedef domain<space_2d, 3> tria_domain;
typedef domain<space_2d, 4> quad_domain;
typedef domain<space_3d, 8> brick_domain;
~~~~~
with the textual identifiers
~~~~~
std::string const name<line_domain>::value = "1D Line domain";
std::string const name<tria_domain>::value = "2D Tria domain";
std::string const name<quad_domain>::value = "2D Quad domain";
std::string const name<brick_domain>::value = "3D Brick domain";
~~~~~

Interpolation functions {#tech_geometry_shapefun}
=======================


Elements {#tech_geometry_elements}
========

