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
These domains are used in NiHu as the reference domains \f$ \mathcal{D} \f$ of the coordinate transforms defining element geometries.

The class is a template with two arguments:
~~~~~
template <class Space, unsigned NumCorners>
class domain;
~~~~~
\c Space describes the coordinate space of the domain, and \c NumCorners describes the number of domain corners.
The class is a self-returning metafunction, and defines its template arguments as a nested type and an enum.

The class further stores its corners (locations in the parameter space), as well as the domain's center in static members.
The class provides static methods to return
- the begin iterator to the domain's corners array ::domain::get_corners
- the domain's specific corner ::domain::get_corner
- the domain's center ::domain::get_center
- and the domain's volume ::domain::get_volume

Each domain is assigned a numeric identifier by the metafunction ::domain_traits::id. The default value of this id is

	10 * dimension + NumCorners

thus the id of \c line_domain is 12, and the id of tria_domain is 23.
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

Core definition {#tech_geometry_shapeset_core}
---------------

Interpolation functions or shape functions \f$ L_i(\xi) \f$ are used to define the element's coordinate transform. \f$ x(\xi) = \sum x_i L_i(\xi) \f$.
They are implemented in the header file shapeset.hpp.

Shape function sets are implemented using the CRTP pattern with traits metafunctions.
All shape set classes are derived from shape_set_base that defines their interface.
The traits of a derived shape set class are implemented in namespace ::shape_set_traits in the form of separate metafunctions. 
The base class template ::shape_set_base gets the necessary type information from the traits metafunctions of the derived classes to define the specific interface.

The interface consists of the following functions:
- shape_set_base::corner_begin and shape_set_base::corner_end return iterators to the nodal corners \f$ \xi_i \f$ of the shape set.
- shape_set_base::corner_at returns a specified corner of the shape set.
- shape_set_base::eval_shape evaluates and returns the shape functions or their derivatives at a local coordinate.

The function template eval_shape can return the shape functions their gradient vectors or second derivatives as follows:
- ::shape_set_base::eval_shape<0> returns the shape functions \f$ L_i(\xi) \f$ in the form of an Eigen vector.
- ::shape_set_base::eval_shape<1> returns the gradient of each function \f$ L_{i,j}(\xi) \f$ in the form of an Eigen matrix, each row containing the gradient of a single shape function.
- ::shape_set_base::eval_shape<2> returns the second derivatives of each function \f$ L_{i,jk}(\xi) \f$ in the form of an Eigen matrix, each row containing the different second derivatives of a shape function.

Efficient shape function evaluation
-----------------------------------

For increased efficiency, NiHu checks _during compilation_ if the shape functions and derivatives are
- really dependent on the local coordinate \f$ \xi \f$ (the general case)
- constant (like for the case of the gradient matrix of a linear interpolation)
- or zero (like the second derivative matrix of a linear interpolation)

For the case of a general \f$ \xi \f$-dependent shape function, the result is computed on the fly and is returned by value.
For the case of a constant function (like the gradient of linear interpolation functions) the shape functions are precomputed during program initialisation, stored in static members, and the function returns a constant reference.
For the case of a zero matrix (like the second derivatives of linear interpolation functions) an efficient Eigen::Zero type matrix is returned.

The value types of the shape function vectors and matrices, as well as the actual return types of function eval_shape are defined in the metafunctions shape_set_traits::shape_value_type and shape_set_traits::shape_return_type.

Further traits
--------------

Similar to class ::space and ::domain, each shape set is assiged a numeric and a textual identifier by the metafunctions ::shape_set_traits::id and ::shape_set_traits::name. The numeric id has a default value of

	100 * domain_traits::id + NumNodes
	
thus, the automatic id of a 2D 3-noded triangle interpolation function set is 2303.

Library reference {#tech_geometry_shapeset_lib}
-----------------

NiHu's component library predefines the following shape functions:
- ::line_0_shape_set constant 1-noded line interpolation (1201)
- ::line_1_shape_set linear 2-noded line interpolation (1202)
- ::line_2_shape_set quadratic 3-noded line interpolation (1203)
- ::tria_0_shape_set constant 1-noded triangle interpolation (2301)
- ::tria_1_shape_set linear 3-noded triangle interpolation (2303)
- ::tria_2_shape_set quadratic 6-noded triangle interpolation (2306)
- ::quad_0_shape_set constant 1-noded quadrangle interpolation (2401)
- ::quad_1_shape_set (bi)linear 4-noded quadrangle interpolation (2404)
- ::quad_28_shape_set quadratic 8-noded quadrangle interpolation (2409)
- ::quad_2_shape_set quadratic 9-noded quadrangle interpolation (2408)
- ::brick_0_shape_set constant 1-noded hexahedron interpolation (3801)
- ::brick_1_shape_set (tri)linear 8-noded hexahedron interpolation (3808)

Example {#tech_geometry_shapeset_example}
-------

Usage of interpolation functions is demonstrated with the following example:
~~~~~
#include "core/shapeset.hpp"
#include "library/quad_2_shape_set.hpp"

int main(void)
{
	double const step = 2e-1;
	Eigen::Matrix<double, 9, 1> nodal_data;
	nodal_data << 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	for (double xi = -1.0; xi < 1.0; xi += step)
	{
		for (double eta = -1.0; eta < 1.0; eta += step)
		{
			quad_2_shape_set::xi_t _xi(xi, eta);
			auto L = quad_2_shape_set::eval_shape<0>(_xi);
			auto value = L.transpose() * nodal_data;
			std::cout << "F(" << _xi.transpose() << ") = " << value << std::endl;
		}
	}
	return 0;
}
~~~~~
The example interpolates a biquadratic function over the standard quad domain and prints out its values.

Elements {#tech_geometry_elements}
========

