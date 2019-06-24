Working with elements {#tut_elements}
=====================================

\page tut_elements

[TOC]

Introduction {#tut_elements_intro}
============

The purpose of this tutorial is to explain how you can work with finite elements in NiHu.

Volume elements {#tut_elements_volume}
===============

Elements are coordinate transforms \f$ {\bf x}({\bf \xi})\f$ that map from an intrinsic domain to the physical coordinate space.
Therefore, elements are defined by an interpolation shape function set \f$ L_i(\xi) \f$ and nodal locations \f$ {\bf x}_i \f$.

A two-dimensional four noded quadrilateral volume element type can be defined as

\snippet elements.cpp Elem type

where NiHu::quad_1_shape_set refers to the four noded linear quadrilateral shape functions.
The second template argument defines the scalar type of the physical coordinate space.

After the element type has been defined, an element instance can be constructed by providing the nodal locations to the constructor.
The nodal locations need to be defined in a matrix of type NiHu::volume_element::coords_t.
In our particular case, this is a 2x4 matrix of double elements, where each column stores a nodal coordinate \f$ {\bf x}_i \f$.
The matrix is filled using Eigen's comma initialiser syntax.

\snippet elements.cpp Elem coordinates

Apparently, our element is a square volume \f$ 0 \le x,y \le 1 \f$.

The local intrinsic coordinate type of the element can be accessed as NiHu::volume_element::xi_t.

\snippet elements.cpp Local coordinate

The element's physical location can be obtained by calling the member NiHu::volume_element::get_x with an argument \f$ \xi \f$.
In our case, the location is a 2x1 double vector.

\snippet elements.cpp location

The element location gradient can be obtained by calling member NiHu::volume_element::get_dx.
In our case, this gradient is a 2x2 matrix, where the first column contains the derivative with respect to the first local coordinate \f$ {\bf x},_{\xi} \f$ and the second colum contains \f$ {\bf x},_{\eta} \f$.
The element Jacobian can be obtained by evaluating the determinant of the gradient matrix.

\snippet elements.cpp gradient

The second derivative of the coordinate mapping is returned by function NiHu::volume_element::get_ddx.
In our case, this is a 2x3 matrix, where the first column contains the second derivative with respect to the fisrt intrinsic coordinate \f$ {\bf x},_{\xi\xi}\f$, the second column contains \f$ {\bf x},_{\xi\eta}\f$, and the third colum contains \f$ {\bf x},_{\eta\eta} \f$.

\snippet elements.cpp second derivative


Surface elements {#tut_elements_surface}
================

Surface elements are handled in a similar fashion.
The example below defines a quadrilateral surface element in the 3-dimensional physical space.
Note that the physical space's dimensionality is automatically deduced from the dimensionality of the intrinsic domain.

\snippet elements.cpp surface elem type

For the case of the surface element, the location gradient is a 3x2 matrix.
The surface elements provides the function NiHu::surface_element::get_normal to return the normal vector.
Note that this normal vector is not the unit normal.
In three dimensions, its definition is \f$ {\bf n} = {\bf x},_{\xi} \times {\bf x},_{\eta}\f$.
The Jacobian can be obtained as the norm of the normal vector.
	
\snippet elements.cpp surface gradient

The second derivative is obtained similarly like for volume elements.

\snippet elements.cpp surface second derivative


