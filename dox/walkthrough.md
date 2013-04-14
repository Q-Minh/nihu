A bottom-up walkthrough {#walkthrough}
=======================

[TOC]

Introduction
============

This bottom-up walkthrough is a short programmers' guide.  If you would like to introduce your shape functions, elements or kernels into the NiHu toolbox, the best practice is to read this guide, understand the structure of the software, and modify the code afterwards.

Main Classes 
============

class Domain
------------

Class ::domain describes a base domain \f$\mathcal{D} = \left\{\xi\right\}\f$ over which interpolation functions \f$L(\xi)\f$, \f$N(\xi)\f$ can be defined.
The class is templated on the number of spatial dimensions of the \f$\xi\f$-space and the number of domain corners.
Additionally, class ::domain defines the whole problem's real scalar type (typically double).
Class ::domain stores and can return the base domain's corner nodes and center.

In the current development state, class ::domain is specialised for
- 1D line domain (::line_domain), used in 2D BEM,
- 2D triangle domain (::tria_domain), used in 3D BEM,
- 2D quadrilateral domain (::quad_domain) used in 3D BEM,
- 3D brick domain (::brick_domain) not used at all.

These specialisations are plain typedefs.


class ShapeSet
--------------

ShapeSet classes define interpolation functions \f$L(\xi)\f$ on the base domain. All ShapeSet classes are derived from the CRTP base class ::shape_set_base. The interface base class ::shape_set_base obtains the derived class's parameters through a traits class ::shape_set_traits. The interface can return a vector of shape functions \f$L(\xi)\f$ for any input coordinate \f$\xi\f$, as well as the shape functions derivative \f$\nabla L(\xi)\f$ with respect to the domain variable \f$\xi\f$. Furthermore, the interface can return begin and end iterators to the interpolation function's nodal locations.

Specialisations of the CRTP base class are templated on base domains. The interface is specialised for two general families of shape sets:
- ::constant_shape_set. A ::constant_shape_set has only one nodal corner located at the ::domain's center. The shape function is constant 1, its derivatives are zero. Further specialisations of the ::constant_shape_set are
	+ ::line_0_shape_set
	+ ::tria_0_shape_set
	+ ::quad_0_shape_set
- Other important family is the family of isoparametric shape sets (::isoparam_shape_set). An isoparametric shape set's corner nodes are located at the base domain's corners. The shape functions are defined so that each shape function has the value 1 at its corner node and zero in other corner nodes. Further specialisations of the ::isoparam_shape_set are
	+ ::line_1_shape_set
	+ ::tria_1_shape_set
	+ ::quad_1_shape_set

- As a special shape set, NiHu implements ::parallelogram_shape_set. A parallelogram shape set has three corner nodes, but it is defined over a quad base domain. This shape set is currently not used in the code.
- Further shape sets can be introduced by implementing their functions that return shape functions, its derivatives and the nodal corner points.
	

class Element
-------------

An element is a geometrical interpolation function set \f$L_k(\xi)\f$ defined over a base domain \f$\mathcal{D}\f$ and a set of nodal coordinates \f$x_k\f$.

