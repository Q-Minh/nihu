A bottom-up walkthrough {#walkthrough}
=======================

[CRTP]:http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

[TOC]


Introduction {#intro}
============

This bottom-up walkthrough is a short programmers' guide.  If you would like to introduce your shape functions, elements or kernels into the NiHu toolbox, the best practice is to read this guide, understand the structure of the software, and modify the code afterwards.

Main Classes {#classes}
============

class Domain {#domain}
------------

Implemented in domain.hpp

Class ::domain describes a base domain \f$\mathcal{D} = \left\{\xi\right\}\f$ over which interpolation functions \f$L(\xi)\f$ can be defined.
The class is templated on the number of spatial dimensions of the \f$\xi\f$-space and the number of domain corners.
Additionally, class ::domain defines the whole problem's real scalar type (typically double).
Class ::domain stores and can return the base domain's corner nodes and center.

In the current development state, class ::domain is specialised for
- 1D line domain (::line_domain), used in 2D BEM,
- 2D triangle domain (::tria_domain), used in 3D BEM,
- 2D quadrilateral domain (::quad_domain) used in 3D BEM,
- 3D brick domain (::brick_domain) not used at all.

These specialisations are plain typedefs.


class ShapeSet {#shapeset}
--------------

Implemented in shapeset.hpp

ShapeSet classes define interpolation functions \f$L(\xi)\f$ on the base domain. All ShapeSet classes are derived from the [CRTP] base class ::shape_set_base. The interface can return a vector of shape functions \f$L(\xi)\f$ for any input coordinate \f$\xi\f$, as well as the shape functions derivative \f$\nabla L(\xi)\f$ with respect to the domain variable \f$\xi\f$. Furthermore, the interface can return begin and end iterators to the interpolation function's nodal locations.

Specialisations of the [CRTP] base class are templated on base domains. The interface is generally specialised for two families of shape sets:
- ::constant_shape_set. A ::constant_shape_set has only one nodal corner located at the ::domain's center. The shape function is constant 1, its derivatives are zero. Further specialisations of the ::constant_shape_set are
	+ ::line_0_shape_set
	+ ::tria_0_shape_set
	+ ::quad_0_shape_set
- ::isoparam_shape_set. An isoparametric shape set's corner nodes are located at the base domain's corners. The shape functions are defined so that each shape function has the value 1 at its corner node and zero in other corner nodes. Further specialisations of the ::isoparam_shape_set are
	+ ::line_1_shape_set
	+ ::tria_1_shape_set
	+ ::quad_1_shape_set

- As a special shape set, ::parallelogram_shape_set is implemented. A parallelogram shape set has three corner nodes, but it is defined over a quad base domain. This shape set is currently not used in the code.
- Further shape sets (quadratic tria and quad elements) that do not fit into the constant or isoparametric families can be introduced separately by implementing their functions that return shape functions, its derivatives and the nodal corner points.
	

class Element {#element}
-------------

Implemented in element.hpp

An element is described by a a set of nodal coordinates \f$x_k\f$ and a geometrical interpolation function set \f$L_k(\xi)\f$ defined over a base domain \f$\mathcal{D}\f$.

Class ::element's template prameters are its shape function set and the number of dimensions of the \f$x\f$ space. The class stores the element's nodal coordinates. The class is capable to compute the element center, as well as the element location, location gradient and unit normal at any internal point.

As elements are usually located inside a mesh, the element class also stores the element's element id and its nodal indices too. These stored quantities are rarely used on element level, but they are essential in forming constant and isoparametric fields from elements.

Class element is specialised for the following frequently used element types with plain typedefs:
- ::line_1_elem is a linear 2 noded line in 2D space
- ::tria_1_elem is a linear 3 noded triangle in 3D space
- ::parallelogram_elem is a linear 3 noded parallelogram elem in 3D space
- ::quad_1_elem is a linear 4 noded quadrilateral element in 3D space

class Field {#field}
-----------

Implemented in field.hpp

A field is an element extended with shape functions \f$N(x)\f$ that describe some physical quantity's variation over the element geometry. Although shape functions are generally independent from the geometrical element representation, in most common cases the shape functions are related to the element geometry. In the present implementation, NiHu supports two kind of shape functions over elements:
- Constant fields describe a constant quantity over an element geometry
- Isoparametric fields describe a quantitiy that is interpolated with the element's geometrical shape functions.

As these fields do not require any additional information, only the element itself, NiHu fields are are simply elements tagged with either constant_field or isoparametric_field tags.
