Defining a custom field {#tut_custom_element}
=======================

\page tut_custom_element Defining a custom field

[metafunction]:\ref tech_meta_metafun

[TOC]

Introduction {#tut_custom_element_intro}
============

The purpose of this tutorial is to demonstrate how a custom element/field type can be introduced into the NiHu toolbox.

We define a custom 3-dimensional four noded quadrangle surface field.
The geometrical representation will be a standard 4-noded linear quadrilateral element (NiHu::quad_1_elem) having its nodes at the four element corners.
The field is going to be defined so that its shape function nodes are located at the Gaussian nodes of the element.
The new field is shown in the figure below.

![Geometrical (red) and field (blue) nodal locations of the four noded quadrilateral element in the reference domain](gaussian_field.svg)

In the following sections we define the new shape function set, introduce the new field type, and present an example how the new field type is used in a collocational or Galerkin type BEM.

The Gaussian shape set {#tut_custom_element_shapeset}
======================

The Gaussian shape set is defined over the quadrilateral domain NiHu::quad_domain located between coordinates \f$(-1,-1)\f$ and \f$(+1,+1)\f$.

The four nodal locations are

\f$ \displaystyle \left(\xi_i,\eta_i\right) = \left(\pm\sqrt{3}/3, \pm\sqrt{3}/3\right)\f$

The shape functions are defined as

\f$ \displaystyle L_1(\xi,\eta) = (1-\sqrt{3}\xi)(1-\sqrt{3}\eta)/4 \\
\displaystyle L_2(\xi,\eta) = (1+\sqrt{3}\xi)(1-\sqrt{3}\eta)/4 \\
\displaystyle L_3(\xi,\eta) = (1+\sqrt{3}\xi)(1+\sqrt{3}\eta)/4 \\
\displaystyle L_4(\xi,\eta) = (1-\sqrt{3}\xi)(1+\sqrt{3}\eta)/4\f$

and their derivatives with respect to both variables are

\f$ \displaystyle L'_{1,\xi}(\xi,\eta) = (-\sqrt{3})(1-\sqrt{3}\eta)/4, \quad L'_{1,\eta}(\xi,\eta) = (1-\sqrt{3}\xi)(-\sqrt{3})/4 \\
\displaystyle L'_{2,\xi}(\xi,\eta) = (+\sqrt{3})(1-\sqrt{3}\eta)/4, \quad L'_{1,\eta}(\xi,\eta) = (1-\sqrt{3}\xi)(-\sqrt{3})/4 \\
\displaystyle L'_{3,\xi}(\xi,\eta) = (+\sqrt{3})(1+\sqrt{3}\eta)/4, \quad L'_{1,\eta}(\xi,\eta) = (1-\sqrt{3}\xi)(-\sqrt{3})/4 \\
\displaystyle L'_{4,\xi}(\xi,\eta) = (-\sqrt{3})(1+\sqrt{3}\eta)/4, \quad L'_{1,\eta}(\xi,\eta) = (1-\sqrt{3}\xi)(-\sqrt{3})/4\f$

The new shape function set will be termed `quad_1_gauss_shape_set`, and is introduced with a forward declaration:

\snippet custom_gaussian_element.cpp Forward declaration

Before defining the shape functions, we first define some basic properties of the shape function set by specialising the traits class template NiHu::shape_set_traits to the new shape set as follows:

\snippet custom_gaussian_element.cpp Shape traits

- We have defined the shape set's domain as the NiHu::quad_domain.
- The number of shape set nodes is 4.
- The polynomial order of the shape functions (the highest power of \f$\xi\f$ or \f$\eta\f$ in the definition of \f$L_i(\xi,\eta)\f$) is 1.
- The polynomial order of the Jacobian (the highest power of \f$\xi\f$ or \f$\eta\f$ in the product of the derivatives \f$L'_{\xi}\cdot L'_{\eta}\f$) is 1 too.

We mention here that the order of the Jacobian is only needed if the shape function set is used as a geometrical interpolation function set.
This will not be the case now, but we keep our code consistent with other shape function definitions.

After having the shape set traits defined, we can define the shape function class itself.
The new shape function class must be derived from the CRTP base NiHu::shape_set_base.
This base class defines a general interface for all shape sets, including convenient type definitions of the variable vector \f$\xi\f$, shape function vector \f$L_i(\xi)\f$ and its gradient matrix \f$\nabla_{\bf \xi} L_i(\xi)\f$.

The new derived class must define three static member functions.
- `eval_shape` evaluates the shape functions.
- `eval_dshape` evaluates the gradient of the shape functions.
- `corner_begin` returns a pointer to the first corner of the shape set.

\snippet custom_gaussian_element.cpp Shape class

The functions returning the shape function set and its derivatives are defined as

\snippet custom_gaussian_element.cpp Shape lsets

The shape function's nodal locations are stored in the static array `m_corners`. The `corner_begin` function returns the address of the array.

\snippet custom_gaussian_element.cpp Shape corners

That's all, we have defined the shape function set.
From now on, it can be used for geometrical interpolation or field interpolation purposes.
Furthermore, when using this shape function set in the collocational BEM context, the shape function nodes defined above are automatically used to generate weakly singular quadratures around the collocation points.


The field {#tut_custom_element_field}
=========

The new field is going to be based on the standard NiHu::quad_1_elem element, extended with our new shape function set.
The field will be termed `quad_1_gauss_field`, and is defined using a simple type definition:

\snippet custom_gaussian_element.cpp Field typedef

Each field type is automatically assigned an integer identifier.
However, we can override the field id definition by specialising the template structure ([metafunction]) NiHu::field_id to our new field type

\snippet custom_gaussian_element.cpp Field id

The defined field id will be used in the function space definition matrix to distinguish between different kind of fields.

\snippet custom_gaussian_element.cpp Field tag

Our new field type is ready to use in collocational, Galerkin or general BEM methods.

Example {#tut_custom_element_example}
=======

We present a simple example, the same as in tutorial \ref tut_laplace_double_integral, that evaluates a double singular integral of the Laplacian kernel on a square domain.
The square domain will be meshed with our new fields.

\snippet custom_gaussian_element.cpp main

The results read as
~~~~~~~~
matrix sum:   1.89281
analytical:   1.8928
log10 error: -5.86068
~~~~~~~~
proving the correctness and applicability of our brand new elements.