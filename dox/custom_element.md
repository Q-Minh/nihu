Custom element {#custom_element}
==============

[TOC]


Introduction {#intro}
============

In this tutorial we introduce a new element type into the NiHu toolbox.
The new element type will be a 3D four noded quadrangle element special in the sense that its nodes are located inside the element, at the Gaussian integration points of the element domain.

The nodal locations are

\f$\xi_1 = \left\{-\sqrt{3}/3, -\sqrt{3}/3\right\} \\
\xi_2 = \left\{+\sqrt{3}/3, -\sqrt{3}/3\right\} \\
\xi_3 = \left\{+\sqrt{3}/3, +\sqrt{3}/3\right\} \\
\xi_4 = \left\{-\sqrt{3}/3, +\sqrt{3}/3\right\}\f$


The shape functions are defined as

\f$L_1 = (1-\sqrt{3}\xi)(1-\sqrt{3}\eta)/4 \\
L_2 = (1+\sqrt{3}\xi)(1-\sqrt{3}\eta)/4 \\
L_3 = (1+\sqrt{3}\xi)(1+\sqrt{3}\eta)/4 \\
L_4 = (1-\sqrt{3}\xi)(1+\sqrt{3}\eta)/4\f$

and their derivatives are

\f$L'_{1,\xi} = (-\sqrt{3})(1-\sqrt{3}\eta)/4, \quad L'_{1,\eta} = (1-\sqrt{3}\xi)(-\sqrt{3})/4 \\
L'_{2,\xi} = (+\sqrt{3})(1-\sqrt{3}\eta)/4, \quad L'_{1,\eta} = (1-\sqrt{3}\xi)(-\sqrt{3})/4 \\
L'_{3,\xi} = (+\sqrt{3})(1+\sqrt{3}\eta)/4, \quad L'_{1,\eta} = (1-\sqrt{3}\xi)(-\sqrt{3})/4 \\
L'_{4,\xi} = (-\sqrt{3})(1+\sqrt{3}\eta)/4, \quad L'_{1,\eta} = (1-\sqrt{3}\xi)(-\sqrt{3})/4\f$



The shape set {#shapeset}
-------------

In order to define the new element, we need to define its interpolation functions first. The new shape function set will be termed quad_1_gauss_shape_set, and is introduced with a forward declaration:

\snippet custom_gaussian_element.hpp Forward declaration

Before defining the shape functions, we first define some basic properties of the shape function set by specialising the traits class template shape_set_traits to the new shape set as follows:

\snippet custom_gaussian_element.hpp Shape traits

We have defined the shape set domain as the quad_domain. The number of shape set nodes is 4, the polynomial order of the shape functions (the highest power of \f$\xi\f$ or \f$\eta\f$ in their definition) is 1, and the polynomial order of the Jacobian is one too.

After having defined the shape set traits, we can define the shape function class itself.
The class must define three static member functions.
eval_shape evaluates the shape functions.
eval_dshape evaluate the gradient of the shape functions.
corner_begin returns a pointer to the first corner of the shape set.

\snippet custom_gaussian_element.hpp Shape class

The functions returning the shape function set and its derivatives are defined as

\snippet custom_gaussian_element.hpp Shape lsets

The shape function's nodal locations are stored in the static array m_corners. The corner_begin function returns the address of the array.

\snippet custom_gaussian_element.hpp Shape corners



