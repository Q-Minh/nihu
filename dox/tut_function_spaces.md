Function Spaces {#tut_function_spaces}
===============

\page tut_function_spaces Function spaces

[mesh definitions]:\ref tut_mesh_building

[TOC]

Introduction {#tut_funspace_intro}
============

Discretised function spaces are meshes extended with shape functions that describe the interpolation scheme of a certain physical quantity over the  elements of the mesh. NiHu provides two methods for building function spaces:
- The first method is to define the function space element-by-element, similar to [mesh definitions].
This method is rarely used, only in academic cases.
- The second method is automatic homogeneous function space view generation from an existing mesh, using a function space generation option.

Function space views {#tut_funspace_view}
====================

Constant and isoparametric views {#tut_funcspace_const_iso_view}
--------------------------------

Any NiHu mesh can be transformed into a piecewise constant function space using the library function NiHu::constant_view
~~~~~~~~~~
auto const &const_fsp = constant_view(my_mesh);
~~~~~~~~~~
Constant view means that the physical quantity is considered to be constant over each element,
and its nodal location is at the element center.

Alternatively, meshes can be viewed as piecewise isoparametric function
spaces, using the function NiHu::isoparametric_view:
~~~~~~~~~~
auto const &isoparam_fsp = isoparametric_view(my_mesh);
~~~~~~~~~~
Isoparametric view means that the physical quantity is interpolated
using the element's geometrical interpolation functions,
so its nodal locations coincide with the geometrical nodes.

\note
Function space view's are light weight (zero weight) objects, they only provide additional type information.
Therefore, it is important to store function space views in references.
This ensures that the mesh is not copied into a new object.
The function space refers to the same mesh that behaves like, or _is viewed as_ a piecewise constant function space.


### Alternative syntax

function space view's may be generated using an alternative, more general syntax:
~~~~~~~~~~
auto const &const_fsp    = create_function_space_view(my_mesh, NiHu::field_option::constant());
auto const &isoparam_fsp = create_function_space_view(my_mesh, NiHu::field_option::isoparametric());
~~~~~~~~~~
where the NiHu::field_option::constant and NiHu::field_option::isoparametric are tag types.
This syntax is more suitable for generic programming, and developers may customise their own field generation options.


Dirac view {#tut_funcspace_dirac_view}
----------

Any function space and function space view (constant, isoparametic, etc.)
can be transformed into a Dirac-view using the NiHu::dirac function.

~~~~~~~~~~
auto const &dirac_fsp = dirac(original_fsp);
~~~~~~~~~~

A Dirac function space behaves as a set of Dirac delta basis functions
located at the original function space's nodes.
For example

~~~~~~~~~~
auto const &dirac_center = dirac(constant_view(my_mesh));
~~~~~~~~~~

generates Dirac delta functions in the element centers of `my_mesh`, while

~~~~~~~~~~
auto const &dirac_nodes = dirac(isoparametric_view(my_mesh));
~~~~~~~~~~

generates Dirac delta functions in the nodal vertices of `my_mesh`.
The Dirac-like behaviour becomes important when we define weighted double integrals
(weighted residual matrices) with Dirac delta function spaces.

\note Using the `auto` keywords becomes inevitable when working with function space views. For example, in the code fragment below:
~~~~~~~~~~
auto my_mesh = create_mesh(nodes, elements, _tria_1_tag(), _quad_1_tag());
auto const &f_sp = dirac(constant_view(my_mesh));
~~~~~~~~~~
the second `auto` resolves to
~~~~~~~~~~
dirac_view<function_space_view<mesh<tmp::vector<NiHu::tria_1_elem, NiHu::quad_1_elem> >, NiHu::field_option::constant> >
~~~~~~~~~~
which is cumbersome to read and interpret.

Vector spaces {#tut_funspace_vector}
=============

A mesh can be extended into a vector function space by passing additional arguments to the view constructors.
The expression
~~~~~~~~~~
auto const &const_fsp = constant_view(my_mesh, _3d());
~~~~~~~~~~
creates a function space that interpolates a vector valued quantity over a mesh with piecewise constant interpolation funcions.


Defining heterogeneous custom function spaces {#tut_funspace_custom}
=============================================

Heterogeneous function spaces are a collection of different element types with different field interpolation shape functions.
NiHu provides a way (function NiHu::create_function_space) to manually define heterogeneous function spaces in a similar manner as defining a mesh.
Although this option is rarely needed, we refer to the tutorial \ref tut_custom_element for further information.

