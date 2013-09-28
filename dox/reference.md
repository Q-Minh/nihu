NiHu Component Library Reference {#reference}
================================

\page reference

[TOC]

Introduction {#reference_intro}
============

This pge is a reference guide of the NiHu Component Library.

Function Space Representations {#reference_func_space}
==============================

Elements {#reference_func_space_elements}
--------

- ::line_1_elem linear 2-noded line element in 2D
- ::line_2_elem quadratic 3-noded line element in 2D
- ::tria_1_elem linear 3-noded triangle element in 3D
- ::tria_2_elem quadratic 6-noded triangle element in 3D
- ::quad_1_elem linear 4-noded quadrangle element in 3D
- ::quad_2_elem quadratic 9-noded quadrangle element in 3D

Fields and FieldViews {#reference_func_space_fields}
---------------------

- ::field_option::constant constant function space view generation option
- ::field_option::isoparametric isoparametric function space view generation option
- ::constant_view constant function space view generation factory function
- ::isoparametric_view isoparametric function space view generation factory function

Kernel Evaluation {#reference_kernel}
=================

- ::laplace_2d_SLP_kernel
- ::laplace_2d_DLP_kernel
- ::laplace_2d_DLPt_kernel
- ::laplace_2d_HSP_kernel

- ::laplace_3d_SLP_kernel
- ::laplace_3d_DLP_kernel
- ::laplace_3d_DLPt_kernel
- ::laplace_3d_HSP_kernel

- ::helmholtz_3d_SLP_kernel
- ::helmholtz_3d_DLP_kernel
- ::helmholtz_3d_DLPt_kernel
- ::helmholtz_3d_HSP_kernel

Numerical Integration {#reference_integration}
=====================

Regularisation
--------------

- Lenoir and Salles 2012
- Matsumoto 2010

