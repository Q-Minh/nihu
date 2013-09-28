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

 linear        | quadratic
---------------|---------------
 ::line_1_elem | ::line_2_elem 
 ::tria_1_elem | ::tria_2_elem 
 ::quad_1_elem | ::quad_2_elem 

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

