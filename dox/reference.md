NiHu Component Library Reference {#reference}
================================

\page reference

[TOC]

Introduction {#reference_intro}
============

This page is a reference guide of the NiHu Component Library.

Function Space Representations {#reference_func_space}
==============================

Elements {#reference_func_space_elements}
--------

 linear        					     | quadratic
-------------------------------------|---------------
 ::line_1_elem 						 | ::line_2_elem 
 ::tria_1_elem 						 | ::tria_2_elem 
 ::quad_1_elem, ::parallelogram_elem | ::quad_2_elem, ::quad_28_elem 

Fields and FieldViews {#reference_func_space_fields}
---------------------

- ::field_option::constant constant function space view generation option
- ::field_option::isoparametric isoparametric function space view generation option
- ::constant_view constant function space view generation factory function
- ::isoparametric_view isoparametric function space view generation factory function

Kernel Evaluation {#reference_kernel}
=================

Laplace equation
----------------

- ::laplace_2d_SLP_kernel
- ::laplace_2d_DLP_kernel
- ::laplace_2d_DLPt_kernel
- ::laplace_2d_HSP_kernel

- ::laplace_3d_SLP_kernel
- ::laplace_3d_DLP_kernel
- ::laplace_3d_DLPt_kernel
- ::laplace_3d_HSP_kernel

Helmholtz equation
------------------

- ::helmholtz_2d_SLP_kernel
- ::helmholtz_2d_DLP_kernel

- ::helmholtz_3d_SLP_kernel
- ::helmholtz_3d_DLP_kernel
- ::helmholtz_3d_DLPt_kernel
- ::helmholtz_3d_HSP_kernel

Linear isotropic elastostatics
------------------------------

- ::elastostatics_2d_U_kernel
- ::elastostatics_2d_T_kernel

- ::elastostatics_3d_U_kernel
- ::elastostatics_3d_T_kernel

Numerical Integration {#reference_integration}
=====================

Specialised singular integrals
------------------------------

Singular integrals implemented in library/laplace_singular_integrals.hpp

- Singular collocational and Galerkin integrals of the normal derivative kernels over constant lines and triangles

- Collocational integral of the Laplace 2D SLP kernel over a constant line
- Collocational integral of the Laplace 2D HSP kernel over a constant line
- Galerkin integral of the Laplace 2D SLP kernel over a constant line
- Galerkin integral of the Laplace 2D SLP kernel over a linear line

- Collocational integral of the Laplace 3D SLP kernel over a constant triangle
- Collocational integral of the Laplace 3D HSP kernel over a constant triangle

Singular integrals implemented in library/helmholtz_singular_integrals.hpp

- Singular collocational and Galerkin integrals of the normal derivative kernels over constant lines and triangles

- Collocational integral of the Helmholtz 2D SLP kernel over a constant line
- Collocational integral of the Helmholtz 3D SLP kernel over a constant triangle
- Collocational integral of the Helmholtz 3D HSP kernel over a constant triangle

- Galerkin integral of the Helmholtz 2D SLP kernel over a constant line

Singular integrals implemented in library/elastostatics_singular_integrals.hpp

- Galerkin integral of 2D elastostatics U kernel over a constant line
- Galerkin integral of 2D elastostatics T kernel over a constant line

Singular integrals implemented in library/lenoir_salles_2012.hpp

- Galerkin integral of the Laplace 3D SLP kernel over a constant triangle

