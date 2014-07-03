Introducing a new family of PDE into the NiHu framework {#tut_custom_kernel}
=======================================================

\page tut_custom_kernel

[metafunction]:\ref tech_meta_metafun

[TOC]

Introduction {#tut_custom_kernel_intro}
============

The purpose of this tutorial is to demonstrate how a new family of partial differential equations (PDE) can be introduced into the NiHu framework.
In a boundary element context, the introduction of a new PDE is equivalent to the implementation of its fundamental solutions.
The demonstrative example is linear isotropic elastostatics.

Theory {#tut_custom_kernel_theory}
======

PDE {#tut_custom_kernel_pde}
---

3D linear isotropic elasticity is governed by the PDE

\f$
\displaystyle
\sigma_{ij,j} = \delta_{ij} u_{k,ki} + \mu \left( u_{i,jj} + u_{j,ij} \right) \\
\displaystyle
\mu u_{i,jj}(x) + \left(\mu + \lambda\right) u_{j,ij}(x) = 0, \quad x \in \Omega \subset \mathbb{R}^{3}
\f$

with the boundary conditions

\f$
\displaystyle
t_i(x) = \bar{t}_i(x), x \in \Gamma_t, \\
u_i(x) = \bar{u}_i(x), x \in \Gamma_u
\f$

where \f$ u_i(x) \f$ denotes the displacement vector field, \f$ t_i(x) \f$ denotes the traction vector field, \f$ \mu \f$ and \f$ \lambda \f$ are the Lamé-coefficients, \f$ \Omega \f$ denotes the solution domain and \f$ \Gamma_{\alpha} \f$ stands for its boundary with prescribed traction or displacement boundary conditions.

BIE {#tut_custom_kernel_bie}
---

The equivalent boundary integral representation of the PDE is

\f$ \displaystyle
\int_{\Gamma} t^*_{ij}(x,y) u_j(y) \mathrm{d} y - \int_{\Gamma} u^*_{ij}(x,y) t_{j}(y) \mathrm{d} y = u_i(x)
\f$

where the displacement and traction fundamental solutions are

\f$ \displaystyle
u^*_{ij} = \frac{(3-4\nu) \delta_{ij} + r,_i r,_j}{16 \pi \mu (1-\nu) r},\\
\displaystyle
t^*_{ij} = \frac{-r,_n ((1-2\nu)\delta_{ij} + 3 r,_i r,_j) + (1-2\nu) (r,_i n_j - r,_j n_i)}{8 \pi (1-\nu) r^2}
\f$

where \f$ \nu \f$ and \f$ \mu \f$ denote the Poisson's number and the shear modulus, respectively.
The displacement kernel contains a weak \f$ O(1/r) \f$ singularity, and can be integrated with blind singular quadratures.
The traction kernel is strongly singular, and its \f$ O(1/r^2) \f$-type singularity needs to be handled in a CPV sense.

Implementing the fundamental solution {#tut_custom_kernel_fundsol}
=====================================

In order to be able to define expressions for the fundamental solutions, we need to define three classes:
- The kernel's test input \f$ x \f$
- The kernel's trial input \f$ y \f$
- The kernel's parameters

For the case of the displacement kernel, both inputs are simple locations.
For the case of the traction kernel, the test input is location, while the trial input is location with normal vector.
The kernel parameter class needs to encapsulate the two material properties:
- The Poisson's ratio \f$ \nu \f$
- The Shear modulus \f$ \mu \f$

As the shear modulus serves only as a scaling factor, it can be omitted in the implementation without the loss of generality.

The parameter class {#tut_custom_kernel_param}
-------------------

First the kernel parameter class is defined that encapsulates a simple `double` as the Poisson's ratio:

\snippet custom_kernel.mex.cpp parameter

The kernel functor {#tut_custom_kernel_ukernelfunctor}
------------------

As the next step, a function object (functor) named `Ukernel` is defined that returns the value of the displacement fundamental solution.

\snippet custom_kernel.mex.cpp ufunctor

As it is clear from the definition, the functor not only returns the kernel's expression but defines the kernel's return type (3d double matrix) as well.

The kernel class {#tut_custom_kernel_ukernelclass}
----------------

The next step is the definition of the final kernel class that is going to be used in numerical integrations.
This is done in three steps.
	1. The kernel class is declared
	2. Compile time properties of the kernel class are defined in its traits class
	3. The kernel class is defined, derived from ::kernel_base using the CRTP pattern.

Step (1): The kernel class is forward declared as

\snippet custom_kernel.mex.cpp udeclare

Step(2): The kernel traits are defined by specialising the traits class ::kernel_traits as follows

\snippet custom_kernel.mex.cpp utraits

The following traits are defined:
- members  `test_input_t` and `trial_input_t` select the kernel's test and trial inputs, as mentioned before.
- member `data_t` denotes the kernel's parameter class, and member `output_t` denotes the kernel's output class.
NiHu is capable of generating combined kernels optimised for parallel evaluations, if the kernel's inputs, parameters and outputs are defined in terms of building blocks (bricks). This feature is not exploited in the present example, meaning that the output type is a single brick wall, defined by the functor class `Ukernel`.
- member `result_dimensions` indicates that the kernel is tensor valued (vector PDE)
- member `far_field_behaviour` indicates that the kernel's far field behaviour is \f$ O(1/r) \f$.
- member `quadrature_family_t` indicates that the kernel is integrated using Gaussian quadratures in the far field.
- members `is_symmetric` and `is_singular` indicate that the kernel is symmetric and singular.

As the kernel is defined as singular, additional kernel traits need to be defined by specialising template ::singular_kernel_traits

\snippet custom_kernel.mex.cpp usingtraits

The singular traits define the singularity type as \f$ O(1/r) \f$.
In 3D, such a singularity is weak, and can be cancelled out by blind quadrature methods.
Member `singular_quadrature_order` defines that the singular quadratures are going to be evaluated with a 7-th order blind singular quadrature rule.

Step (3): Finally, the kernel class can be defined:

\snippet custom_kernel.mex.cpp udefine

The kernel is derived from ::kernel_base, and a constructor is provided that constructs the kernel object from a single Poisson's ratio value.

The traction kernel {#tut_custom_kernel_tractionkernel}
-------------------

The traction kernel is defined similarly. Without detailed explanation, we simply print the corresponding code snippets as follows:

\snippet custom_kernel.mex.cpp tkernel

Compared to the displacement kernel, the only significant difference is the singular behaviour that is defined as \f$ O(1/r^2) \f$.
As this singularity is strong, it can not be handled by blind quadratures.
Instead, we specialise Guiggiani's method to integrate the strong singularity.

