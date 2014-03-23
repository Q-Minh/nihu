Introducing a new family of PDE into the NiHu framework {#tut_custom_kernel}
=======================================================

\page tut_custom_kernel

[metafunction]:\ref tech_meta_metafun

[TOC]

Introduction {#tut_custom_kernel_intro}
============

The purpose of this tutorial is to demonstrate how a new family of partial differential equations (PDE) can be introuced into the NiHu framework.
In a boundary element context, the introduction of a new PDE is equaivalent with the implementation of its fundamental solutions.
The demonstrative example is linear isotropic elasticity.

Theory {#tut_custom_kernel_theory}
======

PDE {#tut_custom_kernel_pde}
---

3D linear isotropic elasticity is governed by the PDE

\f$
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

The equaivalent boundary integral representation of the PDE is

\f$ \displaystyle
\int_{\Gamma} t^*_{ij}(x,y) u_j(y) \mathrm{d} y - \int_{\Gamma} u^*_{ij}(x,y) t_{j}(y) \mathrm{d} y = c(x) u_i(x)
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



