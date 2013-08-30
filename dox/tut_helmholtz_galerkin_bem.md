A Galerkin BEM for the Helmholtz equation {#tut_helmholtz_galerkin_bem}
=========================================

\page tut_helmholtz_galerkin_bem

Introduction {#tut_helmholtz_galerkin_bem_intro}
============

This tutorial demonstrates the usage of NiHu for solving the Helmholtz equation in 3D by a Galerkin type boundary element method (BEM).
The Helmholtz equation is the governing equation of linear time-harmonic wave propagation in a homogeneous medium.
In this tutorial the Helmholtz equation will be applied for the solution of an acoustic wave propagation problem.
In this tutorial we will use NiHu's Matlab interface for the mesh definition and the solution of the resulting matrix equations.

Theory {#tut_helmholtz_galerkin_bem_theory}
======

The Helmholtz integrals {#tut_helmholtz_galerkin_bem_integrals}
-----------------------

The inhomogeneous Helmholtz equation with a Dirac-delta excitation located at the point \f$ {\bf y} \f$ is given in three dimensions as

\f$ \nabla^2 p({\bf x}) + k^2 p({\bf x}) = - \delta({\bf y}) \f$

The sound pressure radiated by the closed vibrating surface \f$ S \f$ is obtained by means of the Helmholtz integrals

\f$ \displaystyle \frac{1}{2} p({\bf x}) = \int_S G({\bf x}, {\bf y}) q({\bf y}) dS_y - \int_S H({\bf x}, {\bf y}) p({\bf y}) dS_y \qquad \text{if} \quad x \in S \f$

\f$ p({\bf x}) = \displaystyle \int_S G({\bf x}, {\bf y}) q({\bf y}) dS_y - \int_S H({\bf x}, {\bf y}) p({\bf y}) dS_y \qquad \text{if} \quad x \in F \f$

where \f$ G({\bf x}, {\bf y}) \f$ and \f$ H({\bf x}, {\bf y}) \f$ denote the fundamental solution of the inhomogeneous Helmholtz equation and its normal derivative with respect to the source point.
The functions are given as

\f$ G({\bf x}, {\bf y}) = \exp(-ikr)/4\pi r \f$

\f$ H({\bf x}) = \displaystyle \frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial n_y} = -G({\bf x}, {\bf y}) \left(\displaystyle \frac{1 + 1 \mathrm{i} k r}{r} \frac{{\bf r} \cdot {\bf n}_y }{r} \right)  \f$

with the distance vector \f$ \bf r \f$ and the distance \f$ r \f$ defined as

\f$ \mathbf{r} = \mathbf{y} - \mathbf{x} \qquad \text{and} \qquad  r = \left| \mathbf{r} \right| . \f$

The normal vector \f$ \mathbf{n}_y \f$ is pointing outward from the surface \f$ S \f$ at the point \f$ \mathbf{y} \f$.

Spatial discretization {#tut_helmholtz_galerkin_bem_discretization}
----------------------

In order to solve the Helmholtz integral equations numerically, the radiating surface \f$ S \f$ is divided into a number of non-overlapping elements.
The sound pressure function \f$ p(\mathbf{x}) \f$ and its normal derivative \f$ q(\mathbf{x}) \f$ are approximated by means of weighting functions \f$ w(\mathbf{x}) \f$ such that

\f$ p(\mathbf{x}) = \displaystyle \sum_{j=1}^n w_j (\mathbf{x}) p_j \qquad \text{and} \qquad q(\mathbf{x}) = \displaystyle \sum_{j=1}^n w_j (\mathbf{x}) q_j , \f$

with \f$ p_j \f$ and \f$ q_j \f$ denoting the corresponding weights and \f$ n \f$ representing the total number of weighting functions, also known as the number of degrees of freedom.

As a next step, the weak form of the Helmholtz integral equation is obtained by multiplying the equations with a set of test functions \f$ t_j(\mathbf{x}) \f$ and integrating the result over the complete surface \f$ S \f$ with respect to the variable \f$ \mathbf{x} \f$.
In the Galerkin formulation, the test function are equivalent to the weighting functions, i.e. \f$ t_j(\mathbf{x}) \equiv w_j(\mathbf{x}) \f$.

This leads to the equation 

\f$ \left< w_i , \left(\left[ \mathcal{M} - \frac{1}{2} \mathcal{I} \right] w_j \right)_S \right>_S p^{(s)}_j = \left< w_i, \left(\mathcal{L} w_j \right)_S \right>_S q^{(s)}_j \f$

For the calculation of the acoustic field of the radiated surface a Dirac-delta test function is applied herein, which gives

\f$ p^{(f)}_i = \left< \delta_i , \left( \mathcal{M} w_j \right)_S \right>_F p^{(s)}_j - \left< \delta_i, \left(\mathcal{L} w_j \right)_S \right>_F q^{(s)}_j \f$

The equations can be rewritten by utilizing matrix notations 

\f$ M^{(s)}_{ij} = \left< w_i , \left(\left[ \mathcal{M} - \frac{1}{2} \mathcal{I} \right] w_j \right)_S \right>_S \qquad \f$
and
\f$ \qquad L^{(s)}_{ij} = \left< w_i, \left(\mathcal{L} w_j \right)_S \right>_S \f$

\f$ M^{(f)}_{ij} = \left< \delta_i , \left( \mathcal{M}  w_j \right)_S \right>_F \qquad \f$
and
\f$ \qquad L^{(f)}_{ij} = \left< \delta_i , \left( \mathcal{L}  w_j \right)_S \right>_F \f$

in the forms

\f$ \mathbf{M}^{(s)} \mathbf{p}^{(s)} = \mathbf{L}^{(s)} \mathbf{q}^{(s)} \f$

\f$ \mathbf{p}^{(f)} = \mathbf{L}^{(f)} \mathbf{q}^{(s)} - \mathbf{M}^{(f)} \mathbf{p}^{(s)} \f$

The vectors \f$ \mathbf{p} \f$ and \f$ \mathbf{q} \f$ denote the column vector of the weights for the degrees of freedom of the pressure and its normal derivative, respectively.
The upper inidices \f$ \cdot^{(s)} \f$ and \f$ \cdot^{(f)} \f$ symbolise the variables on the surface \f$ S \f$ and the field \f$ F \f$ meshes, respectively. 