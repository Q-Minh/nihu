Hypersingular integrals and Guiggiani's method {#tut_guiggiani_hypersingular}
==============================================

\page tut_guiggiani_hypersingular

[Eigen]:http://eigen.tuxfamily.org/index.php?title=Main_Page

[TOC]

Introduction {#tut_guiggiani_intro}
============

This tutorial explains how hypersingular integrals are handled in NiHu.

The integral

\f$
\displaystyle
I = \int_{S} K({\bf x}_0, {\bf x}) N({\bf x}) dS_x,
\f$

is hypersingular when the kernel contains an \f$ O(1/r^3) \f$ type singularity in 3D or an \f$ O(1/r^2) \f$ type singularity in 2D, and the singular point \f$ {\bf x}_0 \f$ is located inside the element domain \f$ S \f$.

Guiggiani's method  {#tut_guiggiani_theory}
==================

The original formulation  {#tut_guiggiani_original}
------------------------

As usual, the integration is performed in intrinsic coordinates:

\f$
\displaystyle
I = \int_{\Sigma} K({\bf \xi}_0, {\bf \xi}) N({\bf \xi}) J({\bf \xi}) d\Sigma,
\f$

A polar coordinate transform is introduced around the singular point in the reference domain with the definition

\f$
\displaystyle \xi = \xi_0 + \rho (\cos \theta, \sin \theta)
\f$

leading to the integral

\f$
\displaystyle
I = \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} K(\rho, \theta) N(\rho, \theta) J(\rho, \theta) \rho d\rho d \theta
= \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} F(\rho, \theta) d \rho d \theta
\f$

The integrand \f$ F(\rho, \theta)\f$ is approximated with its Laurent series around the orgin of the polar coordinate system:

\f$
\displaystyle
F(\rho, \theta) = \frac{F_{-2}(\theta)}{\rho^2} + \frac{F_{-1}(\theta)}{\rho} + O(1)
\f$

and is subtracted and added to the integrand to yield

\f$
\displaystyle
I = \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} F(\rho, \theta) - \frac{F_{-2}(\theta)}{\rho^2} - \frac{F_{-1}(\theta)}{\rho} d\rho d \theta
+ \int_{0}^{2\pi} -\frac{F_{-2}(\theta)}{\bar{\rho}(\theta)} + F_{-1}(\theta) \ln |\bar{\rho}(\theta)| d \theta
\f$

In the last expression both surface and the line integrals are regular, and can be approximated with standard Gaussian quadrature rules.

The method is fully general, it is valid for any type of curved or distorted surface elements and hypersingular kernels. The only limitation of the above formula is that the element should be smooth around the singular point.
For corners and edges, the line integrals are extended with additional terms, but remain regular.


Rong's improvement  {#tut_guiggiani_rong}
------------------

Recently, Rong et al proposed a serious efficiency improvement to Guiggiani's original formulation.
They showed that the intrinsic reference domain \f$ \Sigma \f$ can be chosen such that the derivative

\f$ \displaystyle \left| \lim_{\rho \to 0} \frac{\partial {\bf r}(\rho, \theta)}{\partial \rho} \right| \f$

remains constant over the integration domain.
As this derivative obviously plays an important role in the Laurent coefficients \f$ F_{-2}(\theta), F_{-1}(\theta) \f$ of the hypersingular kernel, the improvement allows for improved accuracy on highly distorted elements.
Furthermore, the new reference domain yields analytical closed form expressions for the single angular integrals.

Implementataion  {#tut_guiggiani_implementation}
===============

Guiggiani's method with Rong's improvement is implemented in class ::guiggiani in file guiggiani_1992.hpp
