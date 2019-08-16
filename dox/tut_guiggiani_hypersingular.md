Hypersingular integrals and Guiggiani's method {#tut_guiggiani_hypersingular}
==============================================

\page tut_guiggiani_hypersingular Hypersingular integrals

[RongMethod]:http://dx.doi.org/10.1016/j.enganabound.2013.10.014
[GuiggianiMethod]:http://dx.doi.org/doi:10.1115/1.2893766

[TOC]

Introduction {#tut_guiggiani_intro}
============

This tutorial explains how strongly singular and hypersingular integrals are handled in NiHu.

The collocational integral

\f$
\displaystyle
I = \int_{S} K({\bf x}_0, {\bf x}) N({\bf x}) dS_x,
\f$

is singular if the singular point \f$ {\bf x}_0 \f$ is located inside the element domain \f$ S \f$.
- If the kernel contains an \f$ O(1/r^2) \f$ type singularity in 3D or an \f$ O(1/r) \f$ type singularity in 2D, then the integral is strongly singular.
- If the kernel contains an \f$ O(1/r^3) \f$ type singularity in 3D or an \f$ O(1/r^2) \f$ type singularity in 2D, then the integral is hypersingular.

Guiggiani's method  {#tut_guiggiani_theory}
==================

[Guiggiani][GuiggianiMethod] presented a unified method in 1992 for the accurate numerical evaluation of collocational strongly and hypersingular integrals.
Recently, [Rong et al][RongMethod] published an improved version of the original generic algorithm.
NiHu implements the improved version.

The original formulation  {#tut_guiggiani_original}
------------------------

As usual, the integration is performed in intrinsic coordinates:

\f$
\displaystyle
I = \int_{\Sigma} K({\bf \xi}_0, {\bf \xi}) N({\bf \xi}) J({\bf \xi}) d\Sigma,
\f$

where \f$ \Sigma \f$ denotes the reference domain of \f$ S \f$, \f$ \xi \f$ denotes the location vector in the local coordinate system, and \f$ \xi_0 \f$ denotes the image of the singular point.
\f$ J(\xi) \f$ is the Jacobian of the coordinate transform.
A further polar coordinate transform is introduced in the local coordinate system around the singular point with the definition

\f$
\displaystyle \xi = \xi_0 + \rho (\cos \theta, \sin \theta)
\f$

leading to the double integral

\f$
\displaystyle
I = \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} K(\rho, \theta) N(\rho, \theta) J(\rho, \theta) \rho d\rho d \theta
= \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} F(\rho, \theta) d \rho d \theta
\f$

The integrand \f$ F(\rho, \theta)\f$ is approximated with its Laurent series around the origin of the polar coordinate system (the singular point):

\f$
\displaystyle
F(\rho, \theta) = \frac{F_{-2}(\theta)}{\rho^2} + \frac{F_{-1}(\theta)}{\rho} + O(1)
\f$

and the truncated expansion is subtracted and added to the integrand to yield

\f$
\displaystyle
I = \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} F(\rho, \theta) - \frac{F_{-2}(\theta)}{\rho^2} - \frac{F_{-1}(\theta)}{\rho} d\rho d \theta
+ \int_{0}^{2\pi} -\frac{F_{-2}(\theta)}{\bar{\rho}(\theta)} + F_{-1}(\theta) \ln |\bar{\rho}(\theta)| d \theta \\
\displaystyle
I = \int_{0}^{2\pi} \int_{0}^{\bar{\rho}(\theta)} O(1) d\rho d \theta
+ \int_{0}^{2\pi} -\frac{F_{-2}(\theta)}{\bar{\rho}(\theta)} + F_{-1}(\theta) \ln |\bar{\rho}(\theta)| d \theta
\f$

In the last expression both the surface and the line integrals are regular, and can be approximated with standard Gaussian quadrature rules.

The method is fully general, it is valid for any type of curved or distorted surface elements and hypersingular kernels.
its only limitation is that the element should be smooth around the singular point.
For corners and edges, the line integrals are extended with additional terms, but remain regular.


Rong's improvement  {#tut_guiggiani_rong}
------------------

Recently, Rong et al proposed an efficiency improvement to Guiggiani's original formulation.
They showed that the intrinsic reference domain \f$ \Sigma \f$ can be chosen such that the derivative

\f$ \displaystyle \left| \lim_{\rho \to 0} \frac{\partial {\bf r}(\rho, \theta)}{\partial \rho} \right| \f$

remains constant (independent of \f$ \theta \f$).
As this derivative obviously plays an important role in the Laurent coefficients \f$ F_{-2}(\theta), F_{-1}(\theta) \f$ of the hypersingular kernel, the improvement yields improved accuracy on highly distorted elements.
Furthermore, the new reference domain yields analytical closed form expressions for the angular line integrals.

Implementation  {#tut_guiggiani_implementation}
==============

Guiggiani's method with Rong's improvement is implemented in class NiHu::guiggiani in file guiggiani_1992.hpp .

Class NiHu::guiggiani is templated on the trial field (element type and shape function type \f$ N \f$) and the kernel type, as well as on the quadrature orders in radial and tangential directions.
The general algorithm can be specialised to a specific kernel type by providing the singularity's Laurent expansion in the class template NiHu::polar_laurent_coeffs.

Note that the template NiHu::polar_laurent_coeffs is not templated on a particular kernel but on the kernel's singularity type.
The reason of this differentiation is that different kernels tipically share singularity types.
For example, the hypersingular Laplace and Helmholtz kernels differ only in their regular parts.
Their singular behaviour, and therefore their Laurent expansion around the singular point are the same.

The class NiHu::polar_laurent_coeffs should implement a static function template called eval that evaluates the laurent coefficients in terms of
- the 2nd order series expansion of the location vector \f$ {\bf r} \approx {\bf r}_1 \rho + {\bf r}_2 \rho^2 \f$
- the 1st order series expansion of the element Jacobian vector \f$ {\bf J} \approx {\bf J}_0 + {\bf J}_1 \rho \f$
- the 1st order series expansion of the shape function \f$ N \approx N_0 + N_1 \rho \f$
- the unit normal vector \f$ {\bf n} \f$ at the singular point

These members are obtained from the general guiggiani object taken as parameter by function eval.


